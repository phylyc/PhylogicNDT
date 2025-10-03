#!/bin/bash -l
#SBATCH --job-name=LeagueModel      # SLURM job name
#SBATCH --cpus-per-task=1           # 1 CPU
#SBATCH --mem=8G                    # 8 GB RAM
#SBATCH --time=04:00:00             # 4h runtime
#SBATCH -p normal                   # partition
#SBATCH --output=logs/%x_%j.log     # log file
#SBATCH --export=NONE               # clean env

set -euo pipefail                   # strict error handling

# pinned python env (condo python unreliable)
#ENV="/data/brastianoslab/dgritsch/envs/phylogic"
#PY="$ENV/bin/python"
#
#export PYTHONNOUSERSITE=1; unset PYTHONPATH
#export PATH="$ENV/bin:$PATH"
#export LD_LIBRARY_PATH="$ENV/lib:${LD_LIBRARY_PATH-}"
#export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1

mkdir -p logs

#cd "${SLURM_SUBMIT_DIR:-$PWD}"

# args: --cohort <name> --comps <file1,file2,...>
usage(){ echo "Usage: $0 --cohort <cohort> --comps <files>"; exit 1; }

COHORT=""; COMPS_RAW=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohort) COHORT="$2"; shift 2;;
    --comps)  COMPS_RAW="$2"; shift 2;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

[[ -z "$COHORT" || -z "$COMPS_RAW" ]] && usage

IFS=',' read -r -a COMPS_ARRAY <<< "$COMPS_RAW"

# validate inputs
abspath() { [[ "$1" = /* ]] && echo "$1" || echo "$PWD/$1"; }

declare -a INFILES=()

for f in "${COMPS_ARRAY[@]}"; do
  af="$(abspath "$f")"
  if [ ! -f $af ] ; then
    echo "Missing file: $af"
    continue
  fi
  [[ -s "$af" ]] || { echo "Missing file: $af"; exit 1; }
  INFILES+=("$af")
done

#OUTDIR="league_comps_clean_${SLURM_JOB_ID:-$$}"
OUTDIR="league_comps_clean"

mkdir -p "$OUTDIR"

# clean comps tsv: normalize headers, drop bad/self rows, renorm probs
clean_one() {
  local in="$1"
  local out="$2"

  # Choose reader based on input suffix
  local -a reader=(cat -- "$in")
  [[ "$in" == *.gz ]] && reader=(gzip -cd -- "$in")

  # Choose writer based on output suffix
  local -a writer=(cat)
  [[ "$out" == *.gz ]] && writer=(gzip -c)

  "${reader[@]}" \
  | awk -F'\t' -v OFS='\t' -v infile="$in" -v outfile="$out" '
    function tolow(s){gsub(/\r/,"",s); return tolower(s)}
    function trim(s){gsub(/^[ \t\r]+|[ \t\r]+$/,"",s); return s}
    BEGIN{ samp_i=event1_i=event2_i=p12_i=p21_i=pu_i=0; kept=0; dropped=0 }
    NR==1{
      for(i=1;i<=NF;i++){
        h=tolow($i)
        if(h=="samp"||h=="sample_id") samp_i=i
        else if(h=="eve1"||h=="event1") event1_i=i
        else if(h=="eve2"||h=="event2") event2_i=i
        else if(h=="p1->2"||h=="prob_event1_before_event2") p12_i=i
        else if(h=="p2->1"||h=="prob_event2_before_event1") p21_i=i
        else if(h=="unknown"||h=="prob_unknown") pu_i=i
      }
      if(!(samp_i&&event1_i&&event2_i&&p12_i&&p21_i&&pu_i)){
        print "bad header" > "/dev/stderr"; exit 2
      }
      print "samp","eve1","eve2","p1->2","p2->1","unknown"; next
    }
    {
      s=trim($samp_i); e1=trim($event1_i); e2=trim($event2_i)
      p12=$p12_i; p21=$p21_i; pu=$pu_i
      te1=tolow(e1); te2=tolow(e2)
      if(e1==""||e2==""||te1~/^(na|nan|\.|none)$/||te2~/^(na|nan|\.|none)$/){dropped++; next}
      if(p12==""||p21==""||pu==""||p12+0!=p12||p21+0!=p21||pu+0!=pu){dropped++; next}
      if(p12<0||p21<0||pu<0||p12>1||p21>1||pu>1){dropped++; next}
      sum=p12+p21+pu; if(sum<=0){dropped++; next}
      if(sum<0.999||sum>1.001){p12/=sum;p21/=sum;pu/=sum}
      if(e1==e2){dropped++; next}
      print s,e1,e2,p12,p21,pu; kept++
    }
    END{
      printf("[CLEAN] %s -> %s kept=%d dropped=%d\n", infile, outfile, kept, dropped) > "/dev/stderr"
    }
  ' \
  | "${writer[@]}" > "$out"
}

declare -a CLEANED=()

for f in "${INFILES[@]}"; do
  base="$(basename "$f")"
  out="$OUTDIR/${base/.tsv.gz/.clean.tsv.gz}"
  clean_one "$f" "$out"
  [[ -s "$out" ]] || { echo "empty cleaned file: $out"; exit 1; }
  CLEANED+=("$out")
done

python /carterlab/phahnel/tools/PhylogicNDT/PhylogicNDT.py \
  LeagueModel \
  --cohort "$COHORT" \
  --comps "${CLEANED[@]}" \
  --n_perms 500 \
  --num_games_against_each_opponent 3 \
  --max_num_snvs 100 \
  --max_num_cnv_arms 50 \
  --max_num_cnv_arm_gains 20 \
  --max_num_cnv_arm_losses 20 \
  --max_num_cnv_focal 50 \
  --max_num_homdel 10 \
  --min_event_prevalence 0.05

echo "LeagueModel finished for cohort: $COHORT "
