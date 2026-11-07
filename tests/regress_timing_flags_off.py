"""Check that the Timing tool gives identical tables with all new flags off, on master and on this branch.

Usage, from a repository root: python tests/regress_timing_flags_off.py <output_dir> <patient_dir> [<patient_dir> ...]
Run it once per checkout, with PYTHONHASHSEED=0 set, then compare the two output folders. Patients come from an
existing somatic_workflow output folder, with a REP0 timing SIF and mut_ccfs inside.
"""
import argparse
import glob
import os
import random
import sys

import numpy as np
import pandas as pd


def main():
    output_dir, patient_dirs = os.path.abspath(sys.argv[1]), sys.argv[2:]
    repo = os.path.abspath('.')
    sys.path.insert(0, repo)
    from SinglePatientTiming import SinglePatientTiming as spt
    for patient_dir in patient_dirs:
        patient_dir = os.path.abspath(patient_dir)
        patient = os.path.basename(patient_dir)
        out = os.path.join(output_dir, patient)
        os.makedirs(out, exist_ok=True)
        sif = pd.read_csv(os.path.join(patient_dir, f'{patient}.REP0.sif.timing.txt.gz'), sep='\t')
        sif['maf_fn'] = [os.path.join(patient_dir, os.path.basename(p)) for p in sif['maf_fn']]
        sif['seg_fn'] = [os.path.join(patient_dir, os.path.basename(p)) for p in sif['seg_fn']]
        sif_path = os.path.join(out, 'input.sif')
        sif.to_csv(sif_path, sep='\t', index=False)
        args = argparse.Namespace(
            artifact_blacklist=None, artifact_whitelist=None, indiv_id=patient, driver_genes_file=None, genome_build='hg19',
            sif=sif_path, sample_data=None, grid_size=101, intersect_cn_trees=False, gistic_fn=os.environ.get('FOCAL_REGIONS'),
            min_supporting_muts=3, run_comparisons=False, sex_aware=False, continuous_cn=False, clonal_unresolved=False,
            segment_homdels=False, homdel_artifact_file=None, homdel_max_length=5000000, homdel_min_posterior=0.95,
            cn_clonal_min_ccf=0.8)
        random.seed(7)
        np.random.seed(7)
        os.chdir(out)
        spt.run_tool(args)
        events = spt.read_timed_events(f'{patient}.timing.tsv', compare_all=True)[:300]
        spt.PhylogicOutput.write_comp_table(patient, spt.compare_event_list(events))
        os.chdir(repo)
        print(patient, 'done', flush=True)


if __name__ == '__main__':
    main()
