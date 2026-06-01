## imports ##
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import csv
from pathlib import Path
import re
import numpy as np
import seaborn as sns
from matplotlib import gridspec
import matplotlib.pyplot as plt
import operator
from scipy.stats import gaussian_kde
import itertools
import copy
import logging

_SINGLE_RESIDUE_EVENT_RE = re.compile(r'^p\.([A-Z\*])(\d+)([A-Z\*])$')
_INT_TOKEN_RE = re.compile(r'(\d+)')
_CODING_SILENT_CLASS = 'SILENT'
_CODING_NONSYN_CLASS = 'NONSYN'

# Canonical column names expected by load_df() and the aliases accepted
# from patient-level comparison tables.
_COLUMN_ALIASES = {
    'samp': 'sample',
    'sample_id': 'sample',
    'sample': 'sample',
    'eve1': 'event1',
    'event1': 'event1',
    'eve2': 'event2',
    'event2': 'event2',
    'p1->2': 'p_event1_win',
    'prob_event1_before_event2': 'p_event1_win',
    'p_event1_win': 'p_event1_win',
    'p2->1': 'p_event2_win',
    'prob_event2_before_event1': 'p_event2_win',
    'p_event2_win': 'p_event2_win',
    'unknown': 'unknown',
    'prob_unknown': 'unknown',
}


def _normalize_comp_columns(df):
    """Rename comparison-table columns to the canonical names expected by
    ``League.load_df``.  Unrecognised columns are left unchanged."""
    rename = {}
    for col in df.columns:
        mapped = _COLUMN_ALIASES.get(col.strip().lower())
        if mapped is not None and mapped not in rename.values():
            rename[col] = mapped
    if rename:
        df = df.rename(columns=rename)
    return df


def _normalize_gene_token(value):
    return str(value).strip().upper()


def _load_hotspot_index(hotspot_fn=None):
    single_residue = set()
    indel_ranges = {}
    if hotspot_fn is None:
        return single_residue, indel_ranges
    try:
        with open(hotspot_fn, 'r', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for row in reader:
                gene = _normalize_gene_token(row.get('Gene', ''))
                residue = str(row.get('Residue', '')).strip()
                htype = str(row.get('Type', '')).strip().lower()
                if not gene or not residue:
                    continue
                if htype == 'single residue':
                    m = re.match(r'^([A-Z\*])(\d+)$', residue)
                    if not m:
                        continue
                    single_residue.add((gene, int(m.group(2))))
                elif htype == 'in-frame indel':
                    if '-' in residue:
                        a, b = residue.split('-', 1)
                        try:
                            start, end = int(a), int(b)
                        except ValueError:
                            continue
                    else:
                        try:
                            start = end = int(residue)
                        except ValueError:
                            continue
                    indel_ranges.setdefault(gene, []).append((min(start, end), max(start, end)))
    except Exception as exc:
        logging.warning('Unable to load hotspot file %s: %s', hotspot_fn, exc)
    return single_residue, indel_ranges


def _legacy_event_gene(event):
    event = str(event)
    if event.split('_')[0] in ['loss', 'gain', 'homdel']:
        return '_'.join(event.split('_')[0:2])
    if event.split('_')[0] in ['Focal']:
        return '_'.join(event.split('_')[1:])
    return event.split('_')[0].split(':')[0]


def _parse_point_event(event):
    event = str(event).strip()
    if '_' in event:
        gene, mutation = event.split('_', 1)
        gene = gene.strip()
        mutation = mutation.strip()
        return gene if gene else None, mutation
    return None, event


def _classify_point_mutation(event):
    gene, mutation = _parse_point_event(event)
    if not mutation:
        return {
            'gene': gene,
            'mutation_token': mutation,
            'is_silent': False,
            'is_nonsyn': True,
            'position': None,
            'point_event_with_gene': bool(gene),
        }
    m = _SINGLE_RESIDUE_EVENT_RE.match(mutation)
    if m:
        ref_aa, pos, alt_aa = m.group(1), int(m.group(2)), m.group(3)
        is_silent = (ref_aa == alt_aa and ref_aa != '*' and alt_aa != '*')
        return {
            'gene': gene,
            'mutation_token': mutation,
            'is_silent': bool(is_silent),
            'is_nonsyn': not bool(is_silent),
            'position': pos,
            'point_event_with_gene': bool(gene),
        }
    ints = [int(x) for x in _INT_TOKEN_RE.findall(mutation)]
    return {
        'gene': gene,
        'mutation_token': mutation,
        'is_silent': False,
        'is_nonsyn': True,
        'position': ints[0] if ints else None,
        'point_event_with_gene': bool(gene),
    }

class Eve_Pair():

    def __init__(self, event1, event2, event1_type, event2_type):
        self.event1 = event1
        self.event_type1 = event1_type
        self.event2 = event2
        self.event_type2 = event2_type
        self.win_rates = {event1:0,event2:0,'unknown':0}
        self._hash = hash(":".join([self.event1, self.event2]))

    def calculate_rates(self):
        total_cooccurrences = sum(self.win_rates.values())
        self.num_cooccur = total_cooccurrences
        if total_cooccurrences > 0:
            self.mut1_win_rate = float(self.win_rates[self.event1]) / total_cooccurrences
            self.mut2_win_rate = float(self.win_rates[self.event2]) / total_cooccurrences
            self.draw_rate = float(self.win_rates['unknown']) / total_cooccurrences

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return self is other

class Season():

    def __init__(self, league, season_number):
        self.league = league
        self.table = {}
        for event in self.league.final_event_list:
            self.table[event] = 0
        self._hash = hash(":".join(self.league.final_event_list) + ":" + str(season_number))

    def __hash__(self):
        return self._hash

    def update_table_from_multinomial_sampling(self, event_pair, multinomial_draw):
        # 2 points for a win, 1 point a draw, 0 points for a loss
        if multinomial_draw == [1, 0, 0]:
            self.table[event_pair[0]] += 2
        elif multinomial_draw == [0, 1, 0]:
            self.table[event_pair[1]] += 2
        else:
            self.table[event_pair[0]] += 1
            self.table[event_pair[1]] += 1

    def get_sorted_league(self):
        self.sorted_league = sorted(self.table.items(), key=operator.itemgetter(1), reverse=True)

    def get_event_positions(self):
        self.league_event_pos = {}
        prev_score = -1
        pos = 1
        pos_counter = 1
        for j,(eve,score) in enumerate(self.sorted_league):
            if score == prev_score:
                self.league_event_pos[eve] = pos
                pos_counter += 1
            else:
                self.league_event_pos[eve] = pos_counter
                prev_score = score
                pos_counter += 1
                pos = pos_counter

    def extra_time(self, event1, event2):
        return 0

    def pk_shootout(self, event1, event2):
        return 0

class League():

    ## params ##
    arm_CNVs = set([])
    for chrom in list(map(str, range(23))) + ["X", "Y"]:
        arm_CNVs.add('loss_' + chrom + 'p')
        arm_CNVs.add('gain_' + chrom + 'p')
        arm_CNVs.add('loss_' + chrom + 'q')
        arm_CNVs.add('gain_' + chrom + 'q')
        arm_CNVs.add('loss_' + chrom)
        arm_CNVs.add('gain_' + chrom)
    arm_CNVs.add('WGD')

    def __init__(
            self, query_res_df, cohort=None, final_event_list=None, keep_only_samples_w_event=None,
            remove_samps_w_event=None, keep_samps=None, remove_samps=None, num_games_against_each_opponent=2, final_samples=None,
            max_num_snvs=20, max_num_cnv_focal=20, max_num_homdel=5, max_num_cnv_arm_losses=15, max_num_cnv_arm_gains=15, max_num_cnv_arms=30, min_event_prevalence=0.05,
            mutation_event_resolution='gene', hotspot_fn=None, drop_silent_mutations=None
    ):

        self.query_res_df = _normalize_comp_columns(query_res_df)
        self.seasons = []
        self.odds = {}
        self.event_pos = []
        self.cohort = cohort
        self.num_games_against_each_opponent = num_games_against_each_opponent
        self.mutation_event_resolution = mutation_event_resolution if mutation_event_resolution in ['gene', 'hotspot_nonsyn'] else 'gene'
        self.drop_silent_mutations = bool(drop_silent_mutations) if drop_silent_mutations is not None else (self.mutation_event_resolution != 'gene')
        self.hotspot_fn = hotspot_fn
        self._hotspot_single_residue, self._hotspot_indel_ranges = _load_hotspot_index(hotspot_fn) if self.mutation_event_resolution != 'gene' else (set(), {})
        self.event_label_map = {}
        logging.info('loading df')
        self.load_df()
        logging.info('subsetting to earliest mut')
        self.subset_to_earliest_point_mut()
        if final_samples is None:
            self.full_n_samps = len(self.events_per_samp)
        else:
            self.full_n_samps = len(final_samples)
        if keep_only_samples_w_event is not None:
            self.remove_samps_wout_event(keep_only_samples_w_event)
        if remove_samps_w_event is not None:
            self.remove_samps_with_event(remove_samps_w_event)
        if keep_samps is not None:
            self.keep_samps(keep_samps)
        if remove_samps is not None:
            self.remove_samps(remove_samps)
        logging.info('updating pairwise probs')
        self.update_pairwise_probs()
        logging.info('getting final event list')
        if final_event_list is not None:
            # Filter to events actually present in the comparison data.
            # Events not in mut_type were never seen in load_df(), so they
            # have no pairwise probabilities and cannot participate.
            unknown_events = [e for e in final_event_list if e not in self.mut_type]
            if unknown_events:
                logging.warning(
                    'Dropping %d event(s) from final_event_list not found in '
                    'comparison data: %s', len(unknown_events), unknown_events)
            self.final_event_list = [e for e in final_event_list if e in self.mut_type]
            self.calc_event_occur(final_event_list=self.final_event_list)
        else:
            self.calc_event_occur()
            self.final_event_list = self.get_final_event_list(
                max_mut=max_num_snvs, max_focal=max_num_cnv_focal, max_homdel=max_num_homdel,
                num_gains_default=max_num_cnv_arm_gains, num_losses_default=max_num_cnv_arm_losses,
                max_arm=max_num_cnv_arms, min_prevalence=min_event_prevalence
            )
        with open(f"{cohort}.final_event_list.txt", "w") as f:
            for e in sorted(self.final_event_list):
                f.write(f"{e}\n")
        self.form_pairs_for_league_model()
        self.update_pairs_for_league_model()
        self.run_league_model_iter(num_seasons=1000)
        self.n_perms = 0 # n_perms updated when updating odds
        self._hash = self.cohort

    def __hash__(self):
        return self._hash

    def add_season(self,season):
        self.seasons.append(season)

    def calc_event_occur(self,samples=None,final_event_list=None):

        self.event_prev = {}
        if samples is not None: sample_list = samples
        else: sample_list = self.events_per_samp.keys()

        if final_event_list is not None:
            for eve in final_event_list:
                if eve not in self.event_prev:
                    self.event_prev[eve] = 0

        for samp in sample_list:
            for eve in self.events_per_samp[samp]:
                if eve not in self.event_prev:
                    self.event_prev[eve] = 0
                self.event_prev[eve] += 1

        self.n_samps = len(sample_list)

    def remove_samps_wout_event(self,event):

        for samp in self.events_per_samp.keys():
            if event not in self.events_per_samp[samp]:
                del self.events_per_samp[samp]
                del self.events_per_samp_full[samp]
                del self.num_comp[samp]
                del self.event_pairs_per_samp_full[samp]

    def remove_samps_with_event(self,event):

        for samp in self.events_per_samp.keys():
            if event in self.events_per_samp[samp]:
                del self.events_per_samp[samp]
                del self.events_per_samp_full[samp]
                del self.num_comp[samp]
                del self.event_pairs_per_samp_full[samp]

    def keep_samps(self,samps):

        n_samps = 0
        output = open(self.cohort+'.final_samples.tsv','w')
        for samp in self.events_per_samp.keys():
            if samp not in samps:
                del self.events_per_samp[samp]
                del self.events_per_samp_full[samp]
                del self.num_comp[samp]
                del self.event_pairs_per_samp_full[samp]
            else:
                output.write(samp+'\n')
                n_samps += 1
        output.close()
        self.full_n_samps = n_samps

    def remove_samps(self,samps):

        for samp in self.events_per_samp.keys():
            if samp in samps:
                del self.events_per_samp[samp]
                del self.events_per_samp_full[samp]
                del self.num_comp[samp]
                del self.event_pairs_per_samp_full[samp]

    def get_samp_matrix_split(self,event1,event2):

        samps_w_both = [samp for samp in self.events_per_samp.keys() if event1 in self.events_per_samp[samp]
                        and event2 in self.events_per_samp[samp]]
        samps_w_event1_only = [samp for samp in self.events_per_samp.keys() if event1 in self.events_per_samp[samp]
                        and event2 not in self.events_per_samp[samp]]
        samps_w_event2_only = [samp for samp in self.events_per_samp.keys() if event1 not in self.events_per_samp[samp]
                        and event2 in self.events_per_samp[samp]]
        samps_w_neither = [samp for samp in self.events_per_samp.keys() if event1 not in self.events_per_samp[samp]
                        and event2 not in self.events_per_samp[samp]]
        return samps_w_both, samps_w_event1_only, samps_w_event2_only, samps_w_neither

    def _is_point_mutation_event(self, event):
        event = str(event)
        if event in self.arm_CNVs or event == 'WGD':
            return False
        if event.split('_')[0] in ['loss', 'gain', 'homdel', 'Focal']:
            return False
        return True

    def _is_hotspot_event(self, gene, mutation_token):
        gene = _normalize_gene_token(gene) if gene is not None else None
        if not gene:
            return False
        parsed = _classify_point_mutation((gene + '_' + mutation_token) if mutation_token else gene)
        position = parsed.get('position')
        if position is not None and (gene, int(position)) in self._hotspot_single_residue:
            return True
        if mutation_token:
            ints = [int(x) for x in _INT_TOKEN_RE.findall(str(mutation_token))]
            if ints:
                mut_start, mut_end = min(ints), max(ints)
                for start, end in self._hotspot_indel_ranges.get(gene, []):
                    if not (mut_end < start or mut_start > end):
                        return True
        return False

    def _harmonize_event(self, event, event_class=None):
        raw_event = str(event)

        # ── If an explicit class annotation was supplied (from the
        #    event1_class / event2_class columns in the comparisons table),
        #    use it directly for point-mutation events instead of running the
        #    internal classification logic.
        if event_class is not None and self._is_point_mutation_event(raw_event):
            event_class = str(event_class).strip()
            parsed = _classify_point_mutation(raw_event)
            gene = parsed.get('gene')
            mutation_token = parsed.get('mutation_token')
            harmonized_event = (gene + ':' + event_class) if gene else raw_event
            is_silent = (event_class.upper() == _CODING_SILENT_CLASS)
            is_hotspot = ('HOTSPOT' in event_class.upper())
            keep = not is_silent or not self.drop_silent_mutations
            info = {
                'raw_event': raw_event,
                'harmonized_event': harmonized_event,
                'family_label': harmonized_event,
                'gene': gene,
                'event_type': 'snv',
                'resolution_mode': 'external_class',
                'contains_hotspot': is_hotspot,
                'coding_class': event_class,
                'keep': keep,
                'mutation_token': mutation_token,
            }
            self.event_label_map[raw_event] = info
            return info

        if self.mutation_event_resolution == 'gene':
            harmonized_event = _legacy_event_gene(raw_event)
            if harmonized_event in self.arm_CNVs:
                event_type = 'arm_level'
            elif 'loss' in harmonized_event or 'gain' in harmonized_event or 'homdel' in harmonized_event:
                event_type = 'focal_level'
            elif harmonized_event == 'WGD':
                event_type = 'WGD'
            else:
                event_type = 'snv'
            info = {
                'raw_event': raw_event,
                'harmonized_event': harmonized_event,
                'family_label': harmonized_event,
                'gene': harmonized_event,
                'event_type': event_type,
                'resolution_mode': self.mutation_event_resolution,
                'contains_hotspot': False,
                'coding_class': None,
                'keep': True,
            }
            self.event_label_map[raw_event] = info
            return info

        if not self._is_point_mutation_event(raw_event):
            harmonized_event = _legacy_event_gene(raw_event)
            if harmonized_event in self.arm_CNVs:
                event_type = 'arm_level'
            elif 'loss' in harmonized_event or 'gain' in harmonized_event or 'homdel' in harmonized_event:
                event_type = 'focal_level'
            elif harmonized_event == 'WGD':
                event_type = 'WGD'
            else:
                event_type = 'snv'
            info = {
                'raw_event': raw_event,
                'harmonized_event': harmonized_event,
                'family_label': harmonized_event,
                'gene': harmonized_event,
                'event_type': event_type,
                'resolution_mode': self.mutation_event_resolution,
                'contains_hotspot': False,
                'coding_class': None,
                'keep': True,
            }
            self.event_label_map[raw_event] = info
            return info

        parsed = _classify_point_mutation(raw_event)
        gene = parsed.get('gene')
        mutation_token = parsed.get('mutation_token')
        is_silent = bool(parsed.get('is_silent'))
        is_hotspot = bool(self._is_hotspot_event(gene, mutation_token))
        if is_hotspot and gene:
            harmonized_event = gene + ':HOTSPOT'
            coding_class = 'HOTSPOT'
            keep = True
        elif is_silent:
            harmonized_event = (gene + ':SILENT') if gene else raw_event
            coding_class = _CODING_SILENT_CLASS
            keep = not self.drop_silent_mutations
        else:
            harmonized_event = (gene + ':NONSYN') if gene else raw_event
            coding_class = _CODING_NONSYN_CLASS
            keep = True
        info = {
            'raw_event': raw_event,
            'harmonized_event': harmonized_event,
            'family_label': harmonized_event,
            'gene': gene,
            'event_type': 'snv',
            'resolution_mode': self.mutation_event_resolution,
            'contains_hotspot': is_hotspot,
            'coding_class': coding_class,
            'keep': keep,
            'mutation_token': mutation_token,
        }
        self.event_label_map[raw_event] = info
        return info

    def load_df(self):

        self.events_per_samp = {}
        self.events_per_samp_full = {}
        self.event_pairs_per_samp_full = {}
        self.num_comp = {}
        self.mut_type = {}
        self.gene_names = {}

        # Check whether the comparisons table carries pre-computed class
        # annotations.  When present we pass them through to
        # _harmonize_event so it uses the external labels instead of its
        # own classification logic.
        _has_class_cols = (
            'event1_class' in self.query_res_df.columns
            and 'event2_class' in self.query_res_df.columns
        )
        if _has_class_cols:
            logging.info(
                'event1_class / event2_class columns detected in comparisons '
                'table – using external class annotations for harmonised '
                'event keys'
            )

        for i, row in self.query_res_df.iterrows():

            samp = row['sample']
            event1 = str(row['event1'])
            event2 = str(row['event2'])
            p_event1_win = float(row['p_event1_win'])
            p_event2_win = float(row['p_event2_win'])
            p_unknown = float(row['unknown'])

            if samp not in self.events_per_samp:
                self.events_per_samp[samp] = set([])
                self.events_per_samp_full[samp] = set([])
                self.event_pairs_per_samp_full[samp] = {}
                self.num_comp[samp] = {}

            # NaN != NaN evaluates True, so this safely skips missing values.
            ev1_class = None
            ev2_class = None
            if _has_class_cols:
                v1 = row['event1_class']
                if v1 == v1 and v1 is not None:      # not NaN / not None
                    ev1_class = str(v1).strip() or None
                v2 = row['event2_class']
                if v2 == v2 and v2 is not None:
                    ev2_class = str(v2).strip() or None

            event1_info = self._harmonize_event(event1, event_class=ev1_class)
            event2_info = self._harmonize_event(event2, event_class=ev2_class)
            event1_gene = event1_info['harmonized_event']
            event2_gene = event2_info['harmonized_event']
            self.gene_names[event1] = event1_gene
            self.gene_names[event2] = event2_gene

            if event1_info['keep']:
                self.events_per_samp[samp].add(event1_gene)
                self.events_per_samp_full[samp].add(event1)
                if event1 not in self.num_comp[samp].keys():
                    self.num_comp[samp][event1] = [0, 0]
                self.num_comp[samp][event1][0] += p_event1_win
                self.mut_type[event1_gene] = event1_info['event_type']

            if event2_info['keep']:
                self.events_per_samp[samp].add(event2_gene)
                self.events_per_samp_full[samp].add(event2)
                if event2 not in self.num_comp[samp].keys():
                    self.num_comp[samp][event2] = [0, 0]
                self.num_comp[samp][event2][1] += p_event2_win
                self.mut_type[event2_gene] = event2_info['event_type']

            if event1_info['keep'] and event2_info['keep']:
                sorted_pair_full = tuple(sorted([event1, event2]))
                self.event_pairs_per_samp_full[samp][sorted_pair_full] = {
                    (event1, event2): p_event1_win,
                    (event2, event1): p_event2_win,
                    'unknown': p_unknown,
                }

    def get_samps_w_event(self,event):
        return [samp for samp in self.events_per_samp.keys() if event in self.events_per_samp[samp]]

    # subsetting to earliest point mutations in cases of multi-hits
    def subset_to_earliest_point_mut(self):

        self.final_events_full = {}
        for samp in self.events_per_samp_full.keys():

            self.final_events_full[samp] = []
            num_hits = {}
            point_muts = {}

            for mut in self.events_per_samp_full[samp]:
                if self.mut_type[self.gene_names[mut]] in ['arm_level','WGD','focal_level']:
                    self.final_events_full[samp].append(mut)
                else:
                    mut_gene = self.gene_names[mut]
                    if mut_gene not in num_hits:
                        num_hits[mut_gene] = 0
                        point_muts[mut_gene] = []
                    num_hits[mut_gene] += 1
                    point_muts[mut_gene].append(mut)

            for mut_gene in num_hits.keys():
                if num_hits[mut_gene] == 1:
                    self.final_events_full[samp].append(point_muts[mut_gene][0])
                else:
                    earliest_mut = sorted([mut for mut in self.num_comp[samp].items() if
                                           self.gene_names[mut[0]] == mut_gene],
                                          key = lambda x:(-x[1][1],x[1][0]))[0][0]
                    self.final_events_full[samp].append(earliest_mut)

    def update_pairwise_probs(self):

        self.event_pairs_per_samp = {}
        for samp in self.event_pairs_per_samp_full.keys():
            self.event_pairs_per_samp[samp] = {}
            for eve_pair in self.event_pairs_per_samp_full[samp].keys():

                eve1 = eve_pair[0]
                eve2 = eve_pair[1]
                eve1_gene = self.gene_names[eve1]
                eve2_gene = self.gene_names[eve2]
                gene_pair = tuple(sorted([eve1_gene,eve2_gene]))
                pair_probs = self.event_pairs_per_samp_full[samp][eve_pair]

                if eve1 in self.final_events_full[samp] and eve2 in self.final_events_full[samp]:

                    self.event_pairs_per_samp[samp][gene_pair] = {(eve1_gene,eve2_gene):pair_probs[(eve1,eve2)],
                                                                  (eve2_gene,eve1_gene):pair_probs[(eve2,eve1)],
                                                                  'unknown':pair_probs['unknown']}

    def get_final_event_list(self,max_mut=20,max_hotspot=20,max_focal=20, max_homdel=5, num_gains_default=15,
                             num_losses_default=15,max_arm=30, min_prevalence=0.05):

        final_event_list = []

        num_arm = 0
        # first, add top occuring gains
        for i,(eve,n_occur) in enumerate(sorted([x for x in self.event_prev.items()
                                                 if self.mut_type[x[0]] in ['arm_level'] and
                                                    'gain' in x[0]], key = lambda y:y[1],reverse=True)):
            if i >= num_gains_default: break
            if n_occur/float(self.n_samps) < min_prevalence: break
            final_event_list.append(eve)
            num_arm += 1

        # then, add top occuring losses
        for i,(eve,n_occur) in enumerate(sorted([x for x in self.event_prev.items()
                                                 if self.mut_type[x[0]] in ['arm_level'] and
                                                    'loss' in x[0]], key = lambda y:y[1],reverse=True)):
            if i >= num_losses_default: break
            if n_occur/float(self.n_samps) < min_prevalence: break
            final_event_list.append(eve)
            num_arm += 1

        # then, add rest of arm events
        for i,(eve,n_occur) in enumerate(sorted([x for x in self.event_prev.items()
                                                 if self.mut_type[x[0]] in ['arm_level'] and
                                                    x[0] not in final_event_list], key = lambda y:y[1],reverse=True)):
            if i+num_arm >= max_arm: break
            if n_occur/float(self.n_samps) < min_prevalence: break
            final_event_list.append(eve)

        # then, add hotspots
        hotspot_indel = sorted([x for x in self.event_prev.items() if self.mut_type[x[0]] in ['snv'] and "HOTSPOT" in x[0]], key = lambda y:y[1], reverse=True)
        for i,(eve,n_occur) in enumerate(hotspot_indel):
            if i >= max_hotspot: break
            if n_occur/float(self.n_samps) < min_prevalence: break
            final_event_list.append(eve)

        # then, add snvs/indels
        snv_indel = sorted([x for x in self.event_prev.items() if self.mut_type[x[0]] in ['snv']], key = lambda y:y[1], reverse=True)
        for i,(eve,n_occur) in enumerate(snv_indel):
            if i >= max_mut: break
            if n_occur/float(self.n_samps) < min_prevalence: break
            final_event_list.append(eve)

        # then, add focal events
        for i,(eve,n_occur) in enumerate(sorted([x for x in self.event_prev.items()
                                                 if self.mut_type[x[0]] in ['focal_level'] and
                                                    'homdel' not in x[0]], key = lambda y:y[1],reverse=True)):
            if i >= max_focal: break
            if n_occur/float(self.n_samps) < min_prevalence: break
            final_event_list.append(eve)

        # then, add homdels
        for i,(eve,n_occur) in enumerate(sorted([x for x in self.event_prev.items()
                                                 if self.mut_type[x[0]] in ['focal_level'] and
                                                    'homdel' in x[0]], key = lambda y:y[1],reverse=True)):
            if i >= max_homdel: break
            if n_occur/float(self.n_samps) < min_prevalence: break
            final_event_list.append(eve)

        # then, add WGD if present
        if 'WGD' in self.event_prev:
            final_event_list.append('WGD')

        final_event_list = list(set(final_event_list))
        return final_event_list

    def form_pairs_for_league_model(self,final_event_list=None):

        self.event_pairs = {}
        if final_event_list is not None:
            events = final_event_list
        else:
            events = self.final_event_list
        for eve1,eve2 in itertools.combinations(events, 2):
            sorted_pair = tuple(sorted([eve1,eve2]))
            self.event_pairs[sorted_pair] = Eve_Pair(sorted_pair[0],sorted_pair[1],self.mut_type.get(eve1, 'snv'),self.mut_type.get(eve2, 'snv'))

    def update_pairs_for_league_model(self,samples=None):

        if samples is not None: sample_list = samples
        else: sample_list = self.event_pairs_per_samp.keys()

        for samp in sample_list:
            for pair in self.event_pairs_per_samp[samp].keys():
                eve1,eve2 = pair
                sorted_pair = tuple(sorted([eve1,eve2]))
                if sorted_pair not in self.event_pairs: continue
                self.event_pairs[sorted_pair].win_rates[eve1] += self.event_pairs_per_samp[samp][pair][(eve1,eve2)]
                self.event_pairs[sorted_pair].win_rates[eve2] += self.event_pairs_per_samp[samp][pair][(eve2,eve1)]
                self.event_pairs[sorted_pair].win_rates['unknown'] += self.event_pairs_per_samp[samp][pair]['unknown']

        for eve_pair in self.event_pairs.keys():
            eve1,eve2 = eve_pair
            if sum(self.event_pairs[eve_pair].win_rates.values()) < 2:
                if self.event_prev[eve1] >= 4 and self.event_prev[eve2] >= 4:
                    self.event_pairs[eve_pair].win_rates['unknown'] += 1
                else:
                    self.event_pairs[eve_pair].win_rates[eve1] += 1
                    self.event_pairs[eve_pair].win_rates[eve2] += 1

        for eve_pair in self.event_pairs.keys():
            self.event_pairs[eve_pair].calculate_rates()

    def run_league_model_iter(self,num_seasons):

        self.event_positions = {}
        for eve in self.final_event_list:
            self.event_positions[eve] = []

        for j,season in enumerate(range(num_seasons)):

            new_season = Season(self,season)
            for k in range(self.num_games_against_each_opponent):
                for event_pair in self.event_pairs.keys():
                    #if event_pair[0] not in final_event_list or event_pair[1] not in final_event_list: continue
                    multinomial_draw = list(np.random.multinomial(1, [self.event_pairs[event_pair].mut1_win_rate,
                                                                      self.event_pairs[event_pair].mut2_win_rate,
                                                                      self.event_pairs[event_pair].draw_rate], size=1)[0])
                    new_season.update_table_from_multinomial_sampling(event_pair, multinomial_draw)

            new_season.get_sorted_league()
            new_season.get_event_positions()
            for eve in new_season.league_event_pos.keys():
                self.event_positions[eve].append(new_season.league_event_pos[eve])
            self.add_season(season)

    def calc_odds(self):

        h1 = len(self.final_event_list) / 2.
        odds_dict = {}
        arr_final_pos = np.array([])

        for event in self.event_positions.keys():
            hist = np.array(self.event_positions[event])
            arr_final_pos = np.concatenate((arr_final_pos, hist), axis=0)
            #odds_early = ( max(float(np.size(np.where(hist < q1))),1.) / (
            #max(float(np.size(np.where(hist >= q1))), 1)) ) / (max(float(np.size(np.where(hist > q4))),1.) / (
            #max(float(np.size(np.where(hist <= q4))), 1)) )

            odds_early = max(float(np.size(np.where(hist < h1))),1.) / max(float(np.size(np.where(hist >= h1))),1.)
            odds_late = max(float(np.size(np.where(hist >= h1))),1.) / max(float(np.size(np.where(hist < h1))),1.)

            #odds_early = max(float(np.size(np.where(hist <= q1))),1.) / max(float(np.size(np.where(hist >= q4))),1.)
            #odds_late = max(float(np.size(np.where(hist >= q4))),1.) / max(float(np.size(np.where(hist <= q1))),1.)
            #odds_late = 1./odds_early
            odds_dict[event] = {'odds_early':odds_early,'odds_late':odds_late}

        return odds_dict

    def init_odds(self,events):
        for eve in events:
            self.odds[eve] = {'odds_early':[],'odds_late':[]}

    def run_permutation(self,num_seasons,samples=None,final_event_list=None):

        self.num_seasons = num_seasons # <-- TODO: don't update every time, limits odds score in each iteration
        self.seasons = [] # <-- need to reset seasons
        self.event_pos = []
        self.calc_event_occur(samples,final_event_list)
        self.form_pairs_for_league_model(final_event_list)
        self.update_pairs_for_league_model(samples)
        self.run_league_model_iter(num_seasons=num_seasons)

    def update_odds(self):

        odds_dict = self.calc_odds()
        for event in odds_dict.keys():
            self.odds[event]['odds_early'].append(odds_dict[event]['odds_early'])
            self.odds[event]['odds_late'].append(odds_dict[event]['odds_late'])
        self.n_perms += 1

    def run_full_run(self,num_seasons,samples=None,final_event_list=None):

        self.run_permutation(num_seasons=num_seasons,samples=samples,final_event_list=final_event_list)
        odds_dict = self.calc_odds()
        self.odds_full_run = {}
        for event in odds_dict.keys():
            self.odds_full_run[event] = {'odds_early':odds_dict[event]['odds_early'],
                                         'odds_late':odds_dict[event]['odds_late']}
        self.full_event_prev = copy.deepcopy(self.event_prev)

    def calc_log_odds_full_run(self):

        self.log_odds_full_run = {}
        for eve in self.final_event_list:
            self.log_odds_full_run[eve] = np.log(np.array(self.odds[eve]['odds_early']))/np.log(10)

    ####################
    ##### plotting #####
    ####################

    ## helper function ##
    def autolabelh(self, rects, labels=None):
        # attach some text labels
        for i, rect in enumerate(rects):
            height = rect.get_width()
            label = str(int(height)) if labels is None else str(labels[i])
            plt.text(height + 1.5, rect.get_xy()[1] + 0.3, str(label), ha='left', va='center', fontsize=4)

    def plot_league_run(self,type='odds'):

        if type == 'odds':
            odds_plot = plt.figure(figsize=(8.0, 10.0))
        elif type == 'pos':
            pos_plot = plt.figure(figsize=(8.0, 10.0))
        else:
            return 0

        # fig4 is the odds plot (no clonal bars)
        sns.set(font_scale=1.0)
        sns.set_style('white')
        gs = gridspec.GridSpec(1, 10,wspace=0.05, hspace=0.05)
        ax0 = plt.subplot(gs[1:8])
        #plt.gcf().subplots_adjust(left=0.15)

        if type == 'odds':
            sorted_medians = [(y[0],y[1]) for y in sorted([(eve,np.median(self.log_odds_full_run[eve]),
                                                            min(self.log_odds_full_run[eve]),
                                                            max(self.log_odds_full_run[eve]))
                                                           for eve in self.log_odds_full_run.keys()],
                                                          key = lambda x:(x[1],x[2],x[3]),reverse=True)]
        elif type == 'pos':
            sorted_medians = [(y[0],y[1]) for y in sorted([(eve,np.median(self.event_positions[eve]),
                                                            min(self.event_positions[eve]),
                                                            max(self.event_positions[eve]))
                                                           for eve in self.final_event_list
                                                           if len(self.event_positions.get(eve, [])) > 0],
                                                          key = lambda x:(x[1],x[2],x[3]))]


        colors_v = []
        for med in sorted_medians:
            if 'loss' in med[0]: colors_v.append("#3498db")
            elif 'homdel' in med[0]: colors_v.append(sns.xkcd_rgb["royal blue"])
            elif med[0] == 'WGD': colors_v.append("black")
            elif 'gain' in med[0]: colors_v.append("#e74c3c")
            else: colors_v.append(sns.xkcd_rgb["faded green"])

        colors_l = ['black'] * len(colors_v)
        if type == 'odds':
            safe_num_seasons = max(self.num_seasons, 1)
            x_in = np.linspace(np.log10(1./safe_num_seasons), np.log10(safe_num_seasons), 10000 + 1)
        elif type == 'pos':
            x_in = np.linspace(1, len(self.final_event_list), 10000 + 1)

        for j, med in enumerate(sorted_medians):

            lw = 1
            if type == 'odds':
                to_plot = -np.log10(np.array(self.odds[med[0]]['odds_early']))
            elif type == 'pos':
                to_plot = self.event_positions[med[0]]

            if len(to_plot) == 0:
                continue
            if max(to_plot) == min(to_plot):
                to_plot = list(to_plot)
                to_plot.append(min(to_plot) - 0.1)
                to_plot.append(max(to_plot) + 0.1)
                to_plot = np.array(to_plot)
            density = gaussian_kde(to_plot)
            density.covariance_factor = lambda: 0.4
            density._compute_covariance()
            y = density(x_in)
            y_new = y / np.sum(y)

            left_sum = 0
            left_idx = 0
            for i in range(len(x_in)):
                left_sum += y_new[i]
                if left_sum >= 0.000005: break
                else: left_idx = i
            right_sum = 0
            right_idx = len(x_in)-1
            for i in range(len(x_in)-1,-1,-1):
                right_sum += y_new[i]
                if right_sum >= 0.000005: break
                else: right_idx = i

            y_new = y_new / max(y_new) * 0.2
            plt.fill_between(x_in[left_idx:right_idx],y_new[left_idx:right_idx] + len(sorted_medians) - j - 1,
                             -y_new[left_idx:right_idx] + len(sorted_medians) - j - 1, color='gray',
                             alpha=0.8, lw=lw)

        sorted_medians.reverse()
        colors_v.reverse()

        ax0.set_yticks([i for i in range(0, len(sorted_medians))])
        ax0.set_yticklabels([str(" ".join(med[0].split("_"))) for med in sorted_medians],
                            fontsize=15. * 35 / max(30, len(self.final_event_list)))
        for ytick, color in zip(ax0.get_yticklabels(), colors_l): ytick.set_color(color)
        plt.title(self.cohort + ', N(samp) = ' + str(self.full_n_samps)+', N(events) = ' +
                  str(len(self.final_event_list)), fontsize=12)
        if type == 'odds':
            plt.xlabel('relative log odds timing', fontsize=12)
            safe_num_seasons = max(self.num_seasons, 1)
            plt.xlim([np.log10(1./safe_num_seasons)-0.1, np.log10(safe_num_seasons)+0.1])
        elif type == 'pos':
            plt.xlabel('event position', fontsize=12)
            plt.xlim([0, len(self.final_event_list)+1])
        plt.ylabel('event', fontsize=12)
        plt.ylim([-0.5, len(sorted_medians) - 0.5])
        ax2 = plt.subplot(gs[8:])
        plt.title('prevalence', fontsize=10)
        sns.set_style('white')
        bars = ax2.barh(np.array(range(len(sorted_medians))),
                        [self.full_event_prev[event[0]] for event in sorted_medians], color=colors_v, align='center')
        labels = [str(int(round(
            [self.full_event_prev[event[0]] for event in
             sorted_medians][i] / float(self.full_n_samps), 2) * 100)) + "%" for i in range(len(sorted_medians))]

        self.autolabelh(bars, labels)
        sns.despine(ax=ax0, left=True, bottom=True)
        plt.axis("off")
        plt.ylim([-0.5, len(sorted_medians) - 0.5])

        if type == 'odds':
            self.odds_plot = odds_plot
        elif type == 'pos':
            self.pos_plot = pos_plot

    ####################
    ###### outputs #####
    ####################
