##########################################################
# Patient - central class to handle and store sample CCF data
##########################################################
from data.Sample import TumorSample
from data.Sample import RNASample
from data.SomaticEvents import SomMutation, CopyNumberEvent
from data.Enums import Genome
from collections import defaultdict, Counter
import os
import sys
import logging
import itertools
import networkx as nx

import numpy as np
from intervaltree import Interval, IntervalTree

# Logsumexp options
import pkgutil

if pkgutil.find_loader('sselogsumexp') is not None:
    from sselogsumexp import logsumexp
else:
    try:
        from scipy.misc import logsumexp
    except ImportError:
        from scipy.special import logsumexp


class Patient:
    """CLASS info

     FUNCTIONS:
         public:

           add_sample() - add new samples from individual to self.sample_list by making a new TumorSample instance

        private:

           _auto_file_type() - guess input file type for ccf data if not manually provided

     PROPERTIES:
         regular variables:
            self.sample_list

         @property methods: none
    """

    def __init__(self, indiv_name='Indiv1',
                 ccf_grid_size=101,
                 driver_genes_file=os.path.join(os.path.dirname(__file__), 'supplement_data/Driver_genes_v1.0.txt'),
                 impute_missing=False,
                 artifact_blacklist=os.path.join(os.path.dirname(__file__), 'supplement_data/Blacklist_SNVs.txt'),
                 artifact_whitelist='',
                 use_indels=False,
                 min_coverage=8,
                 PoN_file=False,
                 genome_build=None):

        # DECLARATIONS
        self.indiv_name = indiv_name
        # :type : list [TumorSample]
        self.sample_list = []
        self.rna_sample_list = []

        self.samples_synchronized = False
        self.rna_samples_synchronized = False
        self.concordant_genes = []
        self.driver_genes = self._parse_driver_g_file(driver_genes_file)
        self.ccf_grid_size = ccf_grid_size

        self.PatientLevel_MutBlacklist = artifact_blacklist
        self.PatientLevel_MutWhitelist = artifact_whitelist

        # Patient configuration settings
        # flag if to impute missing variants as ccf 0
        self.impute_missing = impute_missing
        # min cov is specified both here and passed to tumor sample.
        self.min_coverage = min_coverage
        self.use_indels = use_indels
        self.PoN_file = PoN_file

        self._validate_sample_names()

        # later filled data objects
        self.ND_mutations = []

        self.genome = Genome(build=genome_build) if genome_build is not None else Genome()

        # storing of results
        # Clustering
        self.ClusteringResults = None
        self.MutClusters = None
        self.TruncalMutEvents = None
        self.MCMC_trace = None
        self.k_trace = None
        self.alpha_trace = None

        self.unclustered_muts = []

        self.concordant_cn_tree = {chrom: IntervalTree() for chrom in self.genome.CHROMS}

        # BuildTree
        self.TopTree = None
        self.TreeEnsemble = []
        self.alt_cn_states = []

    def initPatient(self):
        """ accepted input types abs; txt; sqlite3 .db # auto tab if .txt, .tsv or .tab ; abs if .Rdata; sqlite if .db """
        raise NotImplementedError

    def addRNAsample(self, filen, sample_name, input_type='auto', purity=None, timepoint=None):
        """ Accepted input types rsem.genes """
        # make new sample and add to exiting list of samples
        logging.info("Adding expression from RNA Sample: %s", sample_name)
        new_rna_sample = RNASample(filen, input_type, indiv=self.indiv_name, sample_name=sample_name,
                                   timepoint=timepoint, purity=purity)
        self.rna_sample_list.append(new_rna_sample)
        logging.info('Added RNA sample ' + new_rna_sample.sample_name)
        # turn of concordance flag when new sample is added
        self.rna_samples_synchronized = False

    def preprocess_rna_samples(self):

        count_low_exp_genes = 0
        gene_ids = self.rna_sample_list[0].get_gene_ids()
        for gene_id in gene_ids:
            tpm_values = []
            for sample in self.rna_sample_list:
                tpm_values.append(sample.get_tpm_by_gene_id(gene_id))
            if all(i <= 1.0 for i in tpm_values):
                count_low_exp_genes += 1
            else:
                self.concordant_genes.append(gene_id)
        self.rna_samples_synchronized = True
        logging.debug('{} genes have TPM values less a than 1.0 across all timepoints.'.format(count_low_exp_genes))
        logging.debug('{} genes will be used in the analysis.'.format(len(self.concordant_genes)))

    def addSample(self, filen, sample_name, input_type='auto', seg_input_type='auto', grid_size=101, seg_file=None,
                  _additional_muts=None, purity=None, timepoint_value=None, tmb=None):
        """ accepted input types abs; txt; sqlite3 .db # auto tab if .txt, .tsv or .tab ; abs if .Rdata; sqlite if .db"""

        if _additional_muts == []:
            _additional_muts = None
        elif type(_additional_muts) is list:
            _additional_muts = _additional_muts[0]

        # make new sample and add to exiting list of samples
        logging.info("Adding Mutations from Sample: %s", sample_name)
        new_sample = TumorSample(filen, input_type, seg_input_type=seg_input_type, sample_name=sample_name,
                                 artifact_blacklist=self.PatientLevel_MutBlacklist,
                                 artifact_whitelist=self.PatientLevel_MutWhitelist,
                                 ccf_grid_size=grid_size, PoN=self.PoN_file, indiv=self.indiv_name,
                                 use_indels=self.use_indels, min_coverage=self.min_coverage,
                                 _additional_muts=_additional_muts, seg_file=seg_file,
                                 purity=purity, timepoint_value=timepoint_value, tmb=tmb, genome=self.genome)

        self.sample_list.append(new_sample)
        logging.info('Added sample ' + new_sample.sample_name)
        # turn of concordance flag when new sample is added
        self.samples_synchronized = False

    def homogenize_events_across_samples(self):
        # TODO: do not understand this function, find place where it is used
        if len(self.sample_list) == 1:
            UserWarning("Found one sample! Will run 1D clustering. If intended to run ND, please fix!")

    def get_sample_byname(self, sample_name):
        """
        :param str sample_name: sample name to search for
        :return: TumorSample
        """
        for smpl in self.sample_list:
            if smpl.sample_name == sample_name:
                return smpl
        logging.warning("Sample with the name {} is not found in the sample list".format(sample_name))
        return None

    def _validate_sample_names(self, disallowed_strings=['Inner']):
        sample_names = self.sample_names
        # check that no samples have the same name
        if len(sample_names) != len(set(sample_names)):
            logging.error('Several samples appear to have identical names!')
            ValueError('Several samples appear to have identical names! This is not allowed! Please fix!')
        # check that no samples disallowed strings in its name
        for dis_str in disallowed_strings:
            for smpl_name in sample_names:
                if dis_str in smpl_name:
                    logging.error('Disallowed string found in sample name {}'.format(dis_str))
                    ValueError('Several samples appear to have identical names! This is not allowed! Please fix!')

        return True

    @property
    def sample_names(self):
        return [x.sample_name for x in self.sample_list]

    @staticmethod
    def _parse_driver_g_file(filen):
        """ read driver file as one gene per line """
        if not filen:
            return set()  # empty list
        with open(filen, 'r') as drv_file:
            drv = [x.strip() for x in drv_file.read().strip().split('\n')]
        return set(drv)

    def preprocess_samples(self):
        """ Central preprocessing of samples
            Makes sure that there is info on all mutations and everything is concordant on all samples """
        if len(self.sample_list) < 2:
            logging.warning("Only one sample in sample set! Cannot check concordance for multi-D clustering!")
            logging.warning("Cleaning this single sample!")

            # TODO: port the 1D-preprocess
            for sample in self.sample_list:  # because of imputing: has to be done at the end
                sample.concordant_variants = sample.mutations
                # sort concordant variants by presence in driver then var_str/gene name
                sample.concordant_variants.sort(
                    key=lambda x: (str(x).split('_')[0] in self.driver_genes, str(x), x.var_str))
                sample.private_mutations = [mut for mut in sample.concordant_variants if mut.alt_cnt > 0]

        self._validate_sample_names()
        try:
            blacklist = set(set.union(*[set(x.known_blacklisted_mut) for x in self.sample_list]))
        except TypeError:
            logging.info("Blacklist is not specified.")
            blacklist = set()

        # use var stings from sample1 to order all other samples
        var_str_init = self.sample_list[0].mut_varstr
        count_needed = len(self.sample_list)
        full_var_list = list(itertools.chain.from_iterable([x.mut_varstr for x in self.sample_list]))
        vars_present_in_all = set()

        for mut in var_str_init:
            if full_var_list.count(mut) == count_needed and mut not in blacklist:
                vars_present_in_all.add(mut)

        joint_temporarily_removed = set()
        for sample in self.sample_list:
            joint_temporarily_removed.update(sample.temporarily_removed)

        joint_temporarily_removed = set([mut.var_str for mut in joint_temporarily_removed])
        # Only allow one fixed CNV in concordant variants.
        first_cnv = None
        # iterate through all samples and make concordance sets, then sort them
        for sample in self.sample_list:
            if not self.impute_missing:
                # reset any previous joint results ; if imputing leave as is
                sample.concordant_variants = []
            sample.concordant_with_samples = []

            for mut in sample.mutations + list(sample.low_coverage_mutations.values()):

                # TODO: Add this as a flag
                # if mut.mut_category is not None and ("rna" in mut.mut_category.lower() or "utr" in mut.mut_category.lower()):
                #    logging.warning("REMOVING RNA MUTATION:" + str(mut))
                #    sample.artifacts_in_blacklist.append(mut)
                #    blacklist.add(mut.var_str)
                #

                if mut.var_str in blacklist:  # If blacklisted in another sample, remove mutation.
                    logging.info("Found blacklisted mutation still in sample, removing:" + str(mut))
                    if mut.var_str not in sample.known_blacklisted_mut:
                        sample.known_blacklisted_mut.add(mut.var_str)
                    continue

                elif mut.var_str in vars_present_in_all:
                    if mut.var_str not in joint_temporarily_removed:
                        if mut.type == "CNV":
                            if first_cnv is None:
                                first_cnv = mut.var_str
                        if mut.type != "CNV" or first_cnv == mut.var_str:
                            sample.concordant_variants.append(mut)
                    else:
                        continue  # nothing to do for indels in all samples

                elif mut.var_str not in vars_present_in_all:
                    if mut.var_str in joint_temporarily_removed:
                        sample.artifacts_in_blacklist.append(mut)
                        logging.warn(
                            "Forcecalling Failure? Mutation {} not in {} but otherwise graylisted already.".format(
                                str(mut), sample.sample_name))
                        continue

                    elif mut.var_str not in joint_temporarily_removed and self.impute_missing:
                        for mis_sample in [x for x in self.sample_list if x != sample]:
                            if mut not in mis_sample.concordant_variants and mut not in mis_sample.mutations and mut.var_str not in joint_temporarily_removed and ":".join(
                                    map(str,
                                        [mut.chrN, mut.pos, mut.ref, mut.alt])) not in mis_sample.known_blacklisted_mut:
                                logging.info('Imputing missing mutation as 0 CCF, sample %s, var: %s ; %s',
                                             sample.sample_name, mut, mut.var_str)
                                mis_sample.concordant_variants.append(
                                    SomMutation.from_som_mutation_zero(mut, from_sample=mis_sample))
                                mis_sample._mut_varstring_hashtable[mut.var_str] = mis_sample.concordant_variants[-1]
                        sample.concordant_variants.append(mut)

                    else:
                        logging.error(
                            'Mutation missing in some datasets, sample %s, var: %s . This mutation will be skipped. Use --impute, or forcecall the mutations.',
                            sample.sample_name, mut)

        for sample in self.sample_list:  # because of imputing: has to be done at the end
            # sort concordant variants by presence in driver then var_str/gene name
            sample.concordant_variants.sort(
                key=lambda x: (str(x).split('_')[0] in self.driver_genes, str(x), x.var_str))
            # annotate each sample with what samples were used for concordance
            sample.concordant_with_samples = self.sample_list
        # turn on concordance flag when new sample mut are joined
        self.samples_synchronized = True

    def _make_ND_histogram(self):
        if not self.samples_synchronized:
            logging.error("Could not make ND Histogram, please make sure to call preprocess_samples()")
            return False

        combined_ccf = []
        for mut in self.sample_list[0].concordant_variants:
            mut_nd_ccf = [mut.ccf_1d]
            for sample in self.sample_list[1:]:
                mut_nd_ccf.append(sample.get_mut_by_varstr(mut.var_str).ccf_1d)
            combined_ccf.append(mut_nd_ccf)

        return NDHistogram(combined_ccf, [x.var_str for x in self.sample_list[0].concordant_variants])

    def make_ND_histogram(self):
        if not self.samples_synchronized:
            logging.error("Could not make ND Histogram, please make sure to call preprocess_samples()")
            return False
        combined_ccf = []
        for mut in self.sample_list[0].concordant_variants:
            mut_nd_ccf = [mut.ccf_1d]
            for sample in self.sample_list[1:]:
                mut_nd_ccf.append(sample.get_mut_by_varstr(mut.var_str).ccf_1d)
            combined_ccf.append(mut_nd_ccf)

        return NDHistogram(combined_ccf, [x.var_str for x in self.sample_list[0].concordant_variants])

    def cluster_temp_removed(self):
        clust_CCF_results = self.ClusteringResults.clust_CCF_dens
        for mut in self.sample_list[0].low_coverage_mutations.values():
            mut_coincidence = np.ones(len(clust_CCF_results))
            for i, sample in enumerate(self.sample_list):
                try:
                    mut = sample.get_mut_by_varstr(mut.var_str)
                except KeyError:
                    logging.warning(mut.var_str + ' not called across all samples')
                    mut_coincidence.fill(np.nan)
                    break
                for ii, cluster_ccfs in enumerate(clust_CCF_results):
                    cluster_ccf = cluster_ccfs[i]
                    if abs(np.argmax(mut.ccf_1d) - np.argmax(cluster_ccf)) > 70:
                        dot = 0.
                    else:
                        dot = max(sum(np.array(mut.ccf_1d) * cluster_ccf), .0001)
                    mut_coincidence[ii] *= dot
            if np.any(mut_coincidence > 0.):
                cluster_assignment = np.argmax(mut_coincidence) + 1
                for i, sample in enumerate(self.sample_list):
                    mut = sample.get_mut_by_varstr(mut.var_str)
                    mut.cluster_assignment = cluster_assignment
                    mut.clust_ccf = clust_CCF_results[cluster_assignment - 1][i]
            else:
                print('Did not cluster ' + str(mut))
                self.unclustered_muts.append(mut.var_str)
                for sample in self.sample_list:
                    try:
                        mut = sample.low_coverage_mutations[mut.var_str]
                        sample.unclustered_muts.append(mut)
                    except KeyError:
                        print(mut.var_str + ' not found in ' + sample.sample_name)
                # for sample in self.sample_list:
                #     mut = sample.get_mut_by_varstr(mut.var_str)
                #     mut.cluster_assignment = None
                #     mut.clust_ccf = None

    @staticmethod
    def _infer_cn_category(cn_vec, hat_vec=None, eps=0.25, tau=0.5, use_ccf_weight=True):
        """
        Decide 'gain' vs 'loss' from per-sample allele CNs.
        cn_vec: tuple/list of per-sample allele CN (int)
        hat_vec: optional per-sample CCF hats (same length)
        eps: neutrality threshold for |CN-1|
        tau: required consensus fraction among non-neutral samples
        """
        import numpy as np
        cn = np.array(cn_vec, dtype=float)
        if hat_vec is None:
            w = np.ones_like(cn)
        else:
            w = np.array(hat_vec, dtype=float) if use_ccf_weight else np.ones_like(cn)

        # Non-neutral mask
        m = np.abs(cn - 1.0) > eps
        if not np.any(m):
            return None  # neutral overall

        # Signed vote
        D = np.sum(w[m] * (cn[m] - 1.0))

        # Consensus requirement
        signD = 1 if D > 0 else (-1 if D < 0 else 0)
        if signD == 0:
            return None
        agree = np.sum(np.sign(cn[m] - 1.0) == signD)
        frac = agree / float(np.sum(m))
        if frac < tau:
            return None  # ambiguous

        return 'gain' if D > 0 else 'loss'

    def intersect_cn_trees(self):
        """
        Gets copy number events from segment trees and adds them to samples

        """
        def merge_cn_events(event_segs, neighbors, R=frozenset(), X=frozenset()):
            """
            Merges copy number events on a single chromosome if they are adjacent and their ccf values are similar

            Args:
                event_segs: set of CN segments represented as tuple(bands, CNs, CCF_hats, CCF_highs, CCF_lows, allele)
                neighbors: dict mapping seg to set of neighbors (segs with similar CCFs)
                R: only populated in recursive calls
                X: only populated in recursive calls

            Returns:
                Generator for merged segs
            """
            is_max = True
            for s in itertools.chain(event_segs, X):
                if isadjacent(s, R):
                    is_max = False
                    break
            if is_max:
                bands = set.union(*(set(b[0]) for b in R))
                cns = next(iter(R))[1]
                ccf_hat = np.zeros(len(self.sample_list))
                ccf_high = np.zeros(len(self.sample_list))
                ccf_low = np.zeros(len(self.sample_list))
                arm = next(iter(R))[6]
                for seg in R:
                    ccf_hat += np.array(seg[2])
                    ccf_high += np.array(seg[3])
                    ccf_low += np.array(seg[4])
                yield (bands, cns, ccf_hat / len(R), ccf_high / len(R), ccf_low / len(R), arm)
            else:
                for s in event_segs:
                    if isadjacent(s, R):
                        for region in merge_cn_events(event_segs & neighbors[s], neighbors, R=R | {s}, X=X & neighbors[s]):
                            yield region
                        event_segs = event_segs - {s}
                        X = X | {s}

        def isadjacent(s, R):
            """
            Copy number events are adjacent if the max band of one is the same as
            or adjacent to the min band of the other

            """
            if not R:
                return True
            Rchain = list(itertools.chain(*(b[0] for b in R)))
            minR = min(Rchain)
            maxR = max(Rchain)
            mins = min(s[0])
            maxs = max(s[0])
            if mins >= maxR:
                return mins - maxR <= 1 and mins.band[0] == maxR.band[0]
            elif maxs <= minR:
                return minR - maxs <= 1 and maxs.band[0] == minR.band[0]
            else:
                return False

        c_trees = {}
        n_samples = len(self.sample_list)
        for chrom, csize in self.genome.CHROM_DICT.items():
            centromere = self.genome.CENT_DICT[chrom]
            tree = IntervalTree()
            for sample in self.sample_list:
                if sample.CnProfile:
                    tree.update(sample.CnProfile[chrom])
            tree.split_overlaps()
            tree.merge_equals(data_initializer=[], data_reducer=lambda a, c: a + [c])
            c_tree = IntervalTree(filter(lambda s: len(s.data) == n_samples, tree))
            c_trees[chrom] = c_tree
            event_segs = set()
            for seg in c_tree:
                start = seg.begin
                end = seg.end
                arm = "p" if end <= centromere else "q"
                bands = self.get_bands(chrom, start, end)
                cns_a1 = []
                cns_a2 = []
                ccf_hat_a1 = []
                ccf_hat_a2 = []
                ccf_high_a1 = []
                ccf_high_a2 = []
                ccf_low_a1 = []
                ccf_low_a2 = []
                for i, sample in enumerate(self.sample_list):
                    seg_data = seg.data[i][1]
                    cns_a1.append(seg_data['cn_a1'])
                    cns_a2.append(seg_data['cn_a2'])
                    ccf_hat_a1.append(seg_data['ccf_hat_a1'] if seg_data['cn_a1'] != 1 else 0.)
                    ccf_hat_a2.append(seg_data['ccf_hat_a2'] if seg_data['cn_a2'] != 1 else 0.)
                    ccf_high_a1.append(seg_data['ccf_high_a1'] if seg_data['cn_a1'] != 1 else 0.)
                    ccf_high_a2.append(seg_data['ccf_high_a2'] if seg_data['cn_a2'] != 1 else 0.)
                    ccf_low_a1.append(seg_data['ccf_low_a1'] if seg_data['cn_a1'] != 1 else 0.)
                    ccf_low_a2.append(seg_data['ccf_low_a2'] if seg_data['cn_a2'] != 1 else 0.)
                cns_a1 = np.array(cns_a1)
                cns_a2 = np.array(cns_a2)
                if np.all(cns_a1 == 1):
                    pass
                elif np.all(cns_a1 >= 1) or np.all(cns_a1 <= 1):
                    event_segs.add((tuple(bands), tuple(cns_a1), tuple(ccf_hat_a1), tuple(ccf_high_a1), tuple(ccf_low_a1), 'a1', arm))
                else:
                    logging.warning('Seg with inconsistent event: {}:{}:{}'.format(chrom, seg.begin, seg.end))
                if np.all(cns_a2 == 1):
                    pass
                elif np.all(cns_a2 >= 1) or np.all(cns_a2 <= 1):
                    event_segs.add((tuple(bands), tuple(cns_a2), tuple(ccf_hat_a2), tuple(ccf_high_a2), tuple(ccf_low_a2), 'a2', arm))
                else:
                    logging.warning('Seg with inconsistent event: {}:{}:{}'.format(chrom, seg.begin, seg.end))

            neighbors = {s: set() for s in event_segs}
            for seg1, seg2 in itertools.combinations(event_segs, 2):
                s1_hat = np.array(seg1[2])
                s2_hat = np.array(seg2[2])
                if (
                    seg1[6] == seg2[6] and
                    seg1[1] == seg2[1] and
                    np.all(s1_hat >= np.array(seg2[4])) and
                    np.all(s1_hat <= np.array(seg2[3])) and
                    np.all(s2_hat >= np.array(seg1[4])) and
                    np.all(s2_hat <= np.array(seg1[3]))
                ):
                    neighbors[seg1].add(seg2)
                    neighbors[seg2].add(seg1)

            event_cache = []
            if event_segs:
                for bands, cns, ccf_hat, ccf_high, ccf_low, arm in merge_cn_events(event_segs, neighbors):
                    mut_category = self._infer_cn_category(cns, ccf_hat)
                    if mut_category is None:
                        print("[DROPPED]", bands, cns, ccf_hat)
                        continue
                    mut_category = "Focal_" + mut_category
                    # mut_category = 'Focal_gain' if sum(cns) > len(self.sample_list) else 'Focal_loss'
                    a1 = (mut_category, bands) not in event_cache
                    if a1:
                        event_cache.append((mut_category, bands))
                    self._add_cn_event_to_samples(chrom, min(bands), max(bands), arm, cns, mut_category, ccf_hat, ccf_high, ccf_low, a1, dupe=not a1)
        self.concordant_cn_tree = c_trees

    def add_focal_cn_events(self, focal_regions, qval_max=None, focal_max_arm_frac=0.3, min_bp_overlap=1000, padding=100000, require_directional_consistency=False, split_across_centromere=True):
        """
        Add focal CNV events (e.g., from GISTIC2 peaks) as CopyNumberEvents.

        Args:
            focal_regions:
                - str path to a BED-like file: chrom  start  end  [type] [q] [name]
                  (type in {'AMP','DEL'} optional; q optional; extra columns ignored)
                - iterable of dicts, each with keys:
                    {'chrom', 'start', 'end'} and optional {'type','q','name'}
            qval_max: float or None. If set, skip focals with q > qval_max.
            focal_max_arm_frac: float in (0,1]. Skip regions covering >= this fraction of the arm (keep them arm-level).
            min_bp_overlap: minimum bp overlap between focal window and a CN segment to count that segment.
            padding: padding of focal window to check for overlap with CN segment to count that segment.
            require_directional_consistency: if True, for each sample+allele, if overlapping segments
                include both gains and losses (relative to CN=1), treat that sample/allele as neutral for this focal.
            split_across_centromere: if True, split windows crossing the centromere into p and q parts; otherwise skip.
        """

        # --- helpers ---
        def _normalize_focals(focals):
            """Return dict: chrom -> list of dicts {chrom,start,end,type,q,name}."""
            out = defaultdict(list)
            if isinstance(focals, str):
                with open(focals) as fin:
                    for ln in fin:
                        if not ln.strip() or ln.startswith('#') or ln.startswith("track name"):
                            continue
                        t = ln.strip().split()
                        chrom = t[0].lstrip('chr')
                        start = int(t[1])
                        end = int(t[2])
                        ftype = "AMP" if "Any-AP" in t[3] else "DEL" if "Any-DP" in t[3] else t[3].upper() if len(t) > 3 else None
                        q = float(t[4]) if len(t) > 4 and t[4].replace('.', '', 1).isdigit() else None
                        name = t[5] if len(t) > 5 else None
                        out[chrom].append({'chrom': chrom, 'start': start, 'end': end, 'type': ftype, 'q': q, 'name': name})
                return out
            # iterable
            for r in focals:
                d = dict(r)
                d['chrom'] = str(d.get('chrom') or d.get('chr')).lstrip('chr')
                d['start'] = int(d['start'])
                d['end'] = int(d['end'])
                if 'type' in d and d['type'] is not None:
                    d['type'] = str(d['type']).upper()
                if 'q' in d and d['q'] is not None:
                    d['q'] = float(d['q'])
                out[d['chrom']].append(d)
            return out

        def _merge_by_type(regs):
            """Merge overlapping intervals within AMP and within DEL; pass through untyped."""
            def _norm(t):
                if not t: return None
                u = str(t).upper()
                return 'AMP' if u.startswith('AMP') else ('DEL' if u.startswith('DEL') else None)

            out = []
            for T in ('AMP', 'DEL'):
                xs = [dict(r, type=T) for r in regs if _norm(r.get('type')) == T]
                xs.sort(key=lambda d: (int(d['start']), int(d['end'])))
                cur = None
                for r in xs:
                    s, e = int(r['start']), int(r['end'])
                    if cur is None or s > cur['end']:
                        cur = {'chrom': r['chrom'], 'start': s, 'end': e, 'type': T}
                        out.append(cur)
                    else:
                        cur['end'] = max(cur['end'], e)
            # keep untyped as-is
            out += [r for r in regs if _norm(r.get('type')) is None]
            return out

        def _len_weighted_mode(values, weights):
            acc = Counter()
            for v, w in zip(values, weights):
                acc[int(v)] += float(w)
            if not acc:
                return 1
            return max(acc.items(), key=lambda kv: kv[1])[0]

        def _summarize_sample_interval(sample, chrom, s, e, allele_key):
            """
            Summarize allele-specific CN & CCF within [s,e) for one sample.
            Returns (cn, hat, hi, lo), with CCFs zeroed if cn==1 (neutral).
            """
            cn_prof = sample.CnProfile
            if not cn_prof or chrom not in cn_prof:
                return 1, 0.0, 0.0, 0.0

            overlaps = list(cn_prof[chrom].overlap(min(0, s - padding), e + padding))
            if not overlaps:
                return 1, 0.0, 0.0, 0.0

            ws, cns, hats, his, los = [], [], [], [], []
            for iv in overlaps:
                ov = max(0, min(iv.end, e + padding) - max(iv.begin, s - padding))
                if ov < min_bp_overlap:
                    continue
                seg_data = iv.data[1] if isinstance(iv.data, (list, tuple)) else iv.data
                ws.append(ov)
                cns.append(int(seg_data[f"cn_{allele_key}"]))
                hats.append(seg_data[f"ccf_hat_{allele_key}"])
                his.append(seg_data[f"ccf_high_{allele_key}"])
                los.append(seg_data[f"ccf_low_{allele_key}"])

            if not ws:
                return 1, 0.0, 0.0, 0.0

            # optional: enforce directional consistency within this sample+allele
            if require_directional_consistency:
                dirs = {np.sign(c - 1) for c in cns if c != 1}
                if len(dirs) > 1:
                    return 1, 0.0, 0.0, 0.0

            cn_mode = _len_weighted_mode(cns, ws)
            if cn_mode == 1:
                return 1, 0.0, 0.0, 0.0

            wsum = float(sum(ws))
            hat = float(np.dot(hats, ws) / wsum)
            hi = float(np.dot(his, ws) / wsum)
            lo = float(np.dot(los, ws) / wsum)
            return int(cn_mode), hat, hi, lo

        def _split_by_centromere(chrom, start, end, centromere, chrom_size):
            if not split_across_centromere:
                # if it crosses the centromere, skip (keeps method focal-only)
                if start < centromere < end:
                    return []
                arm = 'p' if end <= centromere else 'q'
                return [(start, end, arm)]
            if start < centromere < end:
                return [(start, centromere, 'p'), (centromere, end, 'q')]
            arm = 'p' if end <= centromere else 'q'
            return [(start, end, arm)]

        # --- main ---
        focals_by_chr = _normalize_focals(focal_regions)

        for chrom, csize in self.genome.CHROM_DICT.items():
            centromere = self.genome.CENT_DICT[chrom]
            regs = focals_by_chr.get(chrom, [])
            if not regs:
                continue

            regs = _merge_by_type(regs)

            arm_len = {'p': centromere, 'q': csize - centromere}

            for R in regs:
                # q-value filter if provided in record and threshold set
                if qval_max is not None and "q" in R is not None and R['q'] > qval_max:
                    continue
                s, e = int(R['start']), int(R['end'])
                if e <= s:
                    continue

                subregs = _split_by_centromere(chrom, s, e, centromere, csize)
                for (s0, e0, arm_lbl) in subregs:
                    # focal filter vs arm fraction
                    if (e0 - s0) / float(arm_len[arm_lbl]) >= focal_max_arm_frac:
                        # this window is too large; treat as arm-level elsewhere
                        continue

                    bands = self.get_bands(chrom, s0, e0)

                    # summarize per-allele across samples
                    results = {}
                    for allele_key in ('a1', 'a2'):
                        cn_vec, hat, hi, lo = [], [], [], []
                        for sample in self.sample_list:
                            cn_i, hat_i, hi_i, lo_i = _summarize_sample_interval(sample, chrom, s0, e0, allele_key)
                            cn_vec.append(cn_i)
                            hat.append(hat_i)
                            hi.append(hi_i)
                            lo.append(lo_i)
                        results[allele_key] = (tuple(cn_vec), tuple(hat), tuple(hi), tuple(lo))

                    # emit an event for each allele with any non-neutral evidence
                    for allele_key, (cn_vec, hat, hi, lo) in results.items():
                        if all(cn == 1 for cn in cn_vec):
                            continue

                        if R['type'] in ('AMP', 'AMPLIFICATION'):
                            gistic_cn_category = 'Focal_gain'
                        elif R['type'] in ('DEL', 'DELETION'):
                            gistic_cn_category = 'Focal_loss'
                        else:
                            gistic_cn_category = ""

                        # cn_category = 'Focal_gain' if any(cn > 1 for cn in cn_vec) else 'Focal_loss'
                        cn_category = self._infer_cn_category(cn_vec, hat)
                        if cn_category is None:
                            print("[DROPPED]", bands, cn_vec, hat)
                            continue
                        cn_category = "Focal_" + cn_category

                        if cn_category == gistic_cn_category:
                            self._add_cn_event_to_samples(chrom=chrom, start=min(bands), end=max(bands), arm=arm_lbl, cns=cn_vec, cn_category=cn_category, ccf_hat=hat, ccf_high=hi, ccf_low=lo, a1=(allele_key == 'a1'), dupe=False)

    def get_arm_level_cn_events(self):
        n_samples = len(self.sample_list)
        for chrom, csize in self.genome.CHROM_DICT.items():
            centromere = self.genome.CENT_DICT[chrom]
            tree = IntervalTree()
            for sample in self.sample_list:
                if sample.CnProfile:
                    tree.update(sample.CnProfile[chrom])
            tree.split_overlaps()
            tree.merge_equals(data_initializer=[], data_reducer=lambda a, c: a + [c])
            c_tree = IntervalTree(filter(lambda s: len(s.data) == n_samples, tree))
            event_segs = set()
            for seg in c_tree:
                start = seg.begin
                end = seg.end
                bands = self.get_bands(chrom, start, end)
                cns_a1 = []
                cns_a2 = []
                ccf_hat_a1 = []
                ccf_hat_a2 = []
                ccf_high_a1 = []
                ccf_high_a2 = []
                ccf_low_a1 = []
                ccf_low_a2 = []
                for i, sample in enumerate(self.sample_list):
                    seg_data = seg.data[i][1]
                    cns_a1.append(seg_data['cn_a1'])
                    cns_a2.append(seg_data['cn_a2'])
                    ccf_hat_a1.append(seg_data['ccf_hat_a1'] if seg_data['cn_a1'] != 1 else 0.)
                    ccf_hat_a2.append(seg_data['ccf_hat_a2'] if seg_data['cn_a2'] != 1 else 0.)
                    ccf_high_a1.append(seg_data['ccf_high_a1'] if seg_data['cn_a1'] != 1 else 0.)
                    ccf_high_a2.append(seg_data['ccf_high_a2'] if seg_data['cn_a2'] != 1 else 0.)
                    ccf_low_a1.append(seg_data['ccf_low_a1'] if seg_data['cn_a1'] != 1 else 0.)
                    ccf_low_a2.append(seg_data['ccf_low_a2'] if seg_data['cn_a2'] != 1 else 0.)
                cns_a1 = np.array(cns_a1)
                cns_a2 = np.array(cns_a2)
                if np.all(cns_a1 == 1):
                    pass
                elif np.all(cns_a1 >= 1) or np.all(cns_a1 <= 1):
                    if start < centromere < end:
                        event_segs.add((start, centromere, tuple(bands), 'p', tuple(cns_a1), tuple(ccf_hat_a1), tuple(ccf_high_a1), tuple(ccf_low_a1), 'a1'))
                        event_segs.add((centromere, end, tuple(bands), 'q', tuple(cns_a1), tuple(ccf_hat_a1), tuple(ccf_high_a1), tuple(ccf_low_a1), 'a1'))
                    elif end < centromere:
                        event_segs.add((start, end, tuple(bands), 'p', tuple(cns_a1), tuple(ccf_hat_a1), tuple(ccf_high_a1), tuple(ccf_low_a1), 'a1'))
                    else:
                        event_segs.add((start, end, tuple(bands), 'q', tuple(cns_a1), tuple(ccf_hat_a1), tuple(ccf_high_a1), tuple(ccf_low_a1), 'a1'))
                else:
                    logging.warning('Seg with inconsistent event: {}:{}:{}'.format(chrom, start, end))
                if np.all(cns_a2 == 1):
                    pass
                elif np.all(cns_a2 >= 1) or np.all(cns_a2 <= 1):
                    if start < centromere < end:
                        event_segs.add((start, centromere, tuple(bands), 'p', tuple(cns_a2), tuple(ccf_hat_a2), tuple(ccf_high_a2), tuple(ccf_low_a2), 'a2'))
                        event_segs.add((centromere, end, tuple(bands), 'q', tuple(cns_a2), tuple(ccf_hat_a2), tuple(ccf_high_a2), tuple(ccf_low_a2), 'a2'))
                    elif end < centromere:
                        event_segs.add((start, end, tuple(bands), 'p', tuple(cns_a2), tuple(ccf_hat_a2), tuple(ccf_high_a2), tuple(ccf_low_a2), 'a2'))
                    else:
                        event_segs.add((start, end, tuple(bands), 'q', tuple(cns_a2), tuple(ccf_hat_a2), tuple(ccf_high_a2), tuple(ccf_low_a2), 'a2'))
                else:
                    logging.warning('Seg with inconsistent event: {}:{}:{}'.format(chrom, start, end))
            # neighbors = {s: set() for s in event_segs}
            # for seg1, seg2 in itertools.combinations(event_segs, 2):
            #     s1_hat = np.array(seg1[4])
            #     s2_hat = np.array(seg2[4])
            #     if seg1[2] == seg2[2] and seg1[3] == seg2[3] and all(s1_hat >= np.array(seg2[6])) and all(s1_hat <= np.array(seg2[5])) \
            #             and all(s2_hat >= np.array(seg1[6])) and all(s2_hat <= np.array(seg1[5])):
            #         neighbors[seg1].add(seg2)
            #         neighbors[seg2].add(seg1)

            # Non-pivoted Bron-Kerbosch algorithm enumerating all maximal
            # cliques in the undirected graph whose vertices are event_segs
            # and whose edges are defined by neighbors. This has exponential
            # time complexity. Use pivoted version instead!
            # def _BK(P, neighbors, R=frozenset(), X=frozenset()):
            #     if not P and not X:
            #         yield R
            #     else:
            #         for v in P:
            #             for r in _BK(P & neighbors[v], neighbors, R=R | {v}, X=X & neighbors[v]):
            #                 yield r
            #             P = P - {v}
            #             X = X | {v}

            def build_neighbors(segs):
                G = nx.Graph()
                G.add_nodes_from(segs)
                for s1, s2 in itertools.combinations(segs, 2):
                    s1_hat = np.array(s1[5])
                    s2_hat = np.array(s2[5])
                    if (    all(s1_hat >= np.array(s2[7])) and all(s1_hat <= np.array(s2[6]))
                        and all(s2_hat >= np.array(s1[7])) and all(s2_hat <= np.array(s1[6]))
                    ):
                        G.add_edge(s1, s2)
                return G

            buckets = defaultdict(list)
            for seg in event_segs:
                key = (seg[3], seg[4], seg[8])  # ('p'/'q', local_cn tuple, 'a1'/'a2')
                buckets[key].append(seg)

            # for clique in _BK(event_segs, neighbors):
            for key, segs in buckets.items():
                G = build_neighbors(segs)
                for clique in nx.find_cliques(G):
                    if clique:
                        clique_bands = set(clique[0][2])
                        clique_arm = clique[0][3]
                        arm_len = centromere if clique_arm == 'p' else csize - centromere
                        n_segs = len(clique)
                        clique_len = np.sum([seg[1] - seg[0] for seg in clique])
                        clique_ccf_hat = np.sum([np.array(seg[5]) for seg in clique], axis=0) / n_segs
                        clique_ccf_high = np.sum([np.array(seg[6]) for seg in clique], axis=0) / n_segs
                        clique_ccf_low = np.sum([np.array(seg[7]) for seg in clique], axis=0) / n_segs
                        clique_allele = clique[0][8]
                        local_cn = clique[0][4]
                        # cn_category = 'Arm_gain' if all(np.array(local_cn) > 1) else 'Arm_loss'
                        cn_category = self._infer_cn_category(np.array(local_cn), clique_ccf_hat)
                        if cn_category is None:
                            print("[DROPPED]", clique_bands, local_cn, clique_ccf_hat)
                            continue
                        cn_category = "Arm_" + cn_category
                        if clique_len < arm_len * .5:
                            continue
                        self._add_cn_event_to_samples(chrom, min(clique_bands), max(clique_bands), clique_arm, local_cn, cn_category, clique_ccf_hat, clique_ccf_high, clique_ccf_low, a1=clique_allele == "a1")

    def _add_cn_event_to_samples(self, chrom, start, end, arm, cns, cn_category, ccf_hat, ccf_high, ccf_low, a1=True, dupe=False):
        """
        Adds CN event to sample in low_coverage_mutations attr and mut hashtable

        """
        for i, sample in enumerate(self.sample_list):
            local_cn = cns[i]
            ccf_hat_i = ccf_hat[i] if local_cn != 1. else 0.
            ccf_high_i = ccf_high[i] if local_cn != 1. else 0.
            ccf_low_i = ccf_low[i] if local_cn != 1. else 0.
            cn = CopyNumberEvent(chrom, cn_category, start=start, end=end, ccf_hat=ccf_hat_i, ccf_high=ccf_high_i, ccf_low=ccf_low_i,
                                 local_cn=local_cn, from_sample=sample, arm=arm, a1=a1, dupe=dupe)
            sample.low_coverage_mutations.update({cn.var_str: cn})
            sample.add_muts_to_hashtable(cn)

    def get_bands(self, chrom, start, end):
        """
        Gets cytobands hit by a CN event
        """
        table = self.genome.cytoband_table
        return table.loc[
            table["chromosome"].isin([chrom, "chr" + chrom])
            & (start <= table["end"])
            & (table["start"] <= end)
            ].apply(lambda row: Cytoband(row["chromosome"], row["band"]), axis=1).to_list()


#################################################################################
# NDHistogram() - helper class to store combined histograms of mutations to pass to DP
#################################################################################
class NDHistogram:
    """CLASS info

    FUNCTIONS:

    PROPERTIES:
    """

    def __init__(self, hist_array, labels, phasing=None, ignore_nan=False):
        conv = 1e-40

        # We need to noramlize the array to 1. This is hard.
        # This step does the following:
        # convert to log space after adding a small convolution paramter so we don't get INF and NAN
        # for each mutation
        # normalize the row to 1 in logspace.

        hist = np.asarray(hist_array, dtype=np.float32) + conv
        if (~(hist > 0)).any():
            logging.error("Negative histogram bin or NAN mutation!")
            if ignore_nan:
                logging.warning("Caught ignore nan flag, not exiting, setting nan values to zero")
                hist[np.logical_not(hist > 0)] = conv
            else:
                sys.exit(1)

        n_samples = np.shape(hist)[1]
        for sample in range(n_samples):
            hist[:, :, 0] = conv  ##set zero bins
        self._hist_array = np.apply_over_axes(lambda x, y: np.apply_along_axis(lambda z: z - logsumexp(z), y, x),
                                              np.log(hist), 2)
        ####

        self._labels = labels
        self._label_ids = dict([[y, x] for x, y in enumerate(labels)])
        self._phasing = {} if phasing is None else phasing
        self.n_samples = np.shape(hist_array)[1]
        self.n_bins = np.shape(hist_array)[-1]

    def __getitem__(self, key):
        return self._hist_array[self._label_ids[key]]

    def phase(self, m1, m2):
        return self._phasing[frozenset({m1, m2})]

    def iteritems(self):
        for idx, mut in enumerate(self._hist_array):
            yield self._labels[idx], mut

    @property
    def mutations(self):
        return {}


# TODO: refactor Cytoband class to allow for other mapping of cytobands (other species)
# The cytoband names are the same between hg19 and hg38, so it doesn't matter here.

_cytoband_dict = {}
with open(os.path.dirname(__file__) + '/supplement_data/cytoBand.hg19.txt', 'r') as _f:
    for _i, _line in enumerate(_f):
        _row = _line.strip('\n').split('\t')
        _cytoband_dict[(_row[0], _row[3])] = _i

class Cytoband:
    band_nums = _cytoband_dict

    def __init__(self, chrom, band):
        self.chrom = chrom if chrom.startswith('chr') else 'chr' + chrom
        self.band = band

    def __sub__(self, other):
        if not isinstance(other, Cytoband):
            raise TypeError('Cannot subtract Cytoband with ' + str(type(other)))
        if self.chrom == other.chrom:
            self_num = self.band_nums[(self.chrom, self.band)]
            other_num = self.band_nums[(other.chrom, other.band)]
            return self_num - other_num
        raise ValueError('Cannot subtract cytobands on different chromosomes')

    def __lt__(self, other):
        if not isinstance(other, Cytoband):
            raise TypeError('Cannot compare Cytoband with ' + str(type(other)))
        if self.chrom == other.chrom:
            return self - other < 0
        raise ValueError('Cannot compare cytobands on different chromosomes')

    def __gt__(self, other):
        if not isinstance(other, Cytoband):
            raise TypeError('Cannot compare Cytoband with ' + str(type(other)))
        if self.chrom == other.chrom:
            return self - other > 0
        raise ValueError('Cannot compare cytobands on different chromosomes')

    def __le__(self, other):
        if not isinstance(other, Cytoband):
            raise TypeError('Cannot compare Cytoband with ' + str(type(other)))
        if self.chrom == other.chrom:
            return self - other <= 0
        raise ValueError('Cannot compare cytobands on different chromosomes')

    def __ge__(self, other):
        if not isinstance(other, Cytoband):
            raise TypeError('Cannot compare Cytoband with ' + str(type(other)))
        if self.chrom == other.chrom:
            return self - other >= 0
        raise ValueError('Cannot compare cytobands on different chromosomes')

    def __eq__(self, other):
        if not isinstance(other, Cytoband):
            raise TypeError('Cannot compare Cytoband with ' + str(type(other)))
        return hash(self) == hash(other)

    def __hash__(self):
        return hash((self.chrom, self.band))

    def __repr__(self):
        return str(self.chrom) + self.band

##############################################################
