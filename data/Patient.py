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
        self.hemizygous_chroms = set()  # chromosomes with one allele, filled for male patients by the Timing tool
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
            lines = [x.strip().split(" ")[0] for x in drv_file.read().strip().split('\n')]
            drv = [x for x in lines if not x.startswith("#") and len(x)]
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

    def _read_homdel_artifacts(self, path):
        """Regions where apparent homozygous deletions are known artifacts: chromosome, start, end, other columns ignored."""
        artifacts = defaultdict(list)
        if path:
            with open(path) as handle:
                for line in handle:
                    fields = line.split()
                    if line.startswith('#') or len(fields) < 3 or not (fields[1].isdigit() and fields[2].isdigit()):
                        continue
                    artifacts[fields[0].replace('chr', '')].append((int(fields[1]), int(fields[2])))
        return artifacts

    @staticmethod
    def _in_artifact(start, end, artifacts):
        """A segment is an artifact if the regions jointly cover half of it, or if it holds a region given as one position."""
        covered, reach = 0, start
        for s, e in sorted(artifacts):
            if s == e and start <= s < end:
                return True
            s, e = max(s, reach), min(e, end)
            if e > s:
                covered += e - s
                reach = e
        return covered >= .5 * (end - start)

    def _segment_homdel_call(self, sample, chrom, start, end, artifacts, max_length, max_copies, min_z, min_probes=5,
                             min_overlap=1000, flank_bp=1000000, min_sd=0.1):
        """
        Homozygous deletion in [start, end): one segment with so few copies that most tumor cells must hold none.

        With at most one copy left per cell, c copies mean that at least 1 - c of the cells have none. So max_copies 0.5
        asks for a majority of cells, the same rule as LoY at CCF above 50%. The value must also lie clearly below one
        copy given the segment's noise. The modal state of a whole window cannot do this: windows are megabases long, and
        modal 0/0 can come from rounding at low purity. CCF is the depth of the drop relative to the copy number flanking
        the segment. Returns (ccf_hat, ccf_high, ccf_low) of the deepest segment, or None if no segment qualifies.
        """
        tree = sample.CnProfile[chrom]
        slope = sample.copy_ratio_line()[0]
        if chrom in self.hemizygous_chroms:
            slope = slope / 2.  # one copy changes the copy ratio twice as much on a hemizygous chromosome
        best = None
        for iv in tree.overlap(start, end):
            seg = iv.data[1]
            if min(iv.end, end) - max(iv.begin, start) < min_overlap or iv.length() > max_length:
                continue
            copies = seg.get('total_cn_cont')
            if copies is None or (seg.get('n_probes') or 0) < min_probes or copies > max_copies:
                continue
            if self._in_artifact(iv.begin, iv.end, artifacts.get(chrom, ())):
                continue
            sd = (seg.get('seg_sigma') or 0.) * slope
            sd = max(sd, min_sd) if np.isfinite(sd) else min_sd
            copies = max(copies, 0.)
            if (1. - copies) / sd < min_z or (best is not None and copies >= best[0]):
                continue
            flanks = [(other.data[1].get('total_cn_cont'),
                       min(other.end, iv.end + flank_bp) - max(other.begin, iv.begin - flank_bp))
                      for other in tree.overlap(iv.begin - flank_bp, iv.end + flank_bp) if other != iv]
            flanks = [(value, length) for value, length in flanks if value is not None and value >= .5]
            flank_cn = 1.
            if flanks:
                order = sorted(flanks)
                halfway = .5 * sum(length for _, length in order)
                flank_cn = max(next(value for value, cumulative in zip(
                    (v for v, _ in order), itertools.accumulate(length for _, length in order)) if cumulative >= halfway), 1.)
            ccf = lambda value: float(np.clip(1. - max(value, 0.) / flank_cn, 0., 1.))
            best = (copies, ccf(copies), ccf(copies - 1.96 * sd), ccf(copies + 1.96 * sd))
        return None if best is None else best[1:]

    def _continuous_focal_call(self, sample, chrom, start, end, arm_start, arm_end, kind, artifacts, min_overlap, min_step=.5,
                               min_z=3., min_sd=.1, cn_tol=.2, flank_bp=2000000, max_arm_fraction=.3, min_probes=5, high_level_step=4.,
                               max_gain_length=None, high_level_ploidy=None, two_flanks=False, min_flank_fraction=.25,
                               max_gap_bp=100000):
        """
        A focal gain or loss in [start, end): the most extreme segment in the window, against what lies on either side.

        A focal event must stand out from both flanks. Against the arm alone, a normal segment next to a partial arm loss
        would read as gained, as CCND1 does beside a distal 11q deletion. High amplification is the exception: an amplicon
        can fill its own flanks, so a gain this far above the arm counts as well. The modal columns and cancer.cell.frac
        cannot show any of this, since modal states stay near 1/1 at high amplification and the fraction is the distance to
        the nearest integer for gains. Returns None, or a dict with the allele that moved most, its CCF, the state reached and
        the span of the segment measured; for a gain, extended over adjacent segments at the amplicon's level.

        With two_flanks, a gain that is not high-level needs two real flanks: segments covering min_flank_fraction of
        flank_bp on each side within the arm, so that neither side is the centromere gap or a chromosome end. Otherwise
        the remnant of an arm beside a partial loss passes on its one flank, as proximal 11q does. A high-level gain may
        fill its own flanks there too. Given max_gain_length, a gain's segment can be no longer. Given high_level_ploidy,
        a gain is high-level from that many times the ploidy above the arm, instead of high_level_step copies; losses
        keep high_level_step. The result says whether the gain reached it.
        """
        tree = sample.CnProfile[chrom]
        hemizygous = chrom in self.hemizygous_chroms
        gain_high_level_step = high_level_step
        if high_level_ploidy is not None:
            gain_high_level_step = high_level_ploidy * sample.ploidy_from_segments()
        window = []
        for iv in tree.overlap(start, end):
            if min(iv.end, end) - max(iv.begin, start) < min_overlap or iv.data[1].get('total_cn_cont') is None:
                continue
            if (iv.data[1].get('n_probes') or 0) < min_probes or iv.length() > max_arm_fraction * (arm_end - arm_start):
                continue
            if kind == 'AMP' and max_gain_length is not None and iv.length() > max_gain_length:
                continue
            if kind == 'DEL' and self._in_artifact(iv.begin, iv.end, artifacts.get(chrom, ())):
                continue
            window.append(iv)
        if not window:
            return None
        pick = max if kind == 'AMP' else min
        seg = pick(window, key=lambda iv: iv.data[1]['total_cn_cont'])

        def weighted_median(pairs):
            pairs = [(v, w) for v, w in pairs if v == v and w > 0]
            if not pairs:
                return float('nan')
            values, weights = (np.array(x, dtype=float) for x in zip(*sorted(pairs)))
            return float(values[np.searchsorted(np.cumsum(weights), .5 * weights.sum())])

        def side(lo, hi, cover=False):
            """Length-weighted medians of (total, a1, a2) over one flank, inside the arm, and the length they cover."""
            lo, hi = max(lo, arm_start), min(hi, arm_end)
            rows = [] if hi <= lo else [(iv.data[1].get('total_cn_cont'), sample.allele_copy_number(iv.data[1]), min(iv.end, hi) - max(iv.begin, lo))
                                        for iv in tree.overlap(lo, hi) if iv != seg]
            rows = [row for row in rows if row[0] is not None]
            medians = (weighted_median([(r[0], r[2]) for r in rows]), weighted_median([(r[1][0], r[2]) for r in rows]),
                       weighted_median([(r[1][1], r[2]) for r in rows]))
            return (medians, sum(r[2] for r in rows)) if cover else medians

        sides = [side(seg.begin - flank_bp, seg.begin, cover=True), side(seg.end, seg.end + flank_bp, cover=True)]
        high_level = False
        if kind == 'AMP':
            arm = side(arm_start, arm_end)
            high_level = seg.data[1]['total_cn_cont'] - arm[0] >= gain_high_level_step
        if kind == 'AMP' and two_flanks and not high_level and \
                any(length < min_flank_fraction * flank_bp for _, length in sides):
            return None
        flanks = [flank for flank, _ in sides if flank[0] == flank[0]]
        if not flanks:
            return None
        # The event must clear both sides: the higher flank for a gain and the lower one for a loss.
        flank = max(flanks, key=lambda f: f[0]) if kind == 'AMP' else min(flanks, key=lambda f: f[0])
        if high_level:
            flank = arm
        step = seg.data[1]['total_cn_cont'] - flank[0]
        if (kind == 'AMP') != (step > 0):
            return None
        slope = sample.copy_ratio_line()[0] / (2. if hemizygous else 1.)
        sd = (seg.data[1].get('seg_sigma') or 0.) * slope
        sd = max(sd, min_sd) if np.isfinite(sd) else min_sd
        # One allele has half the room to move, so a step counts from half the size.
        # Noise grows with copy ratio, so the test cannot apply to high amplification, whose size is not in doubt.
        exempt_step = gain_high_level_step if kind == 'AMP' else high_level_step
        if abs(step) < (min_step / 2. if hemizygous else min_step) or (abs(step) / sd < min_z and abs(step) < exempt_step):
            return None
        a1, a2 = sample.allele_copy_number(seg.data[1])
        delta1, delta2 = a1 - flank[1], a2 - flank[2]
        use_a1 = not hemizygous and abs(np.nan_to_num(delta1)) > abs(np.nan_to_num(delta2))
        delta = delta1 if use_a1 else delta2
        if not np.isfinite(delta) or abs(delta) < cn_tol or (delta > 0) != (kind == 'AMP'):
            return None
        reached = (lambda q: float(np.ceil(q - cn_tol))) if kind == 'AMP' else (lambda q: max(float(np.floor(q + cn_tol)), 0.))
        pair = (reached(a1), float(round(a2))) if use_a1 else (0. if hemizygous else float(round(a1)), reached(a2))
        if hemizygous and kind == 'DEL' and a2 >= cn_tol:
            pair = (0., float(a2))  # cells still hold the one allele, so this must not read as nothing left
        span = (int(seg.begin), int(seg.end))
        if kind == 'AMP':
            # An amplicon may be split into several segments. Its span runs over the neighbours, without a gap of more
            # than max_gap_bp, whose copy number lies closer to the measured segment than to the reference it cleared.
            level = flank[0] + .5 * step
            ordered = sorted((iv for iv in tree.overlap(arm_start, arm_end) if iv.data[1].get('total_cn_cont') is not None),
                             key=lambda iv: iv.begin)
            index = ordered.index(seg)
            lo, hi = seg.begin, seg.end
            for j in range(index - 1, -1, -1):
                if ordered[j].end < lo - max_gap_bp or ordered[j].data[1]['total_cn_cont'] < level:
                    break
                lo = min(lo, ordered[j].begin)
            for j in range(index + 1, len(ordered)):
                if ordered[j].begin > hi + max_gap_bp or ordered[j].data[1]['total_cn_cont'] < level:
                    break
                hi = max(hi, ordered[j].end)
            span = (int(lo), int(hi))
        return {'a1': use_a1, 'ccf': float(min(abs(delta), 1.)), 'cn': pair[0] if use_a1 else pair[1], 'pair': pair,
                'high_level': high_level, 'segment': span}

    def add_focal_cn_events(self, focal_regions, qval_max=None, focal_max_arm_frac=0.3, min_bp_overlap=1000, padding=100000, require_directional_consistency=False, split_across_centromere=True, segment_homdels=None, continuous=False,
                            focal_gain_max_length=3000000, focal_amp_ploidy=2.):
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
            segment_homdels: None keeps the window-mode homdel definition. A dict with artifact_file, max_length,
                max_copies and min_z calls homdels per segment instead, see _segment_homdel_call.
            continuous: call focal gains and losses from continuous copy number against the arm, see _continuous_focal_call.
                A gain then needs two real flanks and a segment of at most focal_gain_max_length bp, and is high-level
                at focal_amp_ploidy times the ploidy above the arm.
        """
        homdel_artifacts = self._read_homdel_artifacts(segment_homdels['artifact_file']) if segment_homdels else {}

        # --- helpers ---
        def _normalize_focals(focals):
            """Return dict: chrom -> list of dicts {chrom,start,end,type,q,name}."""
            out = defaultdict(list)
            if isinstance(focals, str):
                with open(focals) as fin:
                    for ln in fin:
                        if not ln.strip() or ln.startswith('#') or ln.startswith("track name") or ln.startswith("chromosome"):
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
            if allele_key == 'a1' and chrom in self.hemizygous_chroms:
                return 0, 0.0, 0.0, 0.0  # no second allele to lose: its absence is the normal state

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

                    if continuous:
                        n_samp = len(self.sample_list)
                        is_deletion = R.get('type') in ('DEL', 'DELETION', None)
                        calls = [None] * n_samp
                        if segment_homdels is not None and is_deletion:
                            calls = [self._segment_homdel_call(sample, chrom, s0, e0, homdel_artifacts, segment_homdels['max_length'],
                                                               segment_homdels['max_copies'], segment_homdels['min_z'])
                                     if sample.CnProfile else None for sample in self.sample_list]
                            if any(call is not None for call in calls):
                                none = (0., 0., 0.)
                                self._add_cn_event_to_samples(chrom=chrom, start=s0, end=e0, arm=arm_lbl,
                                                              cns=tuple(1 if call is None else 0 for call in calls),
                                                              cn_category='Focal_homdel',
                                                              ccf_hat=tuple((call or none)[0] for call in calls),
                                                              ccf_high=tuple((call or none)[1] for call in calls),
                                                              ccf_low=tuple((call or none)[2] for call in calls),
                                                              a1=True, dupe=False)
                        kind = 'DEL' if is_deletion else 'AMP'
                        arm_start, arm_end = (0, centromere) if arm_lbl == 'p' else (centromere, csize)
                        focal = [None if calls[i] is not None or not sample.CnProfile else
                                 self._continuous_focal_call(sample, chrom, s0, e0, arm_start, arm_end, kind, homdel_artifacts, min_bp_overlap,
                                                             max_gain_length=focal_gain_max_length, high_level_ploidy=focal_amp_ploidy,
                                                             two_flanks=True)
                                 for i, sample in enumerate(self.sample_list)]
                        carriers = [call for call in focal if call is not None]
                        if carriers:
                            is_a1 = carriers[0]['a1']
                            ccfs = tuple(0. if call is None else call['ccf'] for call in focal)
                            self._add_cn_event_to_samples(chrom=chrom, start=s0, end=e0, arm=arm_lbl,
                                                          cns=tuple(1 if call is None else call['cn'] for call in focal),
                                                          cn_category='Focal_gain' if kind == 'AMP' else 'Focal_loss',
                                                          ccf_hat=ccfs, ccf_high=ccfs, ccf_low=ccfs, a1=is_a1,
                                                          dupe=(kind == 'AMP') != is_a1,
                                                          state_pairs=tuple(None if call is None else call['pair'] for call in focal),
                                                          high_levels=tuple(None if call is None else call['high_level'] for call in focal),
                                                          segments=tuple(None if call is None else call['segment'] for call in focal))
                        continue

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

                    # ----------------------------------------------------------
                    # Detect homozygous deletion: at least one sample has both
                    # alleles at CN=0 in this focal region. Per-sample CCFs are
                    # preserved, so other samples (neutral, hemizygous loss,
                    # gain) are represented faithfully in the per-sample event.
                    # ----------------------------------------------------------
                    n_samp = len(self.sample_list)
                    cn_a1, hat_a1, hi_a1, lo_a1 = results.get('a1', ((1,)*n_samp, (0.,)*n_samp, (0.,)*n_samp, (0.,)*n_samp))
                    cn_a2, hat_a2, hi_a2, lo_a2 = results.get('a2', ((1,)*n_samp, (0.,)*n_samp, (0.,)*n_samp, (0.,)*n_samp))

                    # Any sample with CN=0 on BOTH alleles → homdel carrier
                    homdel_carriers = {i for i in range(n_samp) if cn_a1[i] == 0 and cn_a2[i] == 0}
                    is_homdel = len(homdel_carriers) > 0
                    hemizygous = chrom in self.hemizygous_chroms

                    if segment_homdels is not None:
                        is_homdel = False
                        calls = [None] * n_samp
                        if R.get('type') in ('DEL', 'DELETION', None):
                            calls = [self._segment_homdel_call(sample, chrom, s0, e0, homdel_artifacts,
                                                               segment_homdels['max_length'], segment_homdels['min_posterior'])
                                     if sample.CnProfile else None for sample in self.sample_list]
                        if any(call is not None for call in calls):
                            none = (0., 0., 0.)
                            self._add_cn_event_to_samples(chrom=chrom, start=s0, end=e0, arm=arm_lbl,
                                                          cns=tuple(1 if call is None else 0 for call in calls),
                                                          cn_category='Focal_homdel',
                                                          ccf_hat=tuple((call or none)[0] for call in calls),
                                                          ccf_high=tuple((call or none)[1] for call in calls),
                                                          ccf_low=tuple((call or none)[2] for call in calls),
                                                          a1=True, dupe=False)

                    if is_homdel and R.get('type') in ('DEL', 'DELETION', None):
                        # Per-sample CN: 0 for carriers, 1 (neutral) for non-carriers.
                        # Carriers contribute their homdel CCF; non-carriers get 0.
                        homdel_cns = tuple(0 if i in homdel_carriers else 1 for i in range(n_samp))
                        homdel_hat = tuple(max(hat_a1[i], hat_a2[i]) if i in homdel_carriers else 0.
                                          for i in range(n_samp))
                        homdel_hi = tuple(max(hi_a1[i], hi_a2[i]) if i in homdel_carriers else 0.
                                         for i in range(n_samp))
                        homdel_lo = tuple(max(lo_a1[i], lo_a2[i]) if i in homdel_carriers else 0.
                                         for i in range(n_samp))
                        self._add_cn_event_to_samples(chrom=chrom, start=s0, end=e0, arm=arm_lbl,
                                                     cns=homdel_cns, cn_category='Focal_homdel',
                                                     ccf_hat=homdel_hat, ccf_high=homdel_hi, ccf_low=homdel_lo,
                                                     a1=True, dupe=False)
                        # Note: do NOT `continue` here — also let the per-allele
                        # loss/gain emission below run for non-carrier samples.
                        # The Focal_homdel and Focal_loss events are independent.

                    # emit an event for each allele with any non-neutral evidence
                    for allele_key, (cn_vec, hat, hi, lo) in results.items():
                        if all(cn == 1 for cn in cn_vec):
                            continue
                        if hemizygous and allele_key == 'a1':
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
                            a1 = allele_key == 'a1'
                            a2 = not a1
                            is_gain = "gain" in cn_category
                            is_loss = not is_gain
                            dupe = is_gain and a2 or is_loss and a1
                            self._add_cn_event_to_samples(chrom=chrom, start=s0, end=e0, arm=arm_lbl, cns=cn_vec, cn_category=cn_category, ccf_hat=hat, ccf_high=hi, ccf_low=lo, a1=a1, dupe=dupe)

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

    def _add_cn_event_to_samples(self, chrom, start, end, arm, cns, cn_category, ccf_hat, ccf_high, ccf_low, a1=True, dupe=False, state_pairs=None,
                                 high_levels=None, segments=None):
        """
        Adds CN event to sample in low_coverage_mutations attr and mut hashtable

        """
        for i, sample in enumerate(self.sample_list):
            local_cn = cns[i]
            # With a state given, one copy can itself be an event: a step down from two copies on the arm.
            carries = local_cn != 1. if state_pairs is None else state_pairs[i] is not None
            ccf_hat_i = ccf_hat[i] if carries else 0.
            ccf_high_i = ccf_high[i] if carries else 0.
            ccf_low_i = ccf_low[i] if carries else 0.
            cn = CopyNumberEvent(chrom, cn_category, start=start, end=end, ccf_hat=ccf_hat_i, ccf_high=ccf_high_i, ccf_low=ccf_low_i,
                                 local_cn=local_cn, from_sample=sample, arm=arm, a1=a1, dupe=dupe)
            cn.state_pair = None if state_pairs is None else state_pairs[i]
            cn.high_level = None if high_levels is None else high_levels[i]
            cn.segment = None if segments is None else segments[i]  # (start, end) of the segment the call measured
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
