from collections import Counter, defaultdict
import numpy as np
import itertools
from intervaltree import IntervalTree, Interval
import scipy.stats
import scipy.interpolate
from scipy.special import logsumexp
import random
import zlib
import data.Enums as Enums

# Mappable chrY blocks in hg19. Neither Y arm reaches the 50% coverage required of an arm, so chrY is called as one
# event over this territory. Short blocks outside it keep signal even in tumors that lost Y, likely from cross-mapping.
Y_CALLABLE_BLOCKS_HG19 = ((6420000, 10070000), (13560000, 19560000), (20800000, 23810000))
SEX_LABELS = {'XY': 'male', 'MALE': 'male', 'M': 'male', 'XX': 'female', 'FEMALE': 'female', 'F': 'female'}
# Smallest uncertainty assumed for an event's CCF, because purity error moves copy number beyond segment noise.
MIN_CCF_SD = .05
Y_COPY_RATIO_FLOOR = .08
# Key suffix of a terminal partial-arm state, kept beside the whole-arm state of the same arm.
TERMINAL = ':terminal'


class TimingOptions(object):
    """
    Opt-in timing behavior. The defaults reproduce the results of the code without these options.

    sex_aware: male X and Y are hemizygous (normal copy number 0/1, or 0/2 after WGD), and chrY gets a whole-chromosome
        event from callable territory. Requires sex.
    continuous_cn: arm events and their CCF come from continuous copy number, so subclonal arm events exist as such.
    clonal_unresolved: clonal events without timing evidence precede subclonal events and stay unknown against clonal ones.
    segment_homdels: focal homozygous deletions come from Patient.add_focal_cn_events(segment_homdels=...), are present
        in WGD tumors, and are timed by the WGD rule or as clonal_unresolved instead of with a uniform distribution.
    cn_event_tol: smallest deviation from the normal copy number that counts as an event when continuous_cn is on,
        and the largest copy number still read as no copy left.
    cn_clonal_min_ccf: copy-number events below this CCF in any sample are subclonal, also when the tumor has no
        subclonal mutation cluster to assign them to. Events whose CCF cannot be told from it stay unknown.
    wgd_override: {sample name: True or False} for samples whose doubled status is known from elsewhere. It replaces
        the arm rule of TimingSample.call_wgd for those samples; samples not in it keep the rule.
    loy_wgd_min_odds: None, or the odds of a loss of Y before the WGD against two losses after it that every sample
        must reach, in loy_odds {sample name: odds}, for a doubled tumor without Y to be timed before its WGD. Below
        it, or without odds, the loss is clonal and unresolved. Parsimony alone is weak when the WGD came early.
    seed: None uses the global random state. A number gives each simulation of detection power its own random state,
        derived from the seed and the event, so a run is reproducible and does not depend on the order events are timed.
    terminal_arm_min_fraction: None, or the share of an arm a terminal run of segments must cover to be called as an
        event of that arm when the whole-arm call misses it, as for 1p36 and distal 11q losses. Needs continuous_cn.
    terminal_arm_min_length: the run must also span this many bp, so that on short arms it is no focal event.
    whole_chromosome_events: an event found on both arms of a chromosome in the same state is one whole-chromosome
        event (gain_17, not gain_17p and gain_17q), so numerical and segmental changes stay apart.
    """
    def __init__(self, sex=None, sex_aware=False, continuous_cn=False, clonal_unresolved=False, segment_homdels=False,
                 y_callable_blocks=Y_CALLABLE_BLOCKS_HG19, cn_event_tol=.2, cn_clonal_min_ccf=.8, wgd_override=None,
                 loy_wgd_min_odds=None, loy_odds=None, seed=None, terminal_arm_min_fraction=None,
                 terminal_arm_min_length=5000000, whole_chromosome_events=False):
        label = SEX_LABELS.get(str(sex).strip().upper())
        if sex_aware and label is None:
            raise ValueError('sex_aware timing requires the sex of the patient, as XY or XX. Got: {!r}'.format(sex))
        self.sex = label
        self.sex_aware = sex_aware
        self.continuous_cn = continuous_cn
        self.clonal_unresolved = clonal_unresolved
        self.segment_homdels = segment_homdels
        self.y_callable_blocks = y_callable_blocks
        self.cn_event_tol = cn_event_tol
        self.cn_clonal_min_ccf = cn_clonal_min_ccf
        self.wgd_override = dict(wgd_override or {})
        self.loy_wgd_min_odds = loy_wgd_min_odds
        self.loy_odds = dict(loy_odds or {})
        self.seed = seed
        if terminal_arm_min_fraction is not None and not continuous_cn:
            raise ValueError('Terminal partial-arm events are called from continuous copy number: set continuous_cn.')
        self.terminal_arm_min_fraction = terminal_arm_min_fraction
        self.terminal_arm_min_length = terminal_arm_min_length
        self.whole_chromosome_events = whole_chromosome_events

    @property
    def male(self):
        return self.sex_aware and self.sex == 'male'

    def is_hemizygous(self, chrN):
        return self.male and str(chrN).replace('chr', '') in ('X', 'Y')

    @property
    def extended_output(self):
        return self.sex_aware or self.continuous_cn or self.clonal_unresolved or self.segment_homdels


def _number(value):
    return np.nan if value is None else float(value)


def _weighted_median(values, weights):
    values, weights = np.asarray(values, dtype=float), np.asarray(weights, dtype=float)
    keep = np.isfinite(values) & (weights > 0)
    if not keep.any():
        return np.nan
    values, weights = values[keep], weights[keep]
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cumulative, 0.5 * cumulative[-1])])


class TimingEngine(object):
    """
    Class for holding samples
    """
    def __init__(self, patient, cn_state_whitelist=Enums.cn_state_whitelist, min_supporting_muts=3, options=None):
        self.options = options if options is not None else TimingOptions()
        self.patient = patient
        self.genome = patient.genome
        self.cn_state_whitelist = cn_state_whitelist
        self.arm_regions = list(itertools.product(self.genome.CHROMS, self.genome.ARMS))
        self.borrow_arm_prior = True
        self.min_supporting_muts = min_supporting_muts
        self.sample_list = []
        for sample in self.patient.sample_list:
            timing_sample = TimingSample(sample, self)
            self.sample_list.append(timing_sample)
        self.concordant_cn_states = {}
        self.timeable_muts = {}
        self.all_cn_events = {}
        self.get_concordant_cn_states()
        self.WGD = None
        self.call_wgd()
        self.mutations = {}
        self.get_mutations()
        self.truncal_cn_events = {}
        self.get_arm_level_cn_events()
        if self.options.whole_chromosome_events:
            self.merge_whole_chromosome_events()
        self.get_focal_cn_events()

    def _clonality(self, ccf, sd=0.):
        """True if the CCF is clearly clonal, False if clearly not, None if it cannot be told from the cutoff."""
        margin = 1.96 * max(sd, MIN_CCF_SD)
        if ccf + margin < self.options.cn_clonal_min_ccf:
            return False
        if ccf - margin < self.options.cn_clonal_min_ccf:
            return None
        return True

    def _loy_before_wgd_unsupported(self, event):
        """True if the WGD rule would place this loss of Y before the WGD while the odds given do not support it."""
        if self.options.loy_wgd_min_odds is None or event.Type != 'Chrom_loss' or str(event.chrN) != 'Y' \
                or event.allelic_cn != 0.:
            return False
        odds = [self.options.loy_odds.get(sample.sample_name) for sample in self.sample_list]
        return not all(value is not None and value >= self.options.loy_wgd_min_odds for value in odds)

    def _drop_wgd_only_states(self):
        """Arm states kept only to time a WGD have no events. Without a WGD they would read as gains against one copy."""
        if not self.options.extended_output:
            return
        for region in list(self.concordant_cn_states):
            if not any(region in s.cn_states and len(s.cn_states[region].cn_events) for s in self.sample_list):
                del self.concordant_cn_states[region]

    def get_concordant_cn_states(self):
        """
        get cn states that are present in all samples (same chromsome arm and copy number) and create multisample cn_state instances
        """
        n_samples = len(self.sample_list)
        purity = [sample.purity for sample in self.sample_list]
        self.concordant_cn_states = {}
        regions = [(chrN, arm, chrN + arm) for chrN, arm in self.arm_regions] + ([('Y', '', 'Y')] if self.options.male else [])
        if self.options.terminal_arm_min_fraction is not None:
            regions += [(chrN, arm, chrN + arm + TERMINAL) for chrN, arm in self.arm_regions]
        for chrN, arm, region in regions:
            states_across_samples = {}  # hold union of all copy number states in all samples
            for sample in self.sample_list:
                # After WGD a 2/2 arm has no events, but it is the main evidence for timing the WGD itself.
                keep_for_wgd = (self.options.extended_output and sample.WGD is not None and region in sample.cn_states
                                and sample.cn_states[region].cn_a2 >= 2)
                if region in sample.cn_states and (len(sample.cn_states[region].cn_events) or keep_for_wgd):
                    # Coerce copy number states to integer states
                    if 0. <= sample.cn_states[region].cn_a1 < 1.:
                        cn_a1 = 0.
                    elif 1. < sample.cn_states[region].cn_a1 <= 2.:
                        cn_a1 = 2.
                    else:
                        cn_a1 = round(sample.cn_states[region].cn_a1)
                    if 0. <= sample.cn_states[region].cn_a2 < 1.:
                        cn_a2 = 0.
                    elif 1. < sample.cn_states[region].cn_a2 <= 2.:
                        cn_a2 = 2.
                    else:
                        cn_a2 = round(sample.cn_states[region].cn_a2)
                    state = (cn_a1, cn_a2)
                    states_across_samples.setdefault(state, [])
                    states_across_samples[state].append(sample.cn_states[region])
            for state in states_across_samples:
                if len(states_across_samples[state]) == n_samples:  # filter states not present in all samples
                    # get all supporting muts from all samples
                    supporting_mut_varstr = {mut.var_str for mut in states_across_samples[state][0].supporting_muts}
                    for sample_state in states_across_samples[state]:
                        supporting_mut_varstr = supporting_mut_varstr & {mut.var_str for mut in sample_state.supporting_muts}
                    supporting_muts = []
                    for mut_varstr in supporting_mut_varstr:
                        mut = self.sample_list[0][mut_varstr]
                        alt_cnt = np.zeros(n_samples)
                        ref_cnt = np.zeros(n_samples)
                        ccf_dist = np.zeros((n_samples, 101))
                        cn_a1 = np.zeros(n_samples)
                        cn_a2 = np.zeros(n_samples)
                        for i, sample in enumerate(self.sample_list):
                            alt_cnt[i] = sample[mut_varstr].alt_cnt
                            ref_cnt[i] = sample[mut_varstr].ref_cnt
                            ccf_dist[i, :] = sample[mut_varstr].ccf_dist
                            cn_a1[i] = sample[mut_varstr].local_cn_a1
                            cn_a2[i] = sample[mut_varstr].local_cn_a2
                        nd_mut = TimingMut(self.sample_list, mut.gene, mut.chrN, mut.pos, mut.alt, mut.ref, alt_cnt,
                                           ref_cnt, cn_a1, cn_a2, ccf_dist, prot_change=mut.prot_change,
                                           cluster_assignment=mut.cluster_assignment)
                        supporting_muts.append(nd_mut)
                        self.timeable_muts[mut_varstr] = nd_mut
                    # The patient-level event must exist in every sample: an event private to one sample is dropped,
                    # as for integer states, so keep the value closest to the shared state. CCF stays per sample.
                    continuous, continuous_sd = None, None
                    sample_states = [s for s in states_across_samples[state] if s.continuous is not None]
                    if sample_states:
                        def closest(i):
                            finite = [s for s in sample_states if np.isfinite(s.continuous[i])]
                            return min(finite, key=lambda s: abs(s.continuous[i] - state[i])) if finite else None
                        chosen = [closest(0), closest(1)]
                        continuous = tuple(np.nan if c is None else c.continuous[i] for i, c in enumerate(chosen))
                        continuous_sd = tuple(0. if c is None or c.continuous_sd is None else c.continuous_sd[i]
                                              for i, c in enumerate(chosen))
                    multisample_state = TimingCNState(self.sample_list, chrN, arm, state, purity,
                                                      supporting_muts=supporting_muts, continuous=continuous,
                                                      continuous_sd=continuous_sd,
                                                      cn_type_prefix=states_across_samples[state][0].cn_type_prefix)
                    limits = [s.absent_tol for s in states_across_samples[state] if s.absent_tol is not None]
                    if limits:
                        multisample_state.absent_tol = max(limits)
                        multisample_state.call_events()
                    multisample_state.terminal = states_across_samples[state][0].terminal
                    self.concordant_cn_states[region] = multisample_state

    def call_wgd(self, concordance_threshold=.8, tol=.2):  # TODO: try other concordance thresholds
        """
        call WGD event for a patient across all samples
        """
        min_pi = np.ones(101)
        for sample in self.sample_list:  # call WGD event in samples first
            sample.fill_mutation_intervaltree(concordant_muts=self.timeable_muts)
            sample.call_arm_level_cn_states(concordant_states=self.concordant_cn_states)
            sample.call_wgd(use_concordant_states=True)
            sample.get_arm_level_cn_events(use_concordant_states=True)
            if sample.concordant_WGD is not None:
                sample.concordant_WGD.get_pi_dist()
                min_pi = np.minimum(min_pi, sample.concordant_WGD.pi_dist)  # only call WGD in patient if sample-level WGD events have similar pi distributions
        if not all([sample.concordant_WGD is not None for sample in self.sample_list]) or sum(min_pi) < concordance_threshold:
            self._drop_wgd_only_states()
            return  # if not samples have WGD event or concordance is below threshold do not call WGD
        regions_supporting_WGD = []
        for chrN, arm in self.arm_regions:
            region = str(chrN) + arm
            if region in self.concordant_cn_states:
                cn_state = self.concordant_cn_states[region]
                if all(sample.cn_states[region].cn_a2 >= 2 for sample in self.sample_list):
                    # make new TimingCNState instances for WGD since patient level states will be called relative to WGD
                    regions_supporting_WGD.append(TimingCNState([self], cn_state.chrN, cn_state.arm,
                        (cn_state.cn_a1, cn_state.cn_a2), cn_state.purity, supporting_muts=cn_state.supporting_muts))
        self.WGD = TimingWGD(supporting_arm_states=regions_supporting_WGD)
        for cn_state in self.concordant_cn_states.values():
            cn_state.call_events(wgd=True)  # call copy number events relative to WGD

    def get_mutations(self):
        for mut_varstr in self.sample_list[0].mut_lookup_table:
            if mut_varstr in self.timeable_muts:
                self.mutations[mut_varstr] = self.timeable_muts[mut_varstr]
            else:
                mut = self.sample_list[0][mut_varstr]
                n_samples = len(self.sample_list)
                alt_cnt = np.zeros(n_samples)
                ref_cnt = np.zeros(n_samples)
                ccf_dist = np.zeros((n_samples, 101))
                cn_a1 = np.zeros(n_samples)
                cn_a2 = np.zeros(n_samples)
                for i, sample in enumerate(self.sample_list):
                    alt_cnt[i] = sample[mut_varstr].alt_cnt
                    ref_cnt[i] = sample[mut_varstr].ref_cnt
                    ccf_dist[i, :] = sample[mut_varstr].ccf_dist
                    cn_a1[i] = sample[mut_varstr].local_cn_a1
                    cn_a2[i] = sample[mut_varstr].local_cn_a2
                nd_mut = TimingMut(self.sample_list, mut.gene, mut.chrN, mut.pos, mut.alt, mut.ref, alt_cnt,
                                   ref_cnt, cn_a1, cn_a2, ccf_dist, prot_change=mut.prot_change,
                                   cluster_assignment=mut.cluster_assignment)
                self.mutations[mut_varstr] = nd_mut

    def get_arm_level_cn_events(self):
        """
        attach truncal copy number events to engine
        """
        self.truncal_cn_events = {gl + chrN + arm: [] for gl, (chrN, arm) in itertools.product(('gain_', 'loss_'), self.arm_regions)}
        self.all_cn_events = {gl + chrN + arm: [] for gl, (chrN, arm) in itertools.product(('gain_', 'loss_'), self.arm_regions)}
        cluster_ccfs = self._get_cluster_ccfs()
        sample_idx = range(len(self.sample_list))
        for cn_state, state in self.concordant_cn_states.items():
            # iterate over the *reference* events (typically Arm_gain and Arm_loss)
            for eve in state.cn_events:
                etype = eve.Type  # 'Arm_gain' or 'Arm_loss'
                ccf_hats = []
                ccf_sds = []
                ccf_tests = []
                for sample in self.sample_list:
                    s_state = sample.cn_states.get(cn_state)
                    if s_state is None:
                        ccf_hats.append(0.0)
                        ccf_sds.append(0.0)
                        ccf_tests.append(0.0)
                        continue
                    # find the matching event by Type and allele; default to 0 if missing
                    match = next((e for e in s_state.cn_events if e.Type == etype and e.allele == eve.allele), None)
                    ccf_hats.append(float(getattr(match, 'ccf_hat', 0.0) or 0.0))
                    ccf_sds.append(float(getattr(match, 'ccf_sd', 0.0) or 0.0))
                    test = getattr(match, 'clonality_ccf', None)
                    ccf_tests.append(ccf_hats[-1] if test is None else float(test))
                ccf_idx = np.clip(np.round(np.array(ccf_hats) * 100).astype(int), 0, 100)
                clonal_concordance = np.prod(cluster_ccfs[1][sample_idx, ccf_idx])
                is_clonal = True
                for c in cluster_ccfs:
                    if np.prod(cluster_ccfs[c][sample_idx, ccf_idx]) > clonal_concordance:
                        is_clonal = False
                        break
                if self.options.continuous_cn or etype.startswith('Chrom_'):
                    # CCF from copy number decides. Mutation clusters are too narrow to judge a CCF of 0.9 against 1.
                    lowest = int(np.argmin(ccf_tests))
                    called = self._clonality(ccf_tests[lowest], ccf_sds[lowest])
                    is_clonal = called is True
                    eve.clonality_ambiguous = called is None
                    eve.ccf_hat = ccf_hats[lowest]
                eve.is_clonal = is_clonal
                self.all_cn_events.setdefault(eve.event_name, []).append(eve)
                if is_clonal:
                    self.truncal_cn_events.setdefault(eve.event_name, []).append(eve)

    def get_focal_cn_events(self, cn_tol=0.20):
        """
        Attach focal CNV events (produced by Patient.add_focal_cn_events) as TimingCNEvent instances.

        Design:
          - build focal events from LOCAL interval mutations only
          - build multi-sample (ND) TimingMut objects across carrier samples
          - do NOT borrow same-arm mutations as direct evidence
          - optionally attach metadata so an arm-level prior can be used later in time_events()

        Notes:
          - This is intentionally closer to the arm-level timing model than the previous focal implementation.
          - For focal events with too few local supporting mutations, we leave pi_dist unset here and let
            time_events() decide whether to use uniform, "not estimable", or an optional borrowed prior.
        """
        n_samples = len(self.sample_list)
        if n_samples == 0:
            return

        # ------------------------------------------------------------------
        # cytoband helpers
        # ------------------------------------------------------------------
        _band_index = {}
        _bands_by_chr = {}
        _band_index_ready = False

        def _ensure_band_index():
            nonlocal _band_index_ready
            if _band_index_ready:
                return
            try:
                for _, row in self.genome.cytoband_table.iterrows():
                    chrom, start, end, band = row.iloc[:4]
                    chrom = str(chrom).lstrip('chr')
                    _band_index[(chrom, str(band))] = (int(start), int(end))
                tmp = defaultdict(list)
                for (chrom, band), (start, end) in _band_index.items():
                    tmp[chrom].append((start, end, band))
                for chrom in tmp:
                    _bands_by_chr[chrom] = sorted(tmp[chrom], key=lambda t: t[0])
            except Exception as e:
                print(f'Could not load cytoBand.txt for band to bp conversion: {e}')
            _band_index_ready = True

        def _band_to_bp(chrom, band_obj):
            if band_obj is None:
                return None
            _ensure_band_index()
            name = getattr(band_obj, 'band', None) or str(band_obj)
            return _band_index.get((str(chrom).lstrip('chr'), name))

        def _band_at_bp(chrom, pos):
            _ensure_band_index()
            arr = _bands_by_chr.get(str(chrom).lstrip('chr'), [])
            for start, end, band in arr:
                if start <= pos < end:
                    return band
            return None

        def _band_range_str(chrom, start_bp, end_bp):
            b1 = _band_at_bp(chrom, start_bp)
            b2 = _band_at_bp(chrom, max(start_bp, end_bp - 1))
            if not b1 or not b2:
                return None
            return b1 if b1 == b2 else f'{b1}-{b2}'

        def _coords(chrom, start, end):
            """
            Normalize event bounds to integer bp coordinates.
            Accepts ints or cytoband-like objects/strings.
            """
            if isinstance(start, (int, np.integer)) and isinstance(end, (int, np.integer)):
                return int(start), int(end)
            s = _band_to_bp(chrom, start)
            e = _band_to_bp(chrom, end)
            if s is None or e is None:
                return None
            return s[0], e[1]

        # ------------------------------------------------------------------
        # CN coercion / matching helpers
        # ------------------------------------------------------------------
        def _scalar(x):
            arr = np.asarray(x).reshape(-1)
            if arr.size == 0:
                return np.nan
            return float(arr[0])

        def _coerce_cn(v):
            if v is None:
                return None
            v = _scalar(v)
            if np.isnan(v):
                return None
            r = round(v)
            return float(r) if abs(v - r) <= cn_tol else float(v)

        def _coerce_pair(a1, a2):
            a1c = _coerce_cn(a1)
            a2c = _coerce_cn(a2)
            if a1c is None or a2c is None:
                return None
            return (a1c, a2c)

        def _best_segment_pair(ts, chrom, start_bp, end_bp, category):
            """
            Pick the best overlapping segment for the focal window in one sample,
            then return its coerced (cn_a1, cn_a2).
            """
            try:
                ov = list(ts.CnProfile[str(chrom)].overlap(start_bp, end_bp))
            except Exception:
                ov = []
            if not ov:
                return None

            def _ovlen(iv):
                return max(0, min(iv.end, end_bp) - max(iv.begin, start_bp))

            def _totcn(iv):
                seg = iv.data[1]
                return float(seg['cn_a1']) + float(seg['cn_a2'])

            is_gain = ('gain' in category.lower()) or ('amp' in category.lower())
            if is_gain:
                best = max(ov, key=lambda iv: (_totcn(iv), _ovlen(iv)))
            else:
                best = min(ov, key=lambda iv: (_totcn(iv), -_ovlen(iv)))

            seg = best.data[1]
            return _coerce_pair(seg['cn_a1'], seg['cn_a2'])

        def _interval_muts_matching_state(ts, chrom, start_bp, end_bp, state):
            """
            Return {var_str -> sample-level TimingMut} for LOCAL mutations in the focal interval
            whose local CN state exactly matches the focal state after coercion.
            """
            out = {}
            tree = ts.mutation_intervaltree.get(str(chrom), IntervalTree())
            for iv in tree[start_bp:end_bp]:
                mut = iv.data  # sample-level TimingMut
                mpair = _coerce_pair(mut.local_cn_a1, mut.local_cn_a2)
                if mpair == state:
                    out[mut.var_str] = mut
            return out

        # ------------------------------------------------------------------
        # group focal events across samples
        # key = (chrom, start, end, arm, is_a1, category)
        # ------------------------------------------------------------------
        groups = {}
        for si, ts in enumerate(self.sample_list):
            for cn in ts.sample.low_coverage_mutations.values():
                if getattr(cn, 'type', None) != 'CNV':
                    continue
                if not getattr(cn, 'cn_category', '').startswith('Focal_'):
                    continue

                bounds = _coords(cn.chrN, cn.start, cn.end)
                if bounds is None:
                    print(f'Skipping focal CNV without resolvable coords: {cn}')
                    continue

                start_bp, end_bp = bounds
                key = (str(cn.chrN), start_bp, end_bp, cn.arm, bool(cn.a1), cn.cn_category)
                rec = groups.setdefault(
                    key,
                    {
                        'ccf': np.zeros(n_samples, dtype=float),
                        'pairs': [None] * n_samples,
                        'high_level': [None] * n_samples,
                    }
                )
                rec['ccf'][si] = float(getattr(cn, 'ccf_hat', 0.0))
                rec['high_level'][si] = getattr(cn, 'high_level', None)
                if self.options.segment_homdels and cn.cn_category == 'Focal_homdel':
                    # The call was made on one segment. The lowest modal state in the window need not agree with it.
                    rec['pairs'][si] = (0., 0.)
                elif self.options.continuous_cn:
                    rec['pairs'][si] = getattr(cn, 'state_pair', None)
                else:
                    rec['pairs'][si] = _best_segment_pair(ts, cn.chrN, start_bp, end_bp, cn.cn_category)

        if not groups:
            return

        cluster_ccfs = self._get_cluster_ccfs()
        sample_idx = np.arange(n_samples)

        # ------------------------------------------------------------------
        # build one TimingCNEvent per grouped focal event
        # ------------------------------------------------------------------
        for (chrom, start_bp, end_bp, arm, is_a1, category), rec in groups.items():

            # carrier samples = samples where the focal event exists and a CN state is resolvable
            carrier_idx = [
                i for i in range(n_samples)
                if rec['pairs'][i] is not None and rec['ccf'][i] > 0.0
            ]
            if not carrier_idx:
                continue

            # consensus state among carriers
            pair_counts = Counter(rec['pairs'][i] for i in carrier_idx)
            consensus_state, _ = pair_counts.most_common(1)[0]

            # keep only carriers whose focal state matches the consensus exactly
            matched_carrier_idx = [i for i in carrier_idx if rec['pairs'][i] == consensus_state]
            if not matched_carrier_idx:
                continue

            # local, event-specific supporting muts:
            # intersection across carrier samples, analogous to arm-level concordant states
            per_sample_local = {}
            for si in matched_carrier_idx:
                ts = self.sample_list[si]
                per_sample_local[si] = _interval_muts_matching_state(
                    ts, chrom, start_bp, end_bp, consensus_state
                )

            shared_varstr = None
            for local_dict in per_sample_local.values():
                keys = set(local_dict.keys())
                shared_varstr = keys if shared_varstr is None else (shared_varstr & keys)

            if shared_varstr is None:
                shared_varstr = set()

            # build ND TimingMuts across carrier samples only
            event_sample_list = [self.sample_list[i] for i in matched_carrier_idx]
            supporting_muts = []

            for var_str in sorted(shared_varstr):
                ref_mut = per_sample_local[matched_carrier_idx[0]][var_str]
                m = len(matched_carrier_idx)

                alt_cnt = np.zeros(m)
                ref_cnt = np.zeros(m)
                ccf_dist = np.zeros((m, 101))
                local_cn_a1 = np.zeros(m)
                local_cn_a2 = np.zeros(m)

                for j, si in enumerate(matched_carrier_idx):
                    mut = per_sample_local[si][var_str]
                    alt_cnt[j] = mut.alt_cnt
                    ref_cnt[j] = mut.ref_cnt
                    ccf_dist[j, :] = mut.ccf_dist
                    local_cn_a1[j] = _scalar(mut.local_cn_a1)
                    local_cn_a2[j] = _scalar(mut.local_cn_a2)

                nd_mut = TimingMut(
                    event_sample_list,
                    ref_mut.gene,
                    ref_mut.chrN,
                    ref_mut.pos,
                    ref_mut.alt,
                    ref_mut.ref,
                    alt_cnt,
                    ref_cnt,
                    local_cn_a1,
                    local_cn_a2,
                    ccf_dist,
                    prot_change=ref_mut.prot_change,
                    cluster_assignment=ref_mut.cluster_assignment
                )
                supporting_muts.append(nd_mut)

                # avoid replacing an existing arm-level ND mutation if already present
                self.timeable_muts.setdefault(var_str, nd_mut)

            purity = [s.purity for s in event_sample_list]
            state = TimingCNState(
                event_sample_list,
                chrom,
                arm,
                consensus_state,
                purity,
                supporting_muts=supporting_muts
            )
            state.call_events(cn_type_prefix='Focal', wgd=(self.WGD is not None))

            # pick the focal event matching category and requested allele channel
            # Use approximate matching because _coerce_cn may keep fractional
            # values (e.g. 4.7) while call_events rounds to the nearest int
            # (e.g. 5.0), causing exact == to fail for high-level CN states.
            def _cn_close(a, b, tol=0.5):
                return abs(float(a) - float(b)) < tol

            target_cn = consensus_state[0] if is_a1 else consensus_state[1]
            eve = None
            best_dist = float('inf')
            for e in state.cn_events:
                if e.Type != category:
                    continue
                d = abs(float(e.allelic_cn) - float(target_cn))
                if d < best_dist and _cn_close(e.allelic_cn, target_cn):
                    eve = e
                    best_dist = d

            if eve is None and self.options.continuous_cn:
                # A partial step keeps a fractional state, so the event nearest in copy number is the one meant.
                candidates = [e for e in state.cn_events if e.Type == category]
                eve = min(candidates, key=lambda e: abs(float(e.allelic_cn) - float(target_cn))) if candidates else None

            if eve is None:
                continue

            eve.start = int(start_bp)
            eve.end = int(end_bp)
            eve.band_range_str = _band_range_str(chrom, int(start_bp), int(end_bp))
            eve.supporting_muts = supporting_muts

            # helpful metadata
            eve.n_local_supporting_muts = len(supporting_muts)
            eve.carrier_sample_idx = tuple(matched_carrier_idx)
            eve.carrier_sample_names = tuple(self.sample_list[i].sample_name for i in matched_carrier_idx)
            eve.local_cn_consensus = consensus_state
            # A high-level amplification in any carrier: its copy number says nothing of its order against a WGD.
            eve.high_level = any(rec['high_level'][i] for i in matched_carrier_idx)

            # optional hook for a borrowed ARM-LEVEL prior later in time_events()
            # This is metadata only; no borrowed evidence is mixed into the local likelihood here.
            if self.borrow_arm_prior:
                arm_type = 'Arm_gain' if 'gain' in category.lower() else 'Arm_loss'
                eve.borrowed_prior_key = (arm_type, str(chrom), arm, float(consensus_state[0]), float(consensus_state[1]))
            else:
                eve.borrowed_prior_key = None

            # clonal / trunkal classification still uses the full-sample focal CCF pattern
            ccf_idx = np.rint(
                np.clip(np.nan_to_num(np.asarray(rec['ccf'], dtype=float), nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0) * 100.0
            ).astype(int)

            clonal_grid = cluster_ccfs.get(1, np.ones((n_samples, 101)))
            clonal_conc = np.prod(clonal_grid[sample_idx, ccf_idx])

            is_clonal = True
            for cluster_id, grid in cluster_ccfs.items():
                if np.prod(grid[sample_idx, ccf_idx]) > clonal_conc:
                    is_clonal = False
                    break

            if self.options.continuous_cn:
                # The event must reach the clonal CCF in every sample, also a sample without it.
                called = self._clonality(float(np.min(rec['ccf'])))
                is_clonal = called is True
                eve.clonality_ambiguous = called is None
                eve.ccf_hat = float(np.min(rec['ccf']))

            eve.is_clonal = is_clonal

            # IMPORTANT:
            # do NOT assign uniform pi_dist here. Leave that decision to time_events().
            eve.pi_dist = None

            self.all_cn_events.setdefault(eve.event_name, []).append(eve)
            if is_clonal:
                self.truncal_cn_events.setdefault(eve.event_name, []).append(eve)

    def merge_whole_chromosome_events(self):
        """
        One event for a change on both arms of a chromosome in the same state, named for the chromosome (gain_17).

        A numerical gain of chromosome 17 otherwise reads as gain_17p and gain_17q, and so counts as the segmental 17q
        gain. In every sample, both arms must hold the same modal state and finite continuous values within the
        tolerance, and the event must be on the same allele of both. The patient-level states cannot decide this: each
        of their values comes from whichever sample is closest to the integer state. The event's mutations are those of
        both arms and its CCF the lower of the two. It is subclonal if either arm is, unknown if neither is but one is
        unknown, and clonal only if both are. Acrocentric chromosomes have no called p arm and keep their q-arm names.
        Terminal partial-arm states are left alone, and so is chrY, which is one event already.
        """
        tol = self.options.cn_event_tol

        def values(state):
            return state.continuous if state.continuous is not None else (state.cn_a1, state.cn_a2)

        def same(p, q):
            if p is None or q is None or (p.cn_a1, p.cn_a2) != (q.cn_a1, q.cn_a2):
                return False
            pairs = list(zip(values(p), values(q)))
            return all(np.isfinite(a) and np.isfinite(b) and abs(a - b) <= tol for a, b in pairs)

        for chrN in self.genome.CHROMS:
            p, q = self.concordant_cn_states.get(chrN + 'p'), self.concordant_cn_states.get(chrN + 'q')
            if not same(p, q) or not p.cn_events or not q.cn_events:
                continue
            if not all(same(sample.cn_states.get(chrN + 'p'), sample.cn_states.get(chrN + 'q'))
                       for sample in self.sample_list):
                continue
            arm_values = [values(p), values(q)]
            weights = np.array([self.genome.CENT_DICT[chrN], self.genome.CHROM_DICT[chrN] - self.genome.CENT_DICT[chrN]],
                               dtype=float)
            weights /= weights.sum()
            continuous = tuple(float(np.dot(weights, [arm_values[0][i], arm_values[1][i]])) for i in range(2))
            sds = [s.continuous_sd if s.continuous_sd is not None else (0., 0.) for s in (p, q)]
            continuous_sd = tuple(float(np.sqrt(np.sum((weights * np.array([sds[0][i], sds[1][i]])) ** 2)))
                                  for i in range(2))
            if p.continuous is None and q.continuous is None:
                # Integer states, as without continuous_cn: the merged events must be called the same way as the arms'.
                continuous = continuous_sd = None
            muts = list(p.supporting_muts) + [mut for mut in q.supporting_muts if mut not in p.supporting_muts]
            state = TimingCNState(p.sample_list, chrN, '', (p.cn_a1, p.cn_a2), p.purity, supporting_muts=muts,
                                  continuous=continuous, continuous_sd=continuous_sd, cn_type_prefix='Chrom')
            state.call_events(wgd=self.WGD is not None)
            state.merged_arms = (p, q)
            # Integer states give both gains of a 2/2 arm the same type and no allele, so each event pairs only once.
            q_left, merged_left = list(q.cn_events), list(state.cn_events)
            for p_event in list(p.cn_events):
                kind = p_event.Type.split('_', 1)[1]
                q_event = next((e for e in q_left if e.Type == p_event.Type and e.allele == p_event.allele), None)
                event = next((e for e in merged_left if e.Type == 'Chrom_' + kind and e.allele == p_event.allele), None)
                if q_event is None or event is None:
                    continue
                q_left.remove(q_event)
                merged_left.remove(event)
                pair = (p_event, q_event)
                lower = min(pair, key=lambda e: -1. if e.ccf_hat is None else e.ccf_hat)
                event.ccf_hat = lower.ccf_hat
                subclonal = any(not e.is_clonal and not e.clonality_ambiguous for e in pair)
                event.is_clonal = bool(p_event.is_clonal and q_event.is_clonal)
                event.clonality_ambiguous = not subclonal and not event.is_clonal
                event.merged_arms = pair
                for events in (self.all_cn_events, self.truncal_cn_events):
                    for arm_event in pair:
                        if arm_event in events.get(arm_event.event_name, []):
                            events[arm_event.event_name].remove(arm_event)
                    if events is self.all_cn_events or event.is_clonal:
                        events.setdefault(event.event_name, []).append(event)

    def _get_cluster_ccfs(self):
        n_samples = len(self.sample_list)
        cluster_ccfs = {}
        for mut in self.mutations.values():
            c = mut.cluster_assignment
            if c is not None:
                cluster_ccfs.setdefault(c, np.zeros((n_samples, 101)))
                cluster_ccfs[c] += np.log(mut.ccf_dist + 1e-10)
        for c in cluster_ccfs:
            cluster_ccfs[c] = np.exp(cluster_ccfs[c] - logsumexp(cluster_ccfs[c], axis=1, keepdims=True))
        return cluster_ccfs

    def time_events(self):
        """
        Timing order:
          1) If WGD exists: time WGD, then mutations relative to WGD, then CN events relative to WGD.
          2) Otherwise:
             a) time arm-level clonal gains from local supporting mutations
             b) build optional arm-level prior lookup
             c) time focal events using local evidence and optional arm prior
             d) assign remaining CN events to subclonal / not-estimable states
          3) assign fallback mutation distributions for untimed mutations
        """
        uniform_dist = np.ones(101, dtype=float) / 101.0
        subclonal_dist = np.zeros(101, dtype=float)
        subclonal_dist[100] = 1.0

        # ------------------------------------------------------------------
        # Small compatibility wrappers in case helper methods were not added
        # ------------------------------------------------------------------
        def _set_not_estimable(ev, reason='insufficient_local_support'):
            if self.options.clonal_unresolved and getattr(ev, 'is_clonal', False):
                ev.pi_dist = None
                ev.timing_status = 'clonal_unresolved'
                ev.timing_source = None
                ev.timing_reason = reason
                return
            if hasattr(ev, 'set_not_estimable'):
                ev.set_not_estimable(reason)
            else:
                ev.pi_dist = None
                if hasattr(ev, 'timing_status'):
                    ev.timing_status = 'not_estimable'
                    ev.timing_source = None
                    ev.timing_reason = reason

        def _set_prior_only(ev, pi_dist, source='arm_prior', reason='borrowed_prior_only'):
            pi_dist = np.asarray(pi_dist, dtype=float)
            pi_dist = pi_dist / pi_dist.sum()
            if hasattr(ev, 'set_prior_only'):
                ev.set_prior_only(pi_dist, source=source, reason=reason)
            else:
                ev.pi_dist = pi_dist
                if hasattr(ev, 'timing_status'):
                    ev.timing_status = 'prior_only'
                    ev.timing_source = source
                    ev.timing_reason = reason

        def _set_status(ev, status=None, source=None, reason=None):
            if hasattr(ev, 'timing_status'):
                ev.timing_status = status
                ev.timing_source = source
                ev.timing_reason = reason

        def _normalize(dist):
            dist = np.asarray(dist, dtype=float)
            s = dist.sum()
            if s <= 0:
                return None
            return dist / s

        def _fuse_prior(post, prior):
            post = _normalize(post)
            prior = _normalize(prior)
            if post is None or prior is None:
                return None
            fused = post * prior
            return _normalize(fused)

        def _is_gain(ev):
            return 'gain' in ev.Type.lower()

        def _is_focal(ev):
            return ev.Type.startswith('Focal_')

        def _is_arm(ev):
            return ev.Type.startswith(('Arm_', 'Chrom_'))

        def _has_support(ev):
            if self.options.clonal_unresolved:
                # Only clonal mutations carry multiplicity evidence. Otherwise the result is the prior alone.
                return sum(1 for mut in ev.supporting_muts if mut.is_clonal) >= self.min_supporting_muts
            return len(ev.supporting_muts) >= self.min_supporting_muts

        def _state_ok(ev):
            return (ev.cn_a1, ev.cn_a2) in self.cn_state_whitelist

        def _time_gain_from_local(ev, fuse_with_prior=None):
            """
            Time a clonal gain from local interval mutations.
            Returns True if successful.
            """
            if not (_is_gain(ev) and _has_support(ev) and _state_ok(ev)):
                return False

            ev.get_pi_dist_for_gain()
            if ev.pi_dist is None:
                return False

            if fuse_with_prior is not None:
                fused = _fuse_prior(ev.pi_dist, fuse_with_prior)
                if fused is not None:
                    ev.pi_dist = fused
                    _set_status(ev, status='estimated', source='local+arm_prior', reason=None)
                else:
                    _set_status(ev, status='estimated', source='local', reason=None)
            else:
                _set_status(ev, status='estimated', source='local', reason=None)

            # time clonal supporting muts relative to this event
            for mut in ev.supporting_muts:
                if mut.is_clonal:
                    mut.get_pi_dist(ev)

            return True

        # ------------------------------------------------------------------
        # WGD path
        # ------------------------------------------------------------------
        if self.WGD is not None:
            self.WGD.get_pi_dist()
            if self.options.clonal_unresolved and self.WGD.n_timed_gains == 0:
                # A WGD with no timed gain has a uniform distribution, which is no evidence to time others against.
                self.WGD.pi_dist = None
                self.WGD.timing_status = 'clonal_unresolved'
            else:
                self.WGD.timing_status = 'estimated'

            for mut in self.mutations.values():
                if mut.is_clonal:
                    mut.get_pi_dist(self.WGD)

            for evs in self.all_cn_events.values():
                for cn_event in evs:
                    if cn_event.clonality_ambiguous:
                        cn_event.pi_dist = None
                        _set_status(cn_event, status='ambiguous_clonality', source='clonality_rule', reason='ccf_near_cutoff')
                    elif not cn_event.is_clonal:
                        cn_event.pi_dist = subclonal_dist.copy()
                        _set_status(cn_event, status='subclonal', source='clonality_rule', reason='event_not_truncal')
                    elif self.WGD.pi_dist is None:
                        _set_not_estimable(cn_event, reason='wgd_not_timed')
                    elif cn_event.Type.endswith('homdel') and self.options.segment_homdels:
                        # No copy left: by parsimony both alleles were lost before the genome doubled.
                        cn_event.get_pi_dist_for_loss(self.WGD)
                        _set_status(cn_event, status='estimated', source='wgd_rule', reason=None)
                    elif cn_event.Type.endswith('homdel'):
                        # Homdel: no local mutations possible; mark as present
                        cn_event.pi_dist = uniform_dist.copy()
                        _set_status(cn_event, status='prior_only', source='homdel',
                                    reason='homdel_no_local_evidence')
                    elif cn_event.Type.endswith('gain') and getattr(cn_event, 'high_level', False):
                        # The parity rule counts copies made by doubling the genome. An amplicon, often extrachromosomal,
                        # gains copies by its own route, so neither parity nor level places it against the WGD.
                        _set_not_estimable(cn_event, reason='amplicon_vs_wgd')
                    elif cn_event.Type.endswith('gain'):
                        cn_event.get_pi_dist_for_higher_gain(self.WGD)
                        if cn_event.pi_dist is not None:
                            _set_status(cn_event, status='estimated', source='wgd_rule', reason=None)
                        else:
                            _set_not_estimable(cn_event, reason='unsupported_gain_state')
                    elif self._loy_before_wgd_unsupported(cn_event):
                        _set_not_estimable(cn_event, reason='loy_wgd_odds_below_min')
                    elif cn_event.Type.endswith('loss'):
                        cn_event.get_pi_dist_for_loss(self.WGD)
                        if cn_event.pi_dist is not None:
                            _set_status(cn_event, status='estimated', source='wgd_rule', reason=None)
                        else:
                            _set_not_estimable(cn_event, reason='unsupported_loss_state')
                    else:
                        _set_not_estimable(cn_event, reason='unsupported_event_type')
                    if self.options.clonal_unresolved and cn_event.timing_source == 'wgd_rule':
                        # Events the rule puts on one side of the WGD share one distribution, which orders neither.
                        cn_event.wgd_side = cn_event.rule_side

        # ------------------------------------------------------------------
        # No-WGD path
        # ------------------------------------------------------------------
        else:
            arm_events = []
            focal_events = []
            other_events = []

            for evs in self.all_cn_events.values():
                for ev in evs:
                    if _is_arm(ev):
                        arm_events.append(ev)
                    elif _is_focal(ev):
                        focal_events.append(ev)
                    else:
                        other_events.append(ev)

            # --------------------------------------------------------------
            # Pass 1: time arm-level events first
            # --------------------------------------------------------------
            arm_prior_lookup = {}

            for ev in arm_events:
                if ev.clonality_ambiguous:
                    ev.pi_dist = None
                    _set_status(ev, status='ambiguous_clonality', source='clonality_rule', reason='ccf_near_cutoff')
                    continue

                if not ev.is_clonal:
                    ev.pi_dist = subclonal_dist.copy()
                    _set_status(ev, status='subclonal', source='clonality_rule', reason='event_not_truncal')
                    continue

                if _time_gain_from_local(ev):
                    arm_prior_lookup[(ev.Type, ev.chrN, ev.arm, ev.cn_a1, ev.cn_a2)] = ev.pi_dist.copy()
                    for arm_event in getattr(ev, 'merged_arms', ()):
                        # Focal events look up the prior of their arm.
                        arm_prior_lookup[(arm_event.Type, ev.chrN, arm_event.arm, ev.cn_a1, ev.cn_a2)] = ev.pi_dist.copy()
                else:
                    # In the no-WGD setting, non-gains / unsupported gains are not event-specifically timeable
                    _set_not_estimable(
                        ev,
                        reason='insufficient_local_support' if not _has_support(ev) else
                        'unsupported_cn_state' if not _state_ok(ev) else
                        'unsupported_event_type'
                    )

            # --------------------------------------------------------------
            # Pass 2: time focal events using local evidence + optional arm prior
            # --------------------------------------------------------------
            for ev in focal_events:
                if ev.clonality_ambiguous:
                    ev.pi_dist = None
                    _set_status(ev, status='ambiguous_clonality', source='clonality_rule', reason='ccf_near_cutoff')
                    continue

                if not ev.is_clonal:
                    ev.pi_dist = subclonal_dist.copy()
                    _set_status(ev, status='subclonal', source='clonality_rule', reason='event_not_truncal')
                    continue

                borrowed_prior = None
                if self.borrow_arm_prior:
                    prior_key = getattr(ev, 'borrowed_prior_key', None)
                    if prior_key is not None:
                        borrowed_prior = arm_prior_lookup.get(prior_key)

                # Preferred path: local likelihood from focal interval mutations
                if _time_gain_from_local(ev, fuse_with_prior=borrowed_prior):
                    continue

                # Homdel events cannot be timed from local mutations (no reads
                # exist in a homozygously deleted region), but they should be
                # present in the output as identifiable events.
                if 'homdel' in ev.Type.lower() and self.options.clonal_unresolved:
                    _set_not_estimable(ev, reason='homdel_no_local_evidence')
                    continue
                if 'homdel' in ev.Type.lower():
                    ev.pi_dist = uniform_dist.copy()
                    _set_status(ev, status='prior_only', source='homdel',
                                reason='homdel_no_local_evidence')
                    continue

                # No local evidence, but explicit borrowing enabled
                if borrowed_prior is not None:
                    _set_prior_only(ev, borrowed_prior, source='arm_prior', reason='borrowed_prior_only')
                    continue

                if self.options.clonal_unresolved:
                    _set_not_estimable(
                        ev,
                        reason='insufficient_local_support' if not _has_support(ev) else
                        'unsupported_cn_state' if not _state_ok(ev) else
                        'unsupported_event_type'
                    )
                    continue

                # Focal events that cannot be locally timed still get a
                # uniform pi_dist so they remain visible as "present but
                # uninformative" in comparisons and downstream output.
                ev.pi_dist = uniform_dist.copy()
                _set_status(
                    ev,
                    status='not_estimable',
                    source='uniform_fallback',
                    reason='insufficient_local_support' if not _has_support(ev) else
                    'unsupported_cn_state' if not _state_ok(ev) else
                    'unsupported_event_type'
                )

            # --------------------------------------------------------------
            # Pass 3: any other CN events (future-proofing)
            # --------------------------------------------------------------
            for ev in other_events:
                if not ev.is_clonal:
                    ev.pi_dist = subclonal_dist.copy()
                    _set_status(ev, status='subclonal', source='clonality_rule', reason='event_not_truncal')
                else:
                    _set_not_estimable(ev, reason='unsupported_event_type')

        # ------------------------------------------------------------------
        # Final fallback for SNVs/indels
        # ------------------------------------------------------------------
        for mut in self.mutations.values():
            if mut.pi_dist is None:
                if mut.is_clonal:
                    mut.pi_dist = uniform_dist.copy()
                    if self.options.clonal_unresolved:
                        mut.timing_status = 'clonal_unresolved'
                else:
                    mut.pi_dist = subclonal_dist.copy()
                    if self.options.clonal_unresolved:
                        mut.timing_status = 'subclonal'
            elif self.options.clonal_unresolved:
                mut.timing_status = 'estimated'


class TimingSample(object):
    """
    class for holding a sample
    """
    def __init__(self, sample, engine, cn_state_whitelist=Enums.cn_state_whitelist):
        self.sample = sample
        self.sample_name = self.sample.sample_name
        self.engine = engine
        self.options = engine.options
        self.cn_state_whitelist = cn_state_whitelist
        self.purity = self.sample.purity
        self.CnProfile = self.sample.CnProfile
        self._chromosomes = self.engine.genome.CHROMS
        self.arm_regions = self.engine.arm_regions
        self.mutation_intervaltree = {chrom: IntervalTree() for chrom in self._chromosomes}
        self.mutations = []
        self.mut_lookup_table = {}
        self.concordant_mutation_intervaltree = {chrom: IntervalTree() for chrom in self._chromosomes}
        self.WGD = None
        self.concordant_WGD = None
        self.timeable_cn_events = None
        self.fill_mutation_intervaltree()
        self.cn_states = None
        self.missing_arms = None
        self.concordant_cn_states = None
        self.call_arm_level_cn_states()
        self.cn_events = {gl + chrN + arm: [] for gl, (chrN, arm) in itertools.product(('gain_', 'loss_'), self.arm_regions)}
        self.concordant_cn_events = {gl + chrN + arm: [] for gl, (chrN, arm) in itertools.product(('gain_', 'loss_'), self.arm_regions)}
        self.call_wgd()
        self.get_arm_level_cn_events()

    def __repr__(self):
        return '<TimingSample object: {}>'.format(self.sample_name)

    def _continuous_arm_state(self, segments, size_threshold, hemizygous=False):
        """Arm state from continuous copy number, or None if no single state holds enough of the arm.

        Grouping on continuous values keeps a state whose segments round to different integers, as happens near CCF 0.5.
        The group grows around the modal state holding the most length, because a median can fall between two states.
        Returns the modal state with the most length in the group, the continuous values and their uncertainty.
        """
        tol = self.options.cn_event_tol
        lengths = np.array([length for _, length in segments], dtype=float)
        rows = [self.sample.allele_copy_number(seg.data[1]) for seg, _ in segments]
        a1, a2 = (np.array([row[i] for row in rows], dtype=float) for i in range(2))
        # One copy changes the copy ratio half as much on a diploid chromosome, split over its two homologs.
        noise = self.sample.copy_ratio_line()[0] / (2. if hemizygous else np.sqrt(2.))
        sd1 = sd2 = np.array([_number(seg.data[1].get('seg_sigma')) * noise for seg, _ in segments], dtype=float)
        modal = [(seg.data[1]['cn_a1'], seg.data[1]['cn_a2']) for seg, _ in segments]
        # Where a value is missing, fall back to the modal state of that segment.
        a1 = np.where(np.isfinite(a1), a1, [m[0] for m in modal])
        a2 = np.where(np.isfinite(a2), a2, [m[1] for m in modal])

        def group(mask):
            by_modal = Counter()
            for index in np.flatnonzero(mask):
                by_modal[modal[index]] += lengths[index]
            largest = by_modal.most_common(1)[0][0]
            seed = np.array([m == largest for m in modal]) & mask
            center = (_weighted_median(a1[seed], lengths[seed]), _weighted_median(a2[seed], lengths[seed]))
            return center, mask & (np.abs(a1 - center[0]) <= tol) & (np.abs(a2 - center[1]) <= tol)

        total = lengths.sum()
        everything = np.ones(len(lengths), dtype=bool)
        center, inside = group(everything)
        if lengths[inside].sum() <= size_threshold * total:
            return None
        if (~inside).any():
            _, second = group(~inside)
            if lengths[second].sum() > size_threshold * total:
                return None
        by_state = Counter()
        for index in np.flatnonzero(inside):
            by_state[modal[index]] += lengths[index]
        weights = lengths[inside] / lengths[inside].sum()
        sds = tuple(float(np.sqrt(np.nansum((weights * sd[inside]) ** 2))) for sd in (sd1, sd2))
        return by_state.most_common(1)[0][0], center, sds

    def __getitem__(self, key):
        return self.mut_lookup_table[key]

    def fill_mutation_intervaltree(self, concordant_muts=None):
        """
        attach local copy number to mutations and attach mutations to intervaltree for easy lookup by position
        """
        if concordant_muts is None:  # concordant_muts are truncal mutations shared by all samples
            self.mutation_intervaltree = {chrom: IntervalTree() for chrom in self._chromosomes}
            self.mutations = []
            self.mut_lookup_table = {}
        else:
            self.concordant_mutation_intervaltree = {chrom: IntervalTree() for chrom in self._chromosomes}
        for mut in itertools.chain(self.sample.mutations, self.sample.low_coverage_mutations.values()):
            if mut.chrN in self.mutation_intervaltree and (concordant_muts is None or mut.var_str in concordant_muts) and mut.type != "CNV":
                segs = self.CnProfile[mut.chrN][mut.pos]
                if segs:
                    cn_seg = next(iter(segs))
                    local_cn_a1 = float(cn_seg.data[1]['cn_a1'])
                    local_cn_a2 = float(cn_seg.data[1]['cn_a2'])
                else:
                    local_cn_a1 = np.nan
                    local_cn_a2 = np.nan
                timing_mut = TimingMut([self], mut.gene, mut.chrN, mut.pos, mut.alt, mut.ref, mut.alt_cnt, mut.ref_cnt,
                                       np.array([local_cn_a1]), np.array([local_cn_a2]), mut.ccf_1d,
                                       prot_change=mut.prot_change, cluster_assignment=mut.cluster_assignment)
                if not np.isnan(local_cn_a1) and not np.isnan(local_cn_a2):
                    timing_mut.get_mult_dist()
                if concordant_muts is None:
                    self.mutations.append(timing_mut)
                    self.mutation_intervaltree[mut.chrN][mut.pos:mut.pos + 1] = timing_mut
                    self.mut_lookup_table[mut.var_str] = timing_mut
                else:
                    self.concordant_mutation_intervaltree[mut.chrN][mut.pos:mut.pos + 1] = timing_mut

    def call_arm_level_cn_states(self, size_threshold=.4, concordant_states=None, tol=.2):
        if concordant_states is None:
            self.cn_states = {}
        else:
            self.concordant_cn_states = {}
        self.missing_arms = []
        full_segtree = self.CnProfile
        for chrN, arm in self.arm_regions:
            centromere = self.engine.genome.CENT_DICT[chrN]
            if arm == 'p':
                arm_segtree = full_segtree[chrN][:centromere]
                start = 0
                end = centromere
            else:
                arm_segtree = full_segtree[chrN][centromere:]
                start = centromere
                end = self.engine.genome.CHROM_DICT[chrN]
            if self.options.male and chrN == 'Y':
                self.missing_arms.append((chrN, arm))  # chrY is called as one event over callable territory
                continue
            true_arm_bp = end - start
            arm_bp = 0
            state_bps = {}
            arm_segments = []
            for seg in arm_segtree:
                if seg.begin < start:
                    seg_len = seg.end - start
                elif seg.end > end:
                    seg_len = end - seg.begin
                else:
                    seg_len = seg.length()
                arm_bp += seg_len
                cn_a1 = seg.data[1]['cn_a1']
                if abs(cn_a1 - round(cn_a1)) <= .2:
                    cn_a1 = round(cn_a1)
                cn_a2 = seg.data[1]['cn_a2']
                if abs(cn_a2 - round(cn_a2)) <= .2:
                    cn_a2 = round(cn_a2)
                state_tuple = (cn_a1, cn_a2)
                state_bps.setdefault(state_tuple, 0)
                state_bps[state_tuple] += seg_len
                arm_segments.append((seg, seg_len))
            states_over_threshold = [cn for cn, bp in state_bps.items() if bp > arm_bp * size_threshold]
            continuous, continuous_sd = None, None
            if arm_bp >= true_arm_bp * .5 and self.options.continuous_cn:
                call = self._continuous_arm_state(arm_segments, size_threshold, self.options.is_hemizygous(chrN))
                states_over_threshold = [] if call is None else [call[0]]
                if call is not None:
                    continuous, continuous_sd = call[1], call[2]
            if arm_bp >= true_arm_bp * .5 and self.options.terminal_arm_min_fraction is not None:
                self._add_terminal_state(chrN, arm, start, end, arm_segments, None if call is None else call[1],
                                         concordant_states)
            if arm_bp < true_arm_bp * .5:  # Only call if at least 50% of arm is in segs
                self.missing_arms.append((chrN, arm))
            elif len(states_over_threshold) == 1:
                cn_a1, cn_a2 = next(iter(states_over_threshold))
                if concordant_states is None:
                    supporting_muts = [interval.data for interval in self.mutation_intervaltree[chrN][start:end]
                                       if interval.data.local_cn_a1 == cn_a1 and interval.data.local_cn_a2 == cn_a2]
                    cn_state = TimingCNState([self], chrN, arm, (cn_a1, cn_a2), [self.purity],
                                             supporting_muts=supporting_muts, continuous=continuous,
                                             continuous_sd=continuous_sd)
                    self.cn_states[chrN + arm] = cn_state
                elif chrN + arm in concordant_states:
                    supporting_muts = [interval.data for interval in self.concordant_mutation_intervaltree[chrN][start:end]
                                      if interval.data.local_cn_a1 == cn_a1 and interval.data.local_cn_a2 == cn_a2]
                    cn_state = TimingCNState([self], chrN, arm, (cn_a1, cn_a2), [self.purity],
                                             supporting_muts=supporting_muts, continuous=continuous,
                                             continuous_sd=continuous_sd)
                    self.concordant_cn_states[chrN + arm] = cn_state
        if self.options.male:
            y_state = self._call_y_state()
            if y_state is not None:
                if concordant_states is None:
                    self.cn_states['Y'] = y_state
                elif 'Y' in concordant_states:
                    self.concordant_cn_states['Y'] = y_state

    def _add_terminal_state(self, chrN, arm, start, end, arm_segments, arm_center, concordant_states):
        """Keep a terminal run that the whole-arm call missed as its own state of this arm, under region + TERMINAL."""
        call = self._terminal_arm_state(arm_segments, arm, start, end, self.options.is_hemizygous(chrN))
        if call is None or call['fraction'] < self.options.terminal_arm_min_fraction or \
                call['span'][1] - call['span'][0] < self.options.terminal_arm_min_length:
            return
        tol = self.options.cn_event_tol
        if arm_center is not None and all(abs(a - b) <= tol for a, b in zip(arm_center, call['center'])):
            return  # the whole-arm call is this state already
        region = chrN + arm + TERMINAL
        lo, hi = call['span']
        cn_a1, cn_a2 = call['state']
        if concordant_states is None:
            tree = self.mutation_intervaltree[chrN]
        elif region in concordant_states:
            tree = self.concordant_mutation_intervaltree[chrN]
        else:
            return
        supporting_muts = [interval.data for interval in tree[lo:hi]
                           if interval.data.local_cn_a1 == cn_a1 and interval.data.local_cn_a2 == cn_a2]
        state = TimingCNState([self], chrN, arm, (cn_a1, cn_a2), [self.purity], supporting_muts=supporting_muts,
                              continuous=call['center'], continuous_sd=call['sd'])
        state.terminal = (lo, hi)
        if concordant_states is None:
            self.cn_states[region] = state
        else:
            self.concordant_cn_states[region] = state

    def _terminal_arm_state(self, segments, arm, start, end, hemizygous=False, min_purity=.9, max_gap=.1):
        """
        The longest run of segments from the telomeric end of an arm that holds one state, if the rest of the arm differs.

        Neuroblastoma's 1p and 11q losses are mostly terminal (1p36, 11q14-qter), and short of 40% of the arm the
        whole-arm call reads the normal remainder. A run holds one state if segments within the tolerance of its
        length-weighted median cover min_purity of its length, so short noisy segments do not break it, and it ends on
        such a segment. The rest of the arm must differ by more than the tolerance on one allele. The run must start
        within max_gap of the arm's length from the arm end. Returns None, or a dict with the modal state, the
        continuous center and its uncertainty, the span in bp, and the share of the arm the span covers.
        """
        if len(segments) < 2:
            return None
        tol = self.options.cn_event_tol
        order = sorted(segments, key=lambda item: item[0].begin, reverse=(arm == 'q'))
        arm_length = float(end - start)
        first = order[0][0]
        if (first.begin - start if arm == 'p' else end - first.end) > max_gap * arm_length:
            return None
        lengths = np.array([length for _, length in order], dtype=float)
        rows = [self.sample.allele_copy_number(seg.data[1]) for seg, _ in order]
        modal = [(seg.data[1]['cn_a1'], seg.data[1]['cn_a2']) for seg, _ in order]
        a1 = np.array([row[0] for row in rows], dtype=float)
        a2 = np.array([row[1] for row in rows], dtype=float)
        a1 = np.where(np.isfinite(a1), a1, [m[0] for m in modal])
        a2 = np.where(np.isfinite(a2), a2, [m[1] for m in modal])
        noise = self.sample.copy_ratio_line()[0] / (2. if hemizygous else np.sqrt(2.))
        sd = np.array([_number(seg.data[1].get('seg_sigma')) * noise for seg, _ in order], dtype=float)
        best = None
        for k in range(1, len(order)):
            center = (_weighted_median(a1[:k], lengths[:k]), _weighted_median(a2[:k], lengths[:k]))
            inside = (np.abs(a1[:k] - center[0]) <= tol) & (np.abs(a2[:k] - center[1]) <= tol)
            if not inside[-1] or lengths[:k][inside].sum() < min_purity * lengths[:k].sum():
                continue
            rest = (_weighted_median(a1[k:], lengths[k:]), _weighted_median(a2[k:], lengths[k:]))
            if all(abs(r - c) <= tol for r, c in zip(rest, center)):
                continue
            best = (k, center, inside)
        if best is None:
            return None
        k, center, inside = best
        inner = order[k - 1][0]
        span = (start, min(inner.end, end)) if arm == 'p' else (max(inner.begin, start), end)
        by_state = Counter()
        for index in np.flatnonzero(inside):
            by_state[modal[index]] += lengths[index]
        weights = lengths[:k][inside] / lengths[:k][inside].sum()
        sds = tuple(float(np.sqrt(np.nansum((weights * sd[:k][inside]) ** 2))) for _ in range(2))
        return {'state': by_state.most_common(1)[0][0], 'center': center, 'sd': sds, 'span': span,
                'fraction': (span[1] - span[0]) / arm_length}

    def _call_y_state(self, min_callable_fraction=.5):
        """Whole-chromosome chrY state over callable territory. None if too little of that territory is segmented."""
        blocks = self.options.y_callable_blocks
        tree = self.CnProfile['Y']
        total, values, lengths, sds = sum(end - start for start, end in blocks), [], [], []
        slope = self.sample.copy_ratio_line()[0]
        for start, end in blocks:
            for seg in tree.overlap(start, end):
                lengths.append(min(seg.end, end) - max(seg.begin, start))
                # Total copy number, because a male Y has one allele and the allelic columns are rescued, not modeled.
                value = _number(seg.data[1].get('total_cn_cont'))
                values.append(value if np.isfinite(value) else float(seg.data[1]['cn_a1']) + float(seg.data[1]['cn_a2']))
                # One Y copy moves the copy ratio twice as much as one copy of a diploid chromosome.
                sds.append(_number(seg.data[1].get('seg_sigma')) * slope / 2.)
        if sum(lengths) < min_callable_fraction * total:
            return None
        y_cn = _weighted_median(values, lengths)
        weights = np.array(lengths, dtype=float) / sum(lengths)
        sd = float(np.sqrt(np.nansum((weights * np.array(sds)) ** 2)))
        state = TimingCNState([self], 'Y', '', (0., float(round(y_cn))), [self.purity], supporting_muts=[],
                              continuous=(0., y_cn), continuous_sd=(0., sd), cn_type_prefix='Chrom')
        # Tumors without Y keep a copy-ratio floor of about 0.07, from cross-mapped reads. Low purity stretches it on the
        # copy-number scale, so the limit for no Y left is set in copy-ratio units there.
        if np.isfinite(slope):
            state.absent_tol = max(self.options.cn_event_tol, Y_COPY_RATIO_FLOOR * slope / 2.)
            state.call_events()
        return state

    def call_wgd(self, use_concordant_states=False):
        """
        call WGD for a sample by looking at supporting (0/2 and 2/2) regions
        """
        # TODO: filter outlier pi distributions
        regions_supporting_WGD = []
        regions_both_arms_gained = []
        supporting_arm_states = []
        if use_concordant_states:
            if self.WGD is None:
                return
            for cn_state in self.concordant_cn_states.values():
                if cn_state.cn_type_prefix == 'Chrom' or cn_state.terminal is not None:
                    continue
                if cn_state.cn_a1 == 0 or cn_state.cn_a1 >= 2 and cn_state.cn_a2 >= 2:  # region supports WGD if 0/2 or 2/2
                    supporting_arm_states.append(TimingCNState([self], cn_state.chrN, cn_state.arm,
                        (cn_state.cn_a1, cn_state.cn_a2), cn_state.purity, supporting_muts=cn_state.supporting_muts))
            self.concordant_WGD = TimingWGD(supporting_arm_states=supporting_arm_states)
            for cn_state in self.concordant_cn_states.values():
                cn_state.call_events(wgd=True)
        else:
            for cn_state in self.cn_states.values():
                if cn_state.cn_type_prefix == 'Chrom' or cn_state.terminal is not None:
                    continue
                if cn_state.cn_a2 >= 2:
                    supporting_arm_states.append(cn_state)
                if (cn_state.cn_a1 == 0 or cn_state.cn_a1 >= 2) and cn_state.cn_a2 >= 2:  # region supports WGD if 0/2 or 2/2
                    regions_supporting_WGD.append(cn_state)
                if cn_state.cn_a1 >= 2 and cn_state.cn_a2 >= 2:
                    regions_both_arms_gained.append(cn_state)
            doubled = len(regions_both_arms_gained) >= 5 and len(regions_supporting_WGD) * 2 >= \
                len(self.arm_regions) - len(self.missing_arms)
            override = self.options.wgd_override.get(self.sample_name)
            if override is not None:
                if override != doubled:
                    print('WGD of sample {} set to {} as given, against {} by the arm rule ({} of {} called arms '
                          'support it)'.format(self.sample_name, override, doubled, len(regions_supporting_WGD),
                                               len(self.arm_regions) - len(self.missing_arms)))
                doubled = override
            if doubled:
                supporting_arm_states = [TimingCNState([self], s.chrN, s.arm, (s.cn_a1, s.cn_a2), s.purity, supporting_muts=s.supporting_muts) for
                                         s in supporting_arm_states]
                self.WGD = TimingWGD(supporting_arm_states=supporting_arm_states)
                for cn_state in self.cn_states.values():
                    cn_state.call_events(wgd=True)

    def get_arm_level_cn_events(self, use_concordant_states=False):
        """
        extract cn events from cn state instances
        """
        if use_concordant_states:
            self.concordant_cn_events = {gl + chrN + arm: [] for gl, (chrN, arm) in itertools.product(('gain_', 'loss_'), self.arm_regions)}
        else:
            self.cn_events = {gl + chrN + arm: [] for gl, (chrN, arm) in itertools.product(('gain_', 'loss_'), self.arm_regions)}
        for state in self.cn_states.values():
            for eve in state.cn_events:
                if use_concordant_states:
                    self.concordant_cn_events.setdefault(eve.event_name, []).append(eve)
                else:
                    self.cn_events.setdefault(eve.event_name, []).append(eve)


class TimingWGD(object):
    """
    Class for holding WGD event
    """
    Type = 'WGD'
    event_name = 'WGD'

    def __init__(self, supporting_arm_states=(), pi_dist=None, min_timeable_muts=3):
        self.supporting_arm_states = supporting_arm_states
        self.pi_dist = pi_dist
        self.min_timeable_muts = min_timeable_muts
        self.supporting_muts = []
        self.timing_status = None
        self.n_timed_gains = 0
        self.get_supporting_muts()

    def get_pi_dist(self):
        """
        Get pi dist by multiplying pi dists for supporting arm gains
        """
        pi_dist = np.zeros(101)
        self.n_timed_gains = 0
        for state in self.supporting_arm_states:
            for eve in state.cn_events:
                if eve.Type == 'Arm_gain':
                    if state.options.clonal_unresolved:
                        n_timeable = sum(1 for mut in eve.supporting_muts if mut.is_clonal)
                    else:
                        n_timeable = len(eve.supporting_muts)
                    if n_timeable < self.min_timeable_muts:
                        continue
                    if eve.pi_dist is None:
                        eve.get_pi_dist_for_gain()
                    if eve.pi_dist is None:
                        continue
                    pi_dist += np.log(eve.pi_dist)
                    self.n_timed_gains += 1
                    # A 2/2 arm has two gains built from the same mutations, which are one piece of evidence.
                    if state.options.extended_output:
                        break
        self.pi_dist = np.exp(pi_dist - logsumexp(pi_dist))

    def get_supporting_muts(self):
        self.supporting_muts = []
        for state in self.supporting_arm_states:
            self.supporting_muts.extend(state.supporting_muts)


class TimingCNState(object):
    """
    Arm events are timed based on copy number state of both arms
    """

    def __init__(self, sample_list, chrN, arm, cn_state, purity, supporting_muts=None, cn_type_prefix='Arm',
                 continuous=None, continuous_sd=None):
        self.sample_list = sample_list
        self.options = getattr(sample_list[0], 'options', None) or TimingOptions()
        self.chrN = chrN
        self.arm = arm
        self.cn_a1, self.cn_a2 = cn_state
        self.purity = purity
        self.supporting_muts = supporting_muts
        self.cn_type_prefix = cn_type_prefix
        self.continuous = continuous  # (a1, a2) continuous copy number, or None to use the integer state
        self.continuous_sd = continuous_sd
        self.absent_tol = None
        self.terminal = None  # (start, end) of a terminal partial-arm run, see TimingSample._terminal_arm_state
        self.haploid = self.options.is_hemizygous(chrN)
        self.cn_events = None
        self.call_events(cn_type_prefix=cn_type_prefix)

    def __repr__(self):
        return '<TimingCNState object: {}{}-{}:{}>'.format(self.chrN, self.arm, self.cn_a1, self.cn_a2)

    @property
    def p2_domain(self):
        if self.cn_a2 == 2.:
            if self.cn_a1 == 0. or self.cn_a1 == 2.:
                return np.linspace(0, 1, 101)
            if self.cn_a1 == 1.:  # for 1/2 even if event happened last p2 should not exceed .5
                return np.linspace(0, .5, 101)

    def _call_events_continuous(self, cn_type_prefix, wgd):
        """
        Events relative to the normal state, with presence and CCF taken from continuous copy number where available.
        A hemizygous chromosome has one allele, so a1 = 0 is its normal state and never an event.

        Each event is the step from the normal state to the next copy number, and ccf_hat is the share of tumor cells
        carrying that step. Its allelic_cn is the copy number the event reached in all cells, which decides its order
        against a WGD: only a value consistent with no copy left can have been lost before doubling. The event keeps the
        modal state, because mutations are matched to it and timed from it.
        """
        baseline = 2. if wgd else 1.
        tol = self.options.cn_event_tol
        continuous = self.continuous if self.continuous is not None else (self.cn_a1, self.cn_a2)
        values = [continuous[i] if np.isfinite(continuous[i]) else state for i, state in enumerate((self.cn_a1, self.cn_a2))]
        sds = self.continuous_sd if self.continuous_sd is not None else (0., 0.)

        def step(value):
            """(allelic_cn, ccf_hat) of the event, or None if value is the normal state."""
            if abs(value - baseline) < tol:
                return None
            if value < baseline:
                if value < (self.absent_tol or tol):
                    return 0., (1. if self.absent_tol else min(baseline - value, 1.))
                # A partial drop below the next integer is that step in all cells, with a further loss in some of them.
                return baseline - 1., min(baseline - value, 1.)
            if value <= baseline + 1:
                return baseline + 1., value - baseline
            return float(np.floor(value + tol)), 1.

        kwargs = dict(chrN=self.chrN, arm=self.arm, supporting_muts=self.supporting_muts)
        copy_number = str(self.cn_a1) + '/' + str(self.cn_a2)
        no_copy = self.absent_tol or tol
        lost_a1 = self.haploid or values[0] < no_copy
        lost_a2 = values[1] < no_copy
        if lost_a1 and lost_a2 and cn_type_prefix == 'Focal':
            event = TimingCNEvent(self.sample_list, self, Type=cn_type_prefix + '_homdel', copy_number='0/0',
                                  allelic_cn=0, ccf_hat=1., **kwargs)
            self.cn_events.append(event)
            return
        alleles = [('a2', values[1], sds[1])] if self.haploid else [('a1', values[0], sds[0]), ('a2', values[1], sds[1])]
        for allele, value, sd in alleles:
            called = step(value)
            if called is None:
                continue
            cn_type = cn_type_prefix + ('_gain' if value > baseline else '_loss')
            event = TimingCNEvent(self.sample_list, self, Type=cn_type, allelic_cn=called[0], copy_number=copy_number,
                                  ccf_hat=called[1], **kwargs)
            event.allele = allele
            event.ccf_sd = sd
            event.continuous_cn = value
            # A gain far above the next copy cannot be subclonal for want of cells: the clonality test reads the excess
            # before the CCF is capped at 1, or the noise of a high gain would make it unknown.
            event.clonality_ccf = value - baseline if value > baseline else called[1]
            self.cn_events.append(event)

    def call_events(self, cn_type_prefix=None, wgd=False):
        """
        create cn event instances for events (gain for 1/2, gain and loss for 0/2, two gains for 2/2)
        Also handles homozygous deletion (0/0) as a single homdel event.
        """
        cn_type_prefix = cn_type_prefix or self.cn_type_prefix
        self.cn_events = []
        if self.haploid or self.continuous is not None:
            self._call_events_continuous(cn_type_prefix, wgd)
            return
        baseline = 2. if wgd else 1.
        if baseline - 1 <= self.cn_a1 < baseline:
            cn_a1 = baseline - 1
        elif baseline < self.cn_a1 <= baseline + 1:
            cn_a1 = baseline + 1
        else:
            cn_a1 = round(self.cn_a1)
        if baseline - 1 <= self.cn_a2 < baseline:
            cn_a2 = baseline - 1
        elif baseline < self.cn_a2 <= baseline + 1:
            cn_a2 = baseline + 1
        else:
            cn_a2 = round(self.cn_a2)

        # Homozygous deletion: both alleles at 0 — only meaningful for focal events.
        # Arm-level (0,0) is almost certainly due to inaccurate CN inference and
        # should be treated as two independent losses (original behavior).
        if cn_a1 == 0 and cn_a2 == 0 and (not wgd or self.options.segment_homdels) and cn_type_prefix == 'Focal':
            cn_type = cn_type_prefix + '_homdel'
            ccf_hat = min(baseline, 1.)
            self.cn_events.append(TimingCNEvent(self.sample_list, self, Type=cn_type, chrN=self.chrN, arm=self.arm,
                                                copy_number='0/0', allelic_cn=0,
                                                supporting_muts=self.supporting_muts, ccf_hat=ccf_hat))
            return

        if cn_a1 != baseline:
            cn_type = cn_type_prefix + ('_gain' if self.cn_a1 > baseline else '_loss')
            ccf_hat = min(self.cn_a1 - baseline, 1.) if self.cn_a1 > baseline else min(baseline - self.cn_a1, 1.)
            self.cn_events.append(TimingCNEvent(self.sample_list, self, Type=cn_type, chrN=self.chrN, arm=self.arm,
                                                copy_number=str(cn_a1) + '/' + str(cn_a2), allelic_cn=cn_a1,
                                                supporting_muts=self.supporting_muts, ccf_hat=ccf_hat))
        if self.cn_a2 != baseline:
            cn_type = cn_type_prefix + ('_gain' if self.cn_a2 > baseline else '_loss')
            ccf_hat = min(self.cn_a2 - baseline, 1.) if self.cn_a2 > baseline else min(baseline - self.cn_a2, 1.)
            self.cn_events.append(TimingCNEvent(self.sample_list, self, Type=cn_type, chrN=self.chrN, arm=self.arm,
                                                copy_number=str(cn_a1) + '/' + str(cn_a2), allelic_cn=cn_a2,
                                                supporting_muts=self.supporting_muts, ccf_hat=ccf_hat))


class TimingCNEvent(object):
    def __init__(self, sample_list, state, Type=None, chrN=None, arm=None, pi_dist=None, copy_number=None,
                 allelic_cn=None, supporting_muts=None, is_clonal=None, cluster_id=None,
                 cn_state_whitelist=Enums.cn_state_whitelist, ccf_hat=None, start=None, end=None, band_range_str=None,
                 timing_status='unknown', timing_source=None, timing_reason=None):
        self.sample_list = sample_list
        self.state = state
        self.Type = Type
        self.chrN = chrN
        self.arm = arm
        self.pi_dist = pi_dist
        self.timing_info = []
        self.copy_number = copy_number
        self.cn_a1, self.cn_a2 = map(float, self.copy_number.split('/'))
        self.allelic_cn = allelic_cn
        self.total_cn = self.cn_a1 + self.cn_a2
        self.supporting_muts = supporting_muts if (supporting_muts is not None) else []
        self.is_clonal = is_clonal
        self.cluster_id = cluster_id
        self.cn_state_whitelist = cn_state_whitelist
        self.gain = None
        self.ccf_hat = ccf_hat
        self.start = start
        self.end = end
        self.band_range_str = band_range_str

        self.timing_status = timing_status
        self.timing_source = timing_source
        self.timing_reason = timing_reason
        self.allele = None
        self.ccf_sd = None
        self.continuous_cn = None
        self.clonality_ambiguous = False
        self.high_level = False  # a focal amplification far above its arm, see Patient._continuous_focal_call
        self.rule_side = None  # 'before' or 'after' the WGD, when the WGD rule placed the event there
        self.wgd_side = None  # rule_side, kept with clonal_unresolved: two events on one side cannot be ordered


    def __repr__(self):
        if self.Type.startswith('Focal_') and self.start is not None and self.end is not None:
            # return '<TimingCNEvent object: {}_{}{}:{}-{}>'.format(self.Type, self.chrN, self.arm, self.start, self.end)
            return '<TimingCNEvent object: {}_{}{}.{}>'.format(self.Type, self.chrN, self.arm, self.band_range_str)
        return '<TimingCNEvent object: {}_{}{}>'.format(self.Type, self.chrN, self.arm)

    @property
    def has_pi_posterior(self):
        return self.pi_dist is not None and np.sum(self.pi_dist) > 0

    def get_effective_pi_dist(self, fallback=None):
        if self.pi_dist is not None:
            return self.pi_dist
        return fallback

    @property
    def is_time_estimable(self):
        return self.timing_status == 'estimated'

    def set_not_estimable(self, reason='insufficient_local_support'):
        self.pi_dist = None
        self.timing_status = 'not_estimable'
        self.timing_source = None
        self.timing_reason = reason

    def set_prior_only(self, pi_dist, source='arm_prior', reason='borrowed_prior_only'):
        self.pi_dist = pi_dist / np.sum(pi_dist)
        self.timing_status = 'prior_only'
        self.timing_source = source
        self.timing_reason = reason

    @property
    def event_name(self):
        if self.Type.startswith('Arm_'):
            return self.Type[4:] + '_' + self.chrN + self.arm
        if self.Type.startswith('Chrom_'):
            return self.Type[6:] + '_' + self.chrN
        if self.Type.startswith('Focal_'):
            # Use coordinate-based naming to avoid key collisions when multiple
            # focal intervals map to the same cytoband(s).
            if self.start is not None and self.end is not None:
                coord = f'{self.chrN}:{self.start}-{self.end}'
            elif self.band_range_str is not None:
                coord = f'{self.chrN}{self.arm}.{self.band_range_str}'
            else:
                coord = f'{self.chrN}{self.arm}'
            return f'{self.Type}_{coord}'
        raise NotImplementedError('ONLY ARM AND FOCAL EVENTS CURRENTLY SUPPORTED')

    @property
    def log_p2_prior(self):
        if self.cn_a2 == 2.:
            if self.cn_a1 == 0. or self.cn_a1 == 2.:
                return np.log(2 * (np.linspace(0, 1, 101) + 1) ** -2)
            if self.cn_a1 == 1.:
                return np.log(3 * (np.linspace(0, .5, 101) + 1) ** -2)
        raise NotImplementedError('Higher gains not implemented')

    def get_p2_dist_for_gain(self):
        """
        get p2 distribution for a copy number gain
        """
        log_p2_posterior = np.zeros(101)
        for mut in self.supporting_muts:  # multiply mutation p2 distributions to get cnv p2 distribution
            mut_p2_dist = np.zeros(101)
            if (mut.local_cn_a1 == self.cn_a1).all() and (mut.local_cn_a2 == self.cn_a2).all() and mut.is_clonal:
                if mut.log_mult_dist is None:
                    mut.get_mult_dist()
                mut_p2_dist += np.sum(mut.log_mult_dist, 1)
            log_p2_posterior += mut_p2_dist - logsumexp(mut_p2_dist)
        log_p2_posterior = log_p2_posterior - logsumexp(log_p2_posterior) + self.log_p2_prior
        p2_real = self._correct_p2(log_p2_posterior - logsumexp(log_p2_posterior))
        return p2_real

    def get_pi_dist_for_gain(self):
        """
        get pi distribution for a copy number gain from p2 distribution
        """
        if (self.cn_a1, self.cn_a2) not in self.cn_state_whitelist:
            return
        p2_real = self.get_p2_dist_for_gain()
        # p2_interp = scipy.interpolate.PchipInterpolator(np.linspace(0, 1, 101), p2_real)
        # pi_domain = np.linspace(0, 1, 101)
        # change of variables
        #TODO: make this more readable
        if self.cn_a1 == 0. or self.cn_a1 == 2.:
            p2_domain_in_pi_space = 2 * np.linspace(0, 1, 101) / (np.linspace(0, 1, 101) + 1)
            pi_dist = np.zeros(101)
            for p2, p2_diff, prob in zip(p2_domain_in_pi_space, np.diff(p2_domain_in_pi_space), p2_real):
                min_ = p2 * 100.
                max_ = (p2 + p2_diff) * 100.
                bins = np.arange(np.floor(min_), np.floor(max_ + 1 if max_.is_integer() else max_ + 2)) # affected bins
                proportions = np.diff(np.clip(bins, min_, max_)) / (p2_diff * 100.)
                pi_dist[bins[:-1].astype(int)] += proportions * prob
            # pi_domain_in_p2_space = pi_domain / (2 - pi_domain)
            # pi_dist = p2_interp(pi_domain_in_p2_space) * 2 / ((2 - pi_domain) ** 2)
        else:
            p2_domain = np.linspace(0, .5, 101)
            p2_domain_in_pi_space = 3 * p2_domain / (p2_domain + 1)
            pi_dist = np.zeros(101)
            p2_linear_interpolator = scipy.interpolate.interp1d(np.linspace(0, 1, 101), p2_real)
            p2_trapezoid = p2_linear_interpolator(p2_domain)
            for p2, p2_diff, prob in zip(p2_domain_in_pi_space, np.diff(p2_domain_in_pi_space), p2_trapezoid):
                min_ = p2 * 100.
                max_ = (p2 + p2_diff) * 100.
                bins = np.arange(np.floor(min_), np.floor(max_ + 1 if max_.is_integer() else max_ + 2))
                proportions = np.diff(np.clip(bins, min_, max_)) / (p2_diff * 100.)
                pi_dist[bins[:-1].astype(int)] += proportions * prob
            # pi_domain_in_p2_space = pi_domain / (3 - pi_domain)
            # pi_dist = p2_interp(pi_domain_in_p2_space) * 3 / ((3 - pi_domain) ** 2)
        pi_dist[100] = 2 * pi_dist[99] - pi_dist[98]
        np.clip(pi_dist, 1e-20, None, out=pi_dist)
        self.pi_dist = pi_dist / sum(pi_dist)

    def _correct_p2(self, log_p2_posterior, polyfit_degree=3):
        """
        Function to correct p2 distribution for detection limit and clonality detection
        """
        full_p2_domain = np.linspace(0, 1, 101)
        p2_simulated = self._simulate_p2() + 1e-20
        p2_simulated_coeffs = np.polyfit(full_p2_domain, p2_simulated, polyfit_degree)
        p2_simulated_poly = np.poly1d(p2_simulated_coeffs)
        p2_transform_derivative = np.poly1d(
            [coeff * (polyfit_degree - i) for i, coeff in enumerate(p2_simulated_coeffs[:-1])])
        p2_real_interp = scipy.interpolate.PchipInterpolator(full_p2_domain, np.exp(log_p2_posterior))
        # Probability density change of variables formula
        p2_real = p2_real_interp(p2_simulated_poly(full_p2_domain)) * p2_transform_derivative(full_p2_domain) + 1e-20
        # p2_simulated_fit = np.clip(p2_simulated_poly(full_p2_domain), 0., 1.)
        # p2_posterior = np.exp(log_p2_posterior)
        # p2_real = np.zeros(101)
        # for p2, p2_diff, prob in zip(p2_simulated_fit, np.diff(p2_simulated_fit), p2_posterior):
        #     min_ = p2 * 100.
        #     max_ = (p2 + p2_diff) * 100.
        #     bins = np.arange(np.floor(min_), np.floor(max_ + 1 if max_.is_integer() else max_ + 2))
        #     proportions = np.diff(np.clip(bins, min_, max_)) / (p2_diff * 100.)
        #     p2_real[bins[:-1].astype(int)] += proportions * prob
        return p2_real

    def _simulate_p2(self, n_iter=100):
        """
        simulate n_iter mutations from ND coverage tracks to get a corrected p2 space
        """
        coverage_list = []
        muts_list = list(self.supporting_muts)
        shuffle, binomial = random.shuffle, np.random.binomial
        seed = getattr(getattr(self.state, 'options', None), 'seed', None)
        if seed is not None:
            # Supporting mutations are sorted first, so the draw depends on the seed and the event only.
            muts_list.sort(key=lambda mut: mut.var_str)
            key = [seed, self.Type, self.chrN, self.arm, self.copy_number, self.allelic_cn, self.start, self.end]
            key += [mut.var_str for mut in muts_list] + [getattr(s, 'sample_name', '') for s in self.sample_list]
            key = '|'.join(map(str, key))
            state = zlib.crc32(key.encode())
            shuffle, binomial = random.Random(state).shuffle, np.random.RandomState(state).binomial
        while len(coverage_list) < n_iter:
            shuffle(muts_list)
            for mut in muts_list:
                coverage = mut.ref_cnt + mut.alt_cnt
                if hasattr(coverage, '__iter__'):
                    coverage_list.append(np.asarray(coverage, dtype=int))
                else:
                    coverage_list.append(int(coverage))
        purity = np.array(self.state.purity)
        germline = 1 if self.state.haploid else 2
        p2_simulated = []
        for p2 in np.linspace(0, 1, 101):
            n_m1_muts = round((1 - p2) * n_iter)
            n_detected_by_mult = {1: 0., 2: 0.}
            for i in range(n_iter):
                cov = coverage_list[i]
                cov = np.asarray(cov, dtype=float)
                mult = 1 if i < n_m1_muts else 2
                expected_af = mult * purity / (germline * (1 - purity) + self.total_cn * purity)
                alt_count = binomial(cov.astype(int), expected_af).astype(float)
                af_mode = np.divide(alt_count, cov, out=np.zeros_like(alt_count, dtype=float), where=(cov > 0))
                ccf_mode = af_mode * (germline * (1 - purity) + purity * self.total_cn) / purity
                if all(ccf_mode >= .86) and all(alt_count >= 3):
                    n_detected_by_mult[mult] += 1.
            if n_detected_by_mult[1] == 0. and n_detected_by_mult[2] == 0.:
                p2_simulated.append(0.)
            else:
                p2_simulated.append(n_detected_by_mult[2] / (n_detected_by_mult[1] + n_detected_by_mult[2]))
        return np.array(p2_simulated)

    def get_pi_dist_for_loss(self, WGD):
        """
        get pi distribution for a copy number loss by comparing to WGD event
        """
        if WGD is None:
            self.pi_dist = None
            return
        if self.allelic_cn == 0.:  # loss before WGD
            pi_dist = 1 - np.cumsum(WGD.pi_dist)
            self.rule_side = 'before'
        elif self.allelic_cn == 1.:  # loss after WGD
            pi_dist = np.cumsum(WGD.pi_dist)
            self.rule_side = 'after'
        else:
            return
        self.pi_dist = pi_dist / sum(pi_dist)

    def get_pi_dist_for_higher_gain(self, WGD):
        """
        get pi distribution for a copy number gain on top of a WGD event
        """
        if WGD is None:
            self.pi_dist = None
            return
        if self.allelic_cn == 3.:  # loss after WGD
            pi_dist = np.cumsum(WGD.pi_dist)
            self.rule_side = 'after'
        elif self.allelic_cn >= 4.:  # loss before WGD
            pi_dist = 1 - np.cumsum(WGD.pi_dist)
            self.rule_side = 'before'
        else:
            return
        self.pi_dist = pi_dist / sum(pi_dist)

    def log_timing_info(self, msg):
        if msg not in self.timing_info:
            self.timing_info += [msg]


class TimingMut(object):
    """
    Class for storing information on snps and indels and for getting multiplicity likelihood distributions
    """
    def __init__(self, sample_list, gene, chrN, pos, alt, ref, alt_cnt, ref_cnt, local_cn_a1, local_cn_a2, ccf_dist,
                 prot_change=None, pi_dist=None, is_clonal=None, clonal_cutoff=.86, cluster_assignment=None):
        self.sample_list = sample_list
        self.gene = gene
        self.chrN = chrN
        self.pos = pos
        self.alt = alt
        self.ref = ref
        self.alt_cnt = alt_cnt
        self.ref_cnt = ref_cnt
        self.local_cn_a1 = local_cn_a1
        self.local_cn_a2 = local_cn_a2
        self.ccf_dist = ccf_dist
        self.cluster_assignment = cluster_assignment
        bins = np.linspace(0, 1, 101)
        ccf_hat = sum(bins * self.ccf_dist)
        if is_clonal is None:
            if self.cluster_assignment is None:
                self.is_clonal = all(ccf_hat > clonal_cutoff) if isinstance(ccf_hat, np.ndarray) else ccf_hat > clonal_cutoff
            else:
                self.is_clonal = self.cluster_assignment == 1
        else:
            self.is_clonal = is_clonal
        self.prot_change = prot_change
        self.pi_dist = pi_dist
        self.mult_lik_dict = {}
        self.log_mult_dist = None
        self.timing_status = None

    def __eq__(self, other):
        if not isinstance(other, TimingMut):
            return False
        return self.var_str == other.var_str

    @property
    def germline_copies(self):
        """Copies of this locus in a normal cell: 1 on a hemizygous chromosome when the sex is known."""
        options = getattr(self.sample_list[0], 'options', None)
        return 1 if options is not None and options.is_hemizygous(self.chrN) else 2

    def __hash__(self):
        return hash((self.var_str, (sample.sample_name for sample in self.sample_list)))

    def __repr__(self):
        return self.var_str

    @property
    def var_str(self):
        return ':'.join(map(str, [self.chrN, self.pos, self.ref, self.alt]))

    @property
    def event_name(self):
        if self.gene and self.prot_change:
            return self.gene + '_' + self.prot_change
        return self.var_str

    def get_multiplicity_likelihoods(self):
        """
        Likelihood = binomial probability of alt count from purity, copy number, and multiplicity
        """
        mult_lik_dict = {}
        n = self.alt_cnt + self.ref_cnt
        k = self.alt_cnt
        purity = np.array([sample.purity for sample in self.sample_list])  # ND mutation
        for m in range(1, int(max(self.local_cn_a2)) + 1):
            # expected allele fraction from multiplicity, purity, and copy number
            p = (m * purity) / ((self.local_cn_a1 + self.local_cn_a2) * purity + (self.germline_copies * (1 - purity)))
            mult_lik_dict[m] = scipy.stats.binom.pmf(k, n, p)
        return mult_lik_dict

    def get_log_mult_1_2_distribution(self):
        """
        p2 distribution for an individual mutation
        """
        mult = np.linspace(self.mult_lik_dict[1], self.mult_lik_dict[2], 101)
        mult = np.nan_to_num(mult, nan=0.0, posinf=0.0, neginf=0.0)
        colsum = mult.sum(axis=0)
        probs = np.empty_like(mult)
        # normalize columns with mass; assume uniform for empty columns
        ok = np.isfinite(colsum) & (colsum > 0)
        probs[:, ok] = mult[:, ok] / colsum[ok]
        probs[:, ~ok] = 1.0 / mult.shape[0]
        probs = np.clip(probs, 1e-300, 1.0)
        return np.log(probs)

    def get_pi_dist(self, matched_gain):
        """
        Compute the timing distribution for this mutation relative to *matched_gain*.

        Per-sample approach: each sample is evaluated independently based on its
        own local copy-number state.  Samples whose CN state does not permit
        multiplicity decomposition are skipped.  The per-sample likelihoods
        (p_before, p_after the gain) are combined across all timeable samples
        via the same log-likelihood product used by the original implementation.

        This replaces the previous requirement that local CN must be identical
        across all samples, which dropped mutations whenever segmental CN
        estimates differed due to noise or subclonal CN changes.
        """
        if matched_gain.pi_dist is None:
            return

        n_samples = len(self.sample_list)
        purity = np.array([s.purity for s in self.sample_list])

        # Per-sample p(before) and p(after) the gain, based on multiplicity.
        p_before = np.full(n_samples, np.nan)
        p_after = np.full(n_samples, np.nan)

        for s in range(n_samples):
            a1_raw = self.local_cn_a1[s]
            a2_raw = self.local_cn_a2[s]
            if np.isnan(a1_raw) or np.isnan(a2_raw):
                continue
            a1 = round(a1_raw)
            a2 = round(a2_raw)
            max_cn = int(a2)
            if max_cn < 2:
                continue  # diploid or haploid — no gain to time against

            # Per-sample multiplicity likelihoods (scalar per multiplicity m)
            n_s = self.alt_cnt[s] + self.ref_cnt[s]
            k_s = self.alt_cnt[s]
            p_s = purity[s]
            total_cn = a1 + a2
            mult_lik = {}
            for m in range(1, max_cn + 1):
                p_m = (m * p_s) / (total_cn * p_s + float(self.germline_copies) * (1.0 - p_s))
                mult_lik[m] = scipy.stats.binom.pmf(k_s, n_s, p_m)

            if a1 == 0 and a2 >= 2 or a1 == a2 != 1:
                lb = sum(mult_lik[i] * i / max_cn for i in range(2, max_cn + 1))
                la = sum(mult_lik[i] * (max_cn - i) / max_cn for i in range(2, max_cn)) + mult_lik[1]
            elif a1 == 1 and a2 >= 2:
                lb = sum(mult_lik[i] * i / max_cn for i in range(1, max_cn + 1))
                la = sum(mult_lik[i] * (max_cn - i) / max_cn for i in range(1, max_cn))
            else:
                continue  # CN state not suitable for timing

            total = lb + la + 1e-12
            p_before[s] = lb / total
            p_after[s] = la / total

        mask = ~np.isnan(p_before)
        if not mask.any():
            return

        pi_cdf = np.cumsum(matched_gain.pi_dist)
        pi_dist = np.outer(1.0 - pi_cdf, p_before[mask]) + np.outer(pi_cdf, p_after[mask])
        pi_dist = pi_dist * (pi_dist >= 0) + 1e-10
        pi_dist = np.sum(np.log(pi_dist / np.sum(pi_dist, 0)), 1)
        self.pi_dist = np.exp(pi_dist - logsumexp(pi_dist))

    def get_mult_dist(self):
        cn_a1 = np.unique(self.local_cn_a1)
        cn_a2 = np.unique(self.local_cn_a2)
        if len(cn_a1) > 1 or len(cn_a2) > 1 or np.isnan(cn_a1[0]) or np.isnan(cn_a2[0]):
            return
        cn_a1 = cn_a1[0]
        cn_a2 = cn_a2[0]
        self.mult_lik_dict = self.get_multiplicity_likelihoods()
        if cn_a1 in (0., 1., 2.) and cn_a2 == 2.:
            self.log_mult_dist = self.get_log_mult_1_2_distribution()
