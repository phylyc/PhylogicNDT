from collections import Counter, defaultdict
import numpy as np
import itertools
from intervaltree import IntervalTree, Interval
import scipy.stats
import scipy.interpolate
from scipy.special import logsumexp
import random
import data.Enums as Enums


class TimingEngine(object):
    """
    Class for holding samples
    """
    def __init__(self, patient, cn_state_whitelist=Enums.cn_state_whitelist, min_supporting_muts=3):
        self.patient = patient
        self.genome = patient.genome
        self.cn_state_whitelist = cn_state_whitelist
        self.arm_regions = list(itertools.product(self.genome.CHROMS, self.genome.ARMS))
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
        self.get_focal_cn_events()

    def get_concordant_cn_states(self):
        """
        get cn states that are present in all samples (same chromsome arm and copy number) and create multisample cn_state instances
        """
        n_samples = len(self.sample_list)
        purity = [sample.purity for sample in self.sample_list]
        self.concordant_cn_states = {}
        for chrN, arm in self.arm_regions:
            region = chrN + arm
            states_across_samples = {}  # hold union of all copy number states in all samples
            for sample in self.sample_list:
                if region in sample.cn_states and len(sample.cn_states[region].cn_events):
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
                    multisample_state = TimingCNState(self.sample_list, chrN, arm, state, purity,
                                                      supporting_muts=supporting_muts)
                    self.concordant_cn_states[chrN + arm] = multisample_state

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
                for sample in self.sample_list:
                    s_state = sample.cn_states.get(cn_state)
                    if s_state is None:
                        ccf_hats.append(0.0)
                        continue
                    # find the matching event by Type; default to 0 if missing
                    match = next((e for e in s_state.cn_events if e.Type == etype), None)
                    ccf_hats.append(float(getattr(match, 'ccf_hat', 0.0) or 0.0))
                ccf_idx = np.clip(np.round(np.array(ccf_hats) * 100).astype(int), 0, 100)
                clonal_concordance = np.prod(cluster_ccfs[1][sample_idx, ccf_idx])
                is_clonal = True
                for c in cluster_ccfs:
                    if np.prod(cluster_ccfs[c][sample_idx, ccf_idx]) > clonal_concordance:
                        is_clonal = False
                        break
                eve.is_clonal = is_clonal
                self.all_cn_events[eve.event_name].append(eve)
                if is_clonal:
                    self.truncal_cn_events[eve.event_name].append(eve)

    def get_focal_cn_events(self):
        """
        Attach focal CNV events (produced by Patient.add_focal_cn_events) as TimingCNEvent instances.
        Minimal: build a consensus (cn_a1, cn_a2) across samples in the window; clonality from focal ccf_hat.
        """
        n_samples = len(self.sample_list)
        if n_samples == 0:
            return

        # --- cytoband cache: (chrom, band_name) -> (start_bp, end_bp)
        _band_index = {}
        _bands_by_chr = {}
        _band_index_ready = False

        def _ensure_band_index():
            nonlocal _band_index_ready
            if _band_index_ready:
                return
            try:
                for _, (chrom, start, end, band, *_) in self.genome.cytoband_table.iterrows():
                    _band_index[(chrom.lstrip('chr'), band)] = (int(start), int(end))
                # build per-chrom sorted lists for bp→band lookup
                tmp = defaultdict(list)
                for (c, b), (s, e) in _band_index.items():
                    tmp[c].append((s, e, b))
                for c in tmp:
                    _bands_by_chr[c] = sorted(tmp[c], key=lambda t: t[0])
                _band_index_ready = True
            except Exception as e:
                print(f'Could not load cytoBand.txt for band to bp conversion: {e}')
                _band_index_ready = True  # avoid retry storms

        def _band_to_bp(chrom, band_obj):
            """Return (start_bp, end_bp) for a Cytoband or band-name string; else None."""
            if band_obj is None:
                return None
            name = getattr(band_obj, 'band', None) or str(band_obj)
            _ensure_band_index()
            return _band_index.get((str(chrom).lstrip('chr'), name))

        def _band_at_bp(chrom, pos):
            arr = _bands_by_chr.get(str(chrom).lstrip('chr'))
            if not arr:
                return None
            for s, e, b in arr:
                if s <= pos < e:
                    return b
            return None

        def _band_range_str(chrom, s_bp, e_bp):
            # inclusive range label using bands at start and end-1
            b1 = _band_at_bp(chrom, s_bp)
            b2 = _band_at_bp(chrom, max(s_bp, e_bp - 1))
            if not b1 or not b2:
                return None
            return b1 if b1 == b2 else f'{b1}-{b2}'

        def _coords(chrom, start, end):
            """
            Normalize event bounds to integer bps.
            Accepts ints or Cytoband/string; returns (start_bp, end_bp) or None if unknown.
            """
            # ints already? just return
            if isinstance(start, (int, np.integer)) and isinstance(end, (int, np.integer)):
                return int(start), int(end)
            # try cytoband → bp
            s = _band_to_bp(chrom, start)
            e = _band_to_bp(chrom, end)
            if s is None or e is None:
                return None
            return s[0], e[1]

        # group focal events across samples
        groups = {}  # (chrom,start,end,arm,is_a1,category) -> {ccf: vec, cn_pairs: list, seg_bounds: list}
        for si, ts in enumerate(self.sample_list):
            for cn in ts.sample.low_coverage_mutations.values():
                if getattr(cn, 'type', None) == 'CNV' and cn.cn_category.startswith('Focal_'):
                    bounds = _coords(cn.chrN, cn.start, cn.end)
                    if bounds is None:
                        print(f'Skipping focal CNV without resolvable coords: {cn}')
                        continue
                    s_bp, e_bp = bounds
                    key = (str(cn.chrN), s_bp, e_bp, cn.arm, cn.a1, cn.cn_category)
                    rec = groups.setdefault(key, {'ccf': np.zeros(n_samples), 'cn_pairs': [None]*n_samples, 'seg_bounds': []})
                    rec['ccf'][si] = cn.ccf_hat
                    cn_a1 = np.nan
                    cn_a2 = np.nan
                    try:
                        ov = list(ts.CnProfile[str(cn.chrN)].overlap(s_bp, e_bp))
                    except Exception:
                        ov = []
                    if ov:
                        def _ovlen(iv):
                            return max(0, min(iv.end, e_bp) - max(iv.begin, s_bp))

                        def _totcn(iv):
                            seg = iv.data[1]
                            return seg['cn_a1'] + seg['cn_a2']

                        is_amp = ('gain' in cn.cn_category.lower()) or ('amp' in cn.cn_category.lower())
                        # prefer extreme total CN; break ties by larger overlap
                        if is_amp:
                            best = max(ov, key=lambda iv: (_totcn(iv), _ovlen(iv)))
                        else:
                            best = min(ov, key=lambda iv: (_totcn(iv), -_ovlen(iv)))
                        seg = best.data[1]
                        cn_a1 = float(seg['cn_a1'])
                        cn_a2 = float(seg['cn_a2'])
                        rec['seg_bounds'].append((best.begin, best.end))
                    rec['cn_pairs'][si] = (cn_a1, cn_a2)

        if not groups:
            return

        cluster_ccfs = self._get_cluster_ccfs()
        sample_idx = range(n_samples)

        def consensus_pair(pairs):
            vals = []
            for a1, a2 in pairs:
                if a1 is None or a2 is None or np.isnan(a1) or np.isnan(a2):
                    continue
                a1r = round(a1) if abs(a1 - round(a1)) <= 0.2 else a1
                a2r = round(a2) if abs(a2 - round(a2)) <= 0.2 else a2
                vals.append((a1r, a2r))
            return Counter(vals).most_common(1)[0][0] if vals else None

        # compare up to ordering and small rounding error
        def _pair_match(m_a1, m_a2, r_a1, r_a2, tol=0.25):
            x = sorted([float(m_a1), float(m_a2)])
            y = sorted([float(r_a1), float(r_a2)])
            return abs(x[0] - y[0]) <= tol and abs(x[1] - y[1]) <= tol

        for (chrom, start, end, arm, is_a1, category), rec in groups.items():
            cn_state = consensus_pair(rec['cn_pairs'])
            if cn_state is None:
                continue
            # gather supporting mutations in the UNION of best-overlap segments (more SNVs, same CN state)
            if rec['seg_bounds']:
                sup_start = min(b for b, _ in rec['seg_bounds'])
                sup_end = max(e for _, e in rec['seg_bounds'])
            else:
                sup_start, sup_end = start, end
            # IMPORTANT: match SNVs to the *per-sample* CN pair, not the cross-sample consensus
            supporting_muts = []
            for si, ts in enumerate(self.sample_list):
                samp_pair = rec['cn_pairs'][si]
                if not samp_pair or any(map(lambda x: x is None or np.isnan(x), samp_pair)):
                    continue
                a1_req, a2_req = samp_pair
                tree = ts.mutation_intervaltree.get(chrom, IntervalTree())
                for iv in tree[sup_start:sup_end]:
                    mut = iv.data  # TimingMut (this sample)
                    # mut.local_cn_a1 / a2 are per-sample scalars here
                    # if float(mut.local_cn_a1) == a1_req and float(mut.local_cn_a2) == a2_req:
                    if _pair_match(mut.local_cn_a1, mut.local_cn_a2, a1_req, a2_req):
                        supporting_muts.append(mut)
            purity = [s.purity for s in self.sample_list]
            state = TimingCNState(self.sample_list, chrom, arm, cn_state, purity, supporting_muts=supporting_muts)
            state.call_events(cn_type_prefix='Focal', wgd=(self.WGD is not None))
            # select allele-matched focal event
            eve = None
            for e in state.cn_events:
                if e.Type == category and ((is_a1 and e.allelic_cn == cn_state[0]) or (not is_a1 and e.allelic_cn == cn_state[1])):
                    eve = e
                    break
            if eve is None:
                continue
            eve.start = int(start)
            eve.end = int(end)
            is_amp = ('gain' in category.lower()) or ('amp' in category.lower())
            # borrow SNV CCF evidence from other AMP events with same CN state on same arm
            if len(supporting_muts) < self.min_supporting_muts and is_amp:
                borrowed = []
                # collect from arm region with the same (per-sample) CN pair
                for si, ts in enumerate(self.sample_list):
                    pair = rec['cn_pairs'][si]
                    if not pair:
                        continue
                    arm_start = 0 if arm=='p' else self.genome.CENT_DICT[chrom]
                    arm_end   = self.genome.CENT_DICT[chrom] if arm=='p' else self.genome.CHROM_DICT[chrom]
                    for iv in ts.mutation_intervaltree.get(chrom, IntervalTree())[arm_start:arm_end]:
                        mut = iv.data
                        if _pair_match(mut.local_cn_a1, mut.local_cn_a2, *pair):
                            borrowed.append(mut)
                            if len(borrowed) >= self.min_supporting_muts:
                                break
                supporting_muts.extend(borrowed)

            # focal-specific fallback: if low evidence, avoid collapsing to "after everything"
            if len(supporting_muts) < self.min_supporting_muts:
                eve.pi_dist = np.ones(101, dtype=float) / 101
            # set cytoband-based label for naming
            eve.band_range_str = _band_range_str(chrom, int(start), int(end))
            # clonality from focal ccf_hat vector
            ccf_idx = np.rint(np.clip(np.nan_to_num(np.asarray(rec['ccf'], float), nan=0.0, posinf=1.0, neginf=0.0),0.0, 1.0) * 100).astype(int)
            clonal_conc = np.prod(cluster_ccfs.get(1, np.ones((n_samples, 101)))[sample_idx, ccf_idx])
            is_clonal = True
            for c, grid in cluster_ccfs.items():
                if np.prod(grid[sample_idx, ccf_idx]) > clonal_conc:
                    is_clonal = False
                    break
            eve.is_clonal = is_clonal
            self.all_cn_events.setdefault(eve.event_name, []).append(eve)
            if is_clonal:
                self.truncal_cn_events.setdefault(eve.event_name, []).append(eve)

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
        time gains, then supporting mutations, then WGD, then supporting mutations, then higher gains and losses
        """
        uniform_dist = np.ones(101) / 101.
        subclonal_dist = np.zeros(101)
        subclonal_dist[100] = 1.
        if self.WGD is not None:
            self.WGD.get_pi_dist()
            for mut in self.mutations.values():
                if mut.is_clonal:
                    mut.get_pi_dist(self.WGD)
            for cn_event_name in self.all_cn_events:
                for cn_event in self.all_cn_events[cn_event_name]:
                    if not cn_event.is_clonal:
                        cn_event.pi_dist = subclonal_dist
                    elif cn_event.Type.endswith('gain'):
                        cn_event.get_pi_dist_for_higher_gain(self.WGD)
                    elif cn_event.Type.endswith('loss'):
                        cn_event.get_pi_dist_for_loss(self.WGD)
        else:
            for cn_event_name, evs in self.all_cn_events.items():
                for cn_event in evs:
                    is_focal_gain = cn_event.Type.startswith('Focal_gain')
                    have_support = len(cn_event.supporting_muts) >= self.min_supporting_muts
                    state_ok = ((cn_event.cn_a1, cn_event.cn_a2) in self.cn_state_whitelist)
                    if not cn_event.is_clonal:
                        # focal amp events with weak clonality get uninformative, not "always-after"
                        cn_event.pi_dist = uniform_dist if is_focal_gain else subclonal_dist
                    elif "gain" in cn_event.Type and have_support and state_ok:
                        cn_event.get_pi_dist_for_gain()
                        for mut in cn_event.supporting_muts:
                            if mut.is_clonal:
                                mut.get_pi_dist(cn_event)
                    else:
                        cn_event.pi_dist = uniform_dist
        for mut in self.mutations.values():
            if mut.pi_dist is None:
                if mut.is_clonal:
                    mut.pi_dist = uniform_dist
                else:
                    mut.pi_dist = subclonal_dist


class TimingSample(object):
    """
    class for holding a sample
    """
    def __init__(self, sample, engine, cn_state_whitelist=Enums.cn_state_whitelist):
        self.sample = sample
        self.sample_name = self.sample.sample_name
        self.engine = engine
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
            true_arm_bp = end - start
            arm_bp = 0
            state_bps = {}
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
            states_over_threshold = [cn for cn, bp in state_bps.items() if bp > arm_bp * size_threshold]
            if arm_bp < true_arm_bp * .5:  # Only call if at least 50% of arm is in segs
                self.missing_arms.append((chrN, arm))
            elif len(states_over_threshold) == 1:
                cn_a1, cn_a2 = next(iter(states_over_threshold))
                if concordant_states is None:
                    supporting_muts = [interval.data for interval in self.mutation_intervaltree[chrN][start:end]
                                       if interval.data.local_cn_a1 == cn_a1 and interval.data.local_cn_a2 == cn_a2]
                    cn_state = TimingCNState([self], chrN, arm, (cn_a1, cn_a2), [self.purity],
                                             supporting_muts=supporting_muts)
                    self.cn_states[chrN + arm] = cn_state
                elif chrN + arm in concordant_states:
                    supporting_muts = [interval.data for interval in self.concordant_mutation_intervaltree[chrN][start:end]
                                      if interval.data.local_cn_a1 == cn_a1 and interval.data.local_cn_a2 == cn_a2]
                    cn_state = TimingCNState([self], chrN, arm, (cn_a1, cn_a2), [self.purity],
                                             supporting_muts=supporting_muts)
                    self.concordant_cn_states[chrN + arm] = cn_state

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
                if cn_state.cn_a1 == 0 or cn_state.cn_a1 >= 2 and cn_state.cn_a2 >= 2:  # region supports WGD if 0/2 or 2/2
                    supporting_arm_states.append(TimingCNState([self], cn_state.chrN, cn_state.arm,
                        (cn_state.cn_a1, cn_state.cn_a2), cn_state.purity, supporting_muts=cn_state.supporting_muts))
            self.concordant_WGD = TimingWGD(supporting_arm_states=supporting_arm_states)
            for cn_state in self.concordant_cn_states.values():
                cn_state.call_events(wgd=True)
        else:
            for cn_state in self.cn_states.values():
                if cn_state.cn_a2 >= 2:
                    supporting_arm_states.append(cn_state)
                if (cn_state.cn_a1 == 0 or cn_state.cn_a1 >= 2) and cn_state.cn_a2 >= 2:  # region supports WGD if 0/2 or 2/2
                    regions_supporting_WGD.append(cn_state)
                if cn_state.cn_a1 >= 2 and cn_state.cn_a2 >= 2:
                    regions_both_arms_gained.append(cn_state)
            if len(regions_both_arms_gained) >= 5 and len(regions_supporting_WGD) * 2 >= \
                    len(self.arm_regions) - len(self.missing_arms):
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
                    self.concordant_cn_events[eve.event_name].append(eve)
                else:
                    self.cn_events[eve.event_name].append(eve)


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
        self.get_supporting_muts()

    def get_pi_dist(self):
        """
        Get pi dist by multiplying pi dists for supporting arm gains
        """
        pi_dist = np.zeros(101)
        for state in self.supporting_arm_states:
            for eve in state.cn_events:
                if eve.Type == 'Arm_gain':
                    if len(eve.supporting_muts) < self.min_timeable_muts:
                        continue
                    if eve.pi_dist is None:
                        eve.get_pi_dist_for_gain()
                    if eve.pi_dist is None:
                        continue
                    pi_dist += np.log(eve.pi_dist)
        self.pi_dist = np.exp(pi_dist - logsumexp(pi_dist))

    def get_supporting_muts(self):
        self.supporting_muts = []
        for state in self.supporting_arm_states:
            self.supporting_muts.extend(state.supporting_muts)


class TimingCNState(object):
    """
    Arm events are timed based on copy number state of both arms
    """

    def __init__(self, sample_list, chrN, arm, cn_state, purity, supporting_muts=None, cn_type_prefix='Arm'):
        self.sample_list = sample_list
        self.chrN = chrN
        self.arm = arm
        self.cn_a1, self.cn_a2 = cn_state
        self.purity = purity
        self.supporting_muts = supporting_muts
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

    def call_events(self, cn_type_prefix='Arm', wgd=False):
        """
        create cn event instances for events (gain for 1/2, gain and loss for 0/2, two gains for 2/2)
        """
        self.cn_events = []
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
                 cn_state_whitelist=Enums.cn_state_whitelist, ccf_hat=None, start=None, end=None, band_range_str=None):
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


    def __repr__(self):
        if self.Type.startswith('Focal_') and self.start is not None and self.end is not None:
            # return '<TimingCNEvent object: {}_{}{}:{}-{}>'.format(self.Type, self.chrN, self.arm, self.start, self.end)
            return '<TimingCNEvent object: {}_{}{}.{}>'.format(self.Type, self.chrN, self.arm, self.band_range_str)
        return '<TimingCNEvent object: {}_{}{}>'.format(self.Type, self.chrN, self.arm)

    @property
    def event_name(self):
        if self.Type.startswith('Arm_'):
            return self.Type[4:] + '_' + self.chrN + self.arm
        if self.Type.startswith('Focal_'):
            # Prefer cytoband naming if available; else fall back to coordinates
            if self.band_range_str is not None:
                coord = f'{self.chrN}{self.arm}.{self.band_range_str}'
            else:
                coord = (f'{self.chrN}{self.arm}:{self.start}-{self.end}'
                         if self.start is not None and self.end is not None
                         else f'{self.chrN}{self.arm}')
            return f'{self.Type}_{coord}'
        raise NotImplementedError('ONLY ARM AND FOCAL EVENTS CURRENTLY SUPPORTED')

    @property
    def log_p2_prior(self):
        if self.cn_a2 == 2.:
            if self.cn_a1 == 0. or self.cn_a2 == 2.:
                return np.log(3 * (np.linspace(0, 1, 101) + 1) ** -2)
            if self.cn_a1 == 1.:
                return np.log(2 * (np.linspace(0, .5, 101) + 1) ** -2)
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
        while len(coverage_list) < n_iter:
            random.shuffle(muts_list)
            for mut in muts_list:
                coverage = mut.ref_cnt + mut.alt_cnt
                if hasattr(coverage, '__iter__'):
                    coverage_list.append(np.asarray(coverage, dtype=int))
                else:
                    coverage_list.append(int(coverage))
        purity = np.array(self.state.purity)
        p2_simulated = []
        for p2 in np.linspace(0, 1, 101):
            n_m1_muts = round((1 - p2) * n_iter)
            n_detected_by_mult = {1: 0., 2: 0.}
            for i in range(n_iter):
                cov = coverage_list[i]
                cov = np.asarray(cov, dtype=float)
                mult = 1 if i < n_m1_muts else 2
                expected_af = mult * purity / (2 * (1 - purity) + self.total_cn * purity)
                alt_count = np.random.binomial(cov.astype(int), expected_af).astype(float)
                af_mode = np.divide(alt_count, cov, out=np.zeros_like(alt_count, dtype=float), where=(cov > 0))
                ccf_mode = af_mode * (2 * (1 - purity) + purity * self.total_cn) / purity
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
        elif self.allelic_cn == 1.:  # loss after WGD
            pi_dist = np.cumsum(WGD.pi_dist)
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
        elif self.allelic_cn >= 4.:  # loss before WGD
            pi_dist = 1 - np.cumsum(WGD.pi_dist)
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

    def __eq__(self, other):
        if not isinstance(other, TimingMut):
            return False
        return self.var_str == other.var_str

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
            p = (m * purity) / ((self.local_cn_a1 + self.local_cn_a2) * purity + (2 * (1 - purity)))
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
        get pi dist for mutation by comparing to a matched gain
        """
        if matched_gain.pi_dist is None:
            return
        cn_a1 = np.unique(self.local_cn_a1)
        cn_a2 = np.unique(self.local_cn_a2)
        if len(cn_a1) > 1 or len(cn_a2) > 1:
            return
        cn_a1 = cn_a1[0]
        cn_a2 = cn_a2[0]
        if not self.mult_lik_dict:
            self.get_mult_dist()
        if not self.mult_lik_dict:
            return
        if cn_a1 == 0. and cn_a2 >= 2. or cn_a1 == cn_a2 != 1.:
            max_cn = int(cn_a2)
            lik_before_gain = sum(self.mult_lik_dict[i] * i / max_cn for i in range(2, max_cn + 1))
            lik_after_gain = sum(self.mult_lik_dict[i] * (max_cn - i) / max_cn for i in range(2, max_cn)) + self.mult_lik_dict[1]
        elif cn_a1 == 1. and cn_a2 >= 2.:
            max_cn = int(cn_a2)
            lik_before_gain = sum(self.mult_lik_dict[i] * i / max_cn for i in range(1, max_cn + 1))
            lik_after_gain = sum(self.mult_lik_dict[i] * (max_cn - i) / max_cn for i in range(1, max_cn))
        else:
            return
        total_lik = lik_before_gain + lik_after_gain + 1e-12 # Add small episilon to avoid division by zero
        p_before_gain = lik_before_gain / total_lik
        p_after_gain = lik_after_gain / total_lik
        pi_cdf = np.cumsum(matched_gain.pi_dist)
        pi_dist = np.outer((1 - pi_cdf), p_before_gain) + np.outer(pi_cdf, p_after_gain)
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
