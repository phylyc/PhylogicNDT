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
                    }
                )
                rec['ccf'][si] = float(getattr(cn, 'ccf_hat', 0.0))
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

            eve.is_clonal = is_clonal

            # IMPORTANT:
            # do NOT assign uniform pi_dist here. Leave that decision to time_events().
            eve.pi_dist = None

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
            return ev.Type.startswith('Arm_')

        def _has_support(ev):
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
            if hasattr(self.WGD, 'timing_status'):
                self.WGD.timing_status = 'estimated'
                self.WGD.timing_source = 'wgd'
                self.WGD.timing_reason = None

            for mut in self.mutations.values():
                if mut.is_clonal:
                    mut.get_pi_dist(self.WGD)

            for evs in self.all_cn_events.values():
                for cn_event in evs:
                    if not cn_event.is_clonal:
                        cn_event.pi_dist = subclonal_dist.copy()
                        _set_status(cn_event, status='subclonal', source='clonality_rule', reason='event_not_truncal')
                    elif cn_event.Type.endswith('gain'):
                        cn_event.get_pi_dist_for_higher_gain(self.WGD)
                        if cn_event.pi_dist is not None:
                            _set_status(cn_event, status='estimated', source='wgd_rule', reason=None)
                        else:
                            _set_not_estimable(cn_event, reason='unsupported_gain_state')
                    elif cn_event.Type.endswith('loss'):
                        cn_event.get_pi_dist_for_loss(self.WGD)
                        if cn_event.pi_dist is not None:
                            _set_status(cn_event, status='estimated', source='wgd_rule', reason=None)
                        else:
                            _set_not_estimable(cn_event, reason='unsupported_loss_state')
                    else:
                        _set_not_estimable(cn_event, reason='unsupported_event_type')

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
                if not ev.is_clonal:
                    ev.pi_dist = subclonal_dist.copy()
                    _set_status(ev, status='subclonal', source='clonality_rule', reason='event_not_truncal')
                    continue

                if _time_gain_from_local(ev):
                    arm_prior_lookup[(ev.Type, ev.chrN, ev.arm, ev.cn_a1, ev.cn_a2)] = ev.pi_dist.copy()
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

                # No local evidence, but explicit borrowing enabled
                if borrowed_prior is not None:
                    _set_prior_only(ev, borrowed_prior, source='arm_prior', reason='borrowed_prior_only')
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
                else:
                    mut.pi_dist = subclonal_dist.copy()


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
