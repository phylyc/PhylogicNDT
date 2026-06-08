import logging
import numpy as np
import pandas as pd
import itertools
from data.Patient import Patient
from SinglePatientTiming import TimingEngine
from output.PhylogicOutput import PhylogicOutput


def run_tool(args):
    logging.debug('Arguments {}'.format(args))

    patient_data = Patient(artifact_blacklist=args.artifact_blacklist,
                           indiv_name=args.indiv_id, artifact_whitelist=args.artifact_whitelist,
                           driver_genes_file=args.driver_genes_file, genome_build=args.genome_build)

    if args.sif:  # if sif file is specified
        sif = pd.read_csv(args.sif, sep="\t")
        for _, row in sif.iterrows():
            patient_data.addSample(
                row['maf_fn'],
                row['sample_id'],
                input_type='post-clustering',
                # seg_input_type='timing_format',
                timepoint_value=row['timepoint'],
                grid_size=args.grid_size,
                _additional_muts=None,
                seg_file=row['seg_fn'],
                purity=row['purity']
            )

    else:  # if sample names/files are specified directly on cmdline

        # sort order on timepoint or order of entry on cmdline if not present
        print(args.sample_data)
        for idx, sample_entry in enumerate(args.sample_data):
            ##for now, assume input order is of the type sample_id\tmaf_fn\tseg_fn\tpurity\ttimepoint
            smpl_spec = sample_entry.strip('\n').split(':')
            sample_id = smpl_spec[0]
            maf_fn = smpl_spec[1]
            seg_fn = smpl_spec[2]
            purity = float(smpl_spec[3])
            timepoint = float(smpl_spec[4])

            patient_data.addSample(maf_fn, sample_id, input_type='post-clustering', # seg_input_type='timing_format',
                                   timepoint_value=timepoint, grid_size=args.grid_size,
                                   _additional_muts=None, seg_file=seg_fn, purity=purity)

    if args.intersect_cn_trees:
        patient_data.intersect_cn_trees()
    else:
        patient_data.get_arm_level_cn_events()
        if args.gistic_fn:
            patient_data.add_focal_cn_events(focal_regions=args.gistic_fn)
    patient_data.preprocess_samples()
    if args.min_supporting_muts < 1:
        raise ValueError('Invalid value for min_supporting_muts')
    timing_engine = TimingEngine.TimingEngine(patient_data, min_supporting_muts=args.min_supporting_muts)
    timing_engine.time_events()
    phylogicoutput = PhylogicOutput()
    phylogicoutput.write_timing_tsv(timing_engine)
    with open(args.driver_genes_file) as f:
        driver_genes = [line.strip() for line in f]
    comps = compare_events(timing_engine, drivers=driver_genes)
    phylogicoutput.write_comp_table(args.indiv_id, comps)
    phylogicoutput.generate_html_from_timing(args.indiv_id, timing_engine, comps, drivers=driver_genes)


def _dedup_cn_events(cn_events_dict):
    """Deduplicate CN events that share the same event_name.

    When a CN state produces two events with the same name (e.g. two
    Arm_gain events on the same arm from a 2/2 state), keep only the
    event with the most supporting mutations so comparisons are not
    silently overwritten downstream.
    """
    seen = {}  # event_name -> (best_event, n_supporting_muts)
    for evs in cn_events_dict.values():
        for ev in evs:
            name = ev.event_name
            n_sup = len(getattr(ev, 'supporting_muts', None) or [])
            prev = seen.get(name)
            if prev is None:
                seen[name] = (ev, n_sup)
            else:
                if n_sup > prev[1]:
                    seen[name] = (ev, n_sup)
    return [ev for ev, _ in seen.values()]


def compare_events(timing_engine, drivers=(), uncertainty_bins=5, uniform_atol=1e-6):
    """Compare event timings using pairwise probabilities from the joint
    distribution over discrete pi bins.

    For each pair of events, the joint probability matrix is computed as the
    outer product of their (independent) pi distributions.  The upper
    triangle (event1 earlier), lower triangle (event1 later), and a
    near-diagonal uncertainty band are summed to obtain the three
    comparison probabilities.

    Events whose timing cannot be estimated (``pi_dist is None``) are
    still included with fully-unknown comparisons ``(0, 0, 1)`` so that
    they remain visible in the comparison table and are counted correctly
    in downstream LeagueModel prevalence.

    Events whose ``pi_dist`` is effectively uniform (e.g., focal losses
    without WGD anchoring, focal homdels, or focal events that fell
    through the uniform fallback in time_events) are also recorded as
    fully-unknown.  A uniform pi_dist carries no locus-specific ordering
    signal: ``np.outer(uniform, anything)`` produces ~50/50 before/after
    splits that the LeagueModel would otherwise mistake for genuine
    pairwise evidence.  Reporting these honestly as ``(0, 0, 1)`` lets
    the LeagueModel's low-count rescue handle them appropriately.

    Note: the ``subclonal_dist`` (a delta at bin 100, "after everything")
    is *not* uniform and remains informative.

    When multiple CN events share the same ``event_name`` (e.g. two
    gains from a 2/2 arm state), only the event with the most supporting
    mutations is kept to avoid silent overwrites in the output dict.

    Args:
        timing_engine: A ``TimingEngine`` instance with timed events.
        drivers: Iterable of driver gene names.
        uncertainty_bins: Half-width of the near-diagonal band (in pi
            bins) that is counted as *unknown*.  Default 5 (≈ 5 %).
        uniform_atol: Absolute tolerance for declaring a pi_dist
            uniform.  A distribution is uniform if every bin is within
            ``uniform_atol`` of the mean after normalization.

    Returns:
        dict mapping ``(event1_name, event2_name)`` to
        ``(p_before, p_after, p_unknown)``.
    """
    def _is_uninformative(pi_dist):
        """True if pi_dist is None, empty, zero-sum, or effectively uniform."""
        if pi_dist is None:
            return True
        p = np.asarray(pi_dist, dtype=float)
        if p.size == 0:
            return True
        s = p.sum()
        if not np.isfinite(s) or s <= 0:
            return True
        p = p / s
        return float(np.max(np.abs(p - p.mean()))) < uniform_atol

    all_events = []
    if timing_engine.WGD:
        all_events.append(timing_engine.WGD)

    all_events.extend(_dedup_cn_events(timing_engine.all_cn_events))
    all_events.extend(
        mut for mut in timing_engine.mutations.values()
        if mut.prot_change
        and mut.gene in drivers
        and (mut.prot_change[0] != mut.prot_change[-1] or not mut.prot_change[-1].isalpha())
    )

    comps = {}

    for eve1, eve2 in itertools.combinations(all_events, 2):
        # Either event has no informative timing → fully unknown.
        if _is_uninformative(eve1.pi_dist) or _is_uninformative(eve2.pi_dist):
            comps[(eve1.event_name, eve2.event_name)] = (0.0, 0.0, 1.0)
            continue

        p1 = np.asarray(eve1.pi_dist, dtype=float)
        p2 = np.asarray(eve2.pi_dist, dtype=float)

        if p1.ndim != 1 or p2.ndim != 1 or len(p1) != len(p2):
            continue

        s1 = p1.sum()
        s2 = p2.sum()
        if s1 <= 0 or s2 <= 0:
            continue

        # normalize defensively
        p1 = p1 / s1
        p2 = p2 / s2

        # joint[i, j] = P(event1 in bin i AND event2 in bin j)
        joint = np.outer(p1, p2)

        # earlier = smaller pi bin
        p_before = np.triu(joint, k=uncertainty_bins + 1).sum()
        p_after = np.tril(joint, k=-(uncertainty_bins + 1)).sum()
        p_unknown = np.trace(joint) if uncertainty_bins == 0 else 0.0

        if uncertainty_bins > 0:
            for d in range(-uncertainty_bins, uncertainty_bins + 1):
                p_unknown += np.trace(joint, offset=d)

        # numerical cleanup
        total = p_before + p_after + p_unknown
        if total > 0:
            p_before /= total
            p_after /= total
            p_unknown /= total

        comps[(eve1.event_name, eve2.event_name)] = (p_before, p_after, p_unknown)

    return comps
