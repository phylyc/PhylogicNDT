import math
import os
import pickle
import tempfile
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import LeagueModel.LeagueModelData as LeagueModelData

# Module-level variable set by _init_worker in each worker process.
_worker_league = None


def _init_worker(league_pickle_path):
    """Initializer for each worker process: load the pre-built League from a
    temp pickle file so that expensive CSV parsing / League construction
    happens only once per worker, not once per chunk."""
    global _worker_league
    with open(league_pickle_path, 'rb') as fh:
        _worker_league = pickle.load(fh)


def _run_perm_chunk(payload):
    """Run a chunk of permutations using the pre-initialised League object."""
    import random
    import numpy as np

    seed = int(payload['seed'])
    random.seed(seed)
    np.random.seed(seed)

    league = _worker_league  # already loaded by _init_worker
    all_samps = sorted(set(league.events_per_samp.keys()))
    n_subset = int(len(all_samps) * float(payload['percent_subset']))
    n_subset = max(1, min(len(all_samps), n_subset)) if all_samps else 0
    out = []
    for _ in range(int(payload['n_perms_chunk'])):
        rand_subset = set(random.sample(all_samps, n_subset)) if n_subset > 0 else set()
        league.run_permutation(
            num_seasons=int(payload['n_seasons']),
            samples=rand_subset,
            final_event_list=payload['final_event_list'],
        )
        out.append(league.calc_odds())
    return out


def _append_odds_dict(league_model_run, odds_dict):
    for event, event_odds in odds_dict.items():
        if event not in league_model_run.odds:
            league_model_run.odds[event] = {'odds_early': [], 'odds_late': []}
        league_model_run.odds[event]['odds_early'].append(event_odds['odds_early'])
        league_model_run.odds[event]['odds_late'].append(event_odds['odds_late'])
    league_model_run.n_perms += 1


def _run_permutations_serial(league_model_run, all_samps, args):
    import random
    import logging

    n_subset = int(len(all_samps) * args.percent_subset)
    n_subset = max(1, min(len(all_samps), n_subset)) if all_samps else 0
    for j in range(args.n_perms):
        logging.info('running iter:' + str(j))
        rand_subset = set(random.sample(all_samps, n_subset)) if n_subset > 0 else set()
        league_model_run.run_permutation(num_seasons=args.n_seasons, samples=rand_subset,
                                         final_event_list=league_model_run.final_event_list)
        league_model_run.update_odds()


def _run_permutations_parallel(league_model_run, args):
    import logging

    n_jobs = max(1, args.n_jobs)
    chunk_size = max(1, args.perm_chunk_size)
    n_chunks = int(math.ceil(float(args.n_perms) / float(chunk_size)))
    final_event_list = list(league_model_run.final_event_list)
    base_seed = int(getattr(args, 'random_seed', getattr(args, 'seed', None) or 1))

    # Pickle the pre-built League to a temp file.  Workers load this once
    # via the initializer instead of re-reading the CSV and re-running the
    # expensive League.__init__ for every chunk.
    saved_df = league_model_run.query_res_df      # large, not needed by workers
    league_model_run.query_res_df = None
    league_pickle_fd, league_pickle_path = tempfile.mkstemp(suffix='.league.pkl')
    try:
        with os.fdopen(league_pickle_fd, 'wb') as fh:
            pickle.dump(league_model_run, fh)
    finally:
        league_model_run.query_res_df = saved_df  # restore for the main process

    payloads = []
    for chunk_id in range(n_chunks):
        n_chunk = min(chunk_size, int(args.n_perms) - chunk_id * chunk_size)
        payloads.append({
            'n_perms_chunk': int(n_chunk),
            'seed': base_seed + chunk_id,
            'n_seasons': int(args.n_seasons),
            'percent_subset': float(args.percent_subset),
            'final_event_list': final_event_list,
        })

    logging.info('running %d permutation(s) across %d worker(s) in %d chunk(s)', int(args.n_perms), n_jobs, n_chunks)
    try:
        mp_ctx = multiprocessing.get_context('fork')
        with ProcessPoolExecutor(max_workers=n_jobs,
                                 mp_context=mp_ctx,
                                 initializer=_init_worker,
                                 initargs=(league_pickle_path,)) as ex:
            futures = {ex.submit(_run_perm_chunk, payload): idx for idx, payload in enumerate(payloads)}
            for fut in as_completed(futures):
                chunk_idx = futures[fut]
                odds_list = fut.result()
                for odds_dict in odds_list:
                    _append_odds_dict(league_model_run, odds_dict)
                logging.info('completed permutation chunk %d/%d', chunk_idx + 1, n_chunks)
    finally:
        try:
            os.unlink(league_pickle_path)
        except OSError:
            pass


def run_league_model(args):
    print(args)
    import pandas as pd
    import pickle as pkl
    import matplotlib.pyplot as plt
    import logging

    if args.comparison_fn is not None:
        full_comp_fn = args.comparison_fn
    elif args.comparison_fn is None and args.comps is not None:
        full_comp_fn = args.cohort + '.comparisons.tsv.gz'
        pd.concat([
            pd.read_csv(fn, sep="\t", header=0, names=['sample', 'event1', 'event2', 'p_event1_win', 'p_event2_win', 'unknown'], low_memory=False)
            for fn in args.comps
        ]).to_csv(full_comp_fn, sep="\t", index=False)
    else:
        logging.error('Please Provide Input Data for League Model')
        return 0

    query_res_df1 = pd.read_csv(full_comp_fn, sep='\t')

    final_samples = None
    if args.final_sample_list is not None:
        final_samples = []
        with open(args.final_sample_list, 'r') as final_sample_list_fn:
            for i, row in enumerate(final_sample_list_fn):
                final_samples.append(row.strip("\n"))

    final_events = None
    if args.final_event_list is not None:
        final_events = []
        with open(args.final_event_list, 'r') as final_event_list_fn:
            for i, row in enumerate(final_event_list_fn):
                final_events.append(row.strip("\n"))

    league_model_run = LeagueModelData.League(
        query_res_df1,
        cohort=args.cohort,
        keep_only_samples_w_event=args.keep_samps_w_event,
        remove_samps_w_event=args.remove_samps_w_event,
        num_games_against_each_opponent=args.num_games_against_each_opponent,
        final_samples=final_samples,
        final_event_list=final_events,
        max_num_snvs=args.max_num_snvs,
        max_num_cnv_arms=args.max_num_cnv_arms,
        max_num_cnv_arm_gains=args.max_num_cnv_arm_gains,
        max_num_cnv_arm_losses=args.max_num_cnv_arm_losses,
        max_num_cnv_focal=args.max_num_cnv_focal,
        max_num_homdel=args.max_num_homdel,
        min_event_prevalence=args.min_event_prevalence,
        mutation_event_resolution=args.mutation_event_resolution,
        hotspot_fn=args.hotspot_fn,
        drop_silent_mutations=args.drop_silent_mutations
    )
    output = open(args.cohort + '.all_events.tsv', 'w')
    output.write('indiv\tevent\n')
    for indiv in league_model_run.events_per_samp_full:
        for event in league_model_run.events_per_samp_full[indiv]:
            output.write(indiv + '\t' + event + '\n')
    output.close()

    if getattr(league_model_run, 'event_label_map', None):
        pd.DataFrame(
            sorted(league_model_run.event_label_map.values(), key=lambda x: (str(x.get('harmonized_event', '')), str(x.get('raw_event', ''))))
        ).to_csv(args.cohort + '.event_label_map.tsv', sep='\t', index=False)

    all_samps = set(league_model_run.events_per_samp.keys())
    league_model_run.run_full_run(num_seasons=0, samples=all_samps, final_event_list=league_model_run.final_event_list)
    league_model_run.init_odds(league_model_run.final_event_list)

    ## output ##
    out_fn = args.cohort + '.matrix_of_comparisons.tsv'
    sorted_final_event_list = sorted(league_model_run.final_event_list)
    rev_list = sorted_final_event_list[::-1]
    with open(out_fn, 'w') as output:
        output.write('event_v_event\t' + '\t'.join(sorted_final_event_list) + '\n')
        for i, event1 in enumerate(rev_list):
            output.write(event1 + '\t')
            for j, event2 in enumerate(sorted_final_event_list):
                if event1 == event2:
                    output.write("na")
                else:
                    event_pair = tuple(sorted([event2, event1]))
                    if j <= len(sorted_final_event_list) - i - 1:
                        output.write(','.join(map(str, [league_model_run.event_pairs[event_pair].win_rates[event1],
                                                        league_model_run.event_pairs[event_pair].win_rates[event2],
                                                        league_model_run.event_pairs[event_pair].win_rates[
                                                            'unknown']])))
                    else:
                        output.write(','.join(map(str, [league_model_run.event_pairs[event_pair].win_rates[event2],
                                                        league_model_run.event_pairs[event_pair].win_rates[event1],
                                                        league_model_run.event_pairs[event_pair].win_rates[
                                                            'unknown']])))
                if j < len(sorted_final_event_list) - 1:
                    output.write('\t')
            output.write('\n')

    sorted_all_samps = sorted(all_samps)
    n_jobs = int(max(1, args.n_jobs))
    if int(args.n_perms) > 0:
        if n_jobs > 1 and int(args.n_perms) > 1:
            _run_permutations_parallel(league_model_run, args)
        else:
            _run_permutations_serial(league_model_run, sorted_all_samps, args)

    # Run a final set of seasons on the full (un-subsetted) data to
    # populate event_positions for the position plot and set num_seasons
    # for the odds plot x-axis.  In parallel mode workers modify their own
    # copies, leaving the main process stale from run_full_run(num_seasons=0).
    league_model_run.run_permutation(
        num_seasons=int(args.n_seasons),
        samples=all_samps,
        final_event_list=league_model_run.final_event_list,
    )

    league_model_run.calc_log_odds_full_run()

    ## plotting ##
    league_model_run.plot_league_run(type='odds')
    fig = league_model_run.odds_plot
    if fig is not None:
        plt.savefig(args.cohort + '.log_odds.pdf', transparent=True, bbox_inches="tight")
    # plt.savefig(args.cohort + '.log_odds.png', dpi=300, bbox_inches="tight")
    league_model_run.plot_league_run(type='pos')
    fig = league_model_run.pos_plot
    if fig is not None:
        plt.savefig(args.cohort + '.positions.pdf', transparent=True, bbox_inches="tight")
    # plt.savefig(args.cohort + '.positions.png', dpi=300, bbox_inches="tight")

    ## output ##
    out_fn = args.cohort + '.prevalence.tsv'
    with open(out_fn, 'w') as output:
        output.write('cohort\tevent\tn_occur\tevent_split\ttype\tn_samp\n')
        for eve in league_model_run.final_event_list:
            if args.keep_samps_w_event is not None:
                output.write(args.cohort + '\t' + eve + '\t' + str(
                    league_model_run.full_event_prev[eve]) + '\t' + args.keep_samps_w_event +
                             '\twith\t' + str(league_model_run.full_n_samps) + '\n')
            elif args.remove_samps_w_event is not None:
                output.write(args.cohort + '\t' + eve + '\t' + str(
                    league_model_run.full_event_prev[eve]) + '\t' + args.remove_samps_w_event +
                             '\twithout\t' + str(league_model_run.full_n_samps) + '\n')
            else:
                output.write(
                    args.cohort + '\t' + eve + '\t' + str(league_model_run.full_event_prev[eve]) + '\tNone\tna\t' +
                    str(league_model_run.full_n_samps) + '\n')

    out_fn = args.cohort + '.full_prevalence.tsv'
    with open(out_fn, 'w') as output:
        output.write('cohort\tevent\tn_occur\tn_samp\n')
        for eve in league_model_run.full_event_prev:
            output.write(args.cohort + '\t' + eve + '\t' + str(league_model_run.full_event_prev[eve]) + '\t' + str(
                league_model_run.full_n_samps) + '\n')

    out_fn = args.cohort + '.log_odds.tsv'
    with open(out_fn, 'w') as output:
        output.write('cohort\tevent\tevent_split\ttype\tperm_run\tlog_odds_early\n')
        for eve in league_model_run.log_odds_full_run.keys():
            for j, odds in enumerate(league_model_run.log_odds_full_run[eve]):
                if args.keep_samps_w_event is not None:
                    output.write(args.cohort + '\t' + eve + '\t' + str(j) + '\t' + str(
                        league_model_run.log_odds_full_run[eve][j]) +
                                 '\t' + args.keep_samps_w_event + '\twith\t' + '\n')
                elif args.remove_samps_w_event is not None:
                    output.write(args.cohort + '\t' + eve + '\t' + str(j) + '\t' + str(
                        league_model_run.log_odds_full_run[eve][j]) +
                                 '\t' + args.remove_samps_w_event + '\twithout\t' + '\n')
                else:
                    output.write(args.cohort + '\t' + eve + '\t' + str(j) + '\t' + str(
                        league_model_run.log_odds_full_run[eve][j]) +
                                 '\tNone\tna\t' + '\n')

    pkl.dump(league_model_run, open(args.cohort + 'fig_pkl.pkl', 'wb'))


def run_league_model_ipython_notebook(args):
    import pandas as pd
    import pickle as pkl
    import random
    import matplotlib.pyplot as plt
    import logging

    if args.comparison_fn is not None:
        full_comp_fn = args.comparison_fn
    elif args.comparison_fn is None and args.comps is not None:
        full_comp_fn = args.cohort + '.comparisons.tsv'
        output = open(full_comp_fn, 'w')
        output.write('\t'.join(['sample', 'event1', 'event2', 'p_event1_win', 'p_event2_win', 'unknown']) + '\n')
        for fn in args.comps:
            with open(fn, 'r') as input_fn:
                for i, row in enumerate(input_fn):
                    if i == 0: continue
                    output.write(row)
        output.close()
    else:
        logging.error('Please Provide Input Data for League Model')
        return 0

    query_res_df1 = pd.read_csv(full_comp_fn, sep='\t')
    final_samples = None
    if args.final_sample_list is not None:
        final_samples = []
        with open(args.final_sample_list, 'r') as final_sample_list_fn:
            for i, row in enumerate(final_sample_list_fn):
                final_samples.append(row.strip("\n"))
    league_model_run = LeagueModelData.League(query_res_df1, cohort=args.cohort,
                                                       keep_only_samples_w_event=args.keep_samps_w_event,
                                                       remove_samps_w_event=args.remove_samps_w_event,
                                                       num_games_against_each_opponent=args.num_games_against_each_opponent,
                                                       final_samples=final_samples,
                                                       mutation_event_resolution=args.mutation_event_resolution,
                                                       hotspot_fn=getattr(args, 'hotspot_fn', None),
                                                       drop_silent_mutations=getattr(args, 'drop_silent_mutations', None))
    output = open(args.cohort + '.all_events.tsv', 'w')
    output.write('indiv\tevent\n')
    for indiv in league_model_run.events_per_samp_full:
        for event in league_model_run.events_per_samp_full[indiv]:
            output.write(indiv + '\t' + event + '\n')
    output.close()

    all_samps = set(league_model_run.events_per_samp.keys())
    league_model_run.run_full_run(num_seasons=0, samples=all_samps, final_event_list=league_model_run.final_event_list)
    league_model_run.init_odds(league_model_run.final_event_list)

    ## output ##
    out_fn = args.cohort + '.matrix_of_comparisons.tsv'
    sorted_final_event_list = sorted(league_model_run.final_event_list)
    rev_list = sorted_final_event_list[::-1]
    with open(out_fn, 'w') as output:
        output.write('event_v_event\t' + '\t'.join(sorted_final_event_list) + '\n')
        for i, event1 in enumerate(rev_list):
            output.write(event1 + '\t')
            for j, event2 in enumerate(sorted_final_event_list):
                if event1 == event2:
                    output.write("na")
                else:
                    event_pair = tuple(sorted([event2, event1]))
                    if j <= len(sorted_final_event_list) - i - 1:
                        output.write(','.join(map(str, [league_model_run.event_pairs[event_pair].win_rates[event1],
                                                        league_model_run.event_pairs[event_pair].win_rates[event2],
                                                        league_model_run.event_pairs[event_pair].win_rates[
                                                            'unknown']])))
                    else:
                        output.write(','.join(map(str, [league_model_run.event_pairs[event_pair].win_rates[event2],
                                                        league_model_run.event_pairs[event_pair].win_rates[event1],
                                                        league_model_run.event_pairs[event_pair].win_rates[
                                                            'unknown']])))
                if j < len(sorted_final_event_list) - 1:
                    output.write('\t')
            output.write('\n')

    sorted_all_samps = sorted(all_samps)
    n_jobs = int(max(1, args.n_jobs))
    if int(args.n_perms) > 0:
        if n_jobs > 1 and int(args.n_perms) > 1:
            _run_permutations_parallel(league_model_run, args)
        else:
            _run_permutations_serial(league_model_run, sorted_all_samps, args)

    # Run a final set of seasons on the full (un-subsetted) data to
    # populate event_positions for the position plot and set num_seasons
    # for the odds plot x-axis.
    league_model_run.run_permutation(
        num_seasons=int(args.n_seasons),
        samples=all_samps,
        final_event_list=league_model_run.final_event_list,
    )

    league_model_run.calc_log_odds_full_run()

    ## plotting ##
    league_model_run.plot_league_run(type='odds')
    odds_plot = league_model_run.odds_plot
    plt.savefig(args.cohort + '.log_odds.pdf', transparent=True)
    # plt.savefig(args.cohort + '.log_odds.png')
    league_model_run.plot_league_run(type='pos')
    pos_plot = league_model_run.pos_plot
    plt.savefig(args.cohort + '.positions.pdf', transparent=True)
    # plt.savefig(args.cohort + '.positions.png')

    ## output ##
    out_fn = args.cohort + '.prevalence.tsv'
    with open(out_fn, 'w') as output:
        output.write('cohort\tevent\tn_occur\tevent_split\ttype\tn_samp\n')
        for eve in league_model_run.final_event_list:
            if args.keep_samps_w_event is not None:
                output.write(args.cohort + '\t' + eve + '\t' + str(
                    league_model_run.full_event_prev[eve]) + '\t' + args.keep_samps_w_event +
                             '\twith\t' + str(league_model_run.full_n_samps) + '\n')
            elif args.remove_samps_w_event is not None:
                output.write(args.cohort + '\t' + eve + '\t' + str(
                    league_model_run.full_event_prev[eve]) + '\t' + args.remove_samps_w_event +
                             '\twithout\t' + str(league_model_run.full_n_samps) + '\n')
            else:
                output.write(
                    args.cohort + '\t' + eve + '\t' + str(league_model_run.full_event_prev[eve]) + '\tNone\tna\t' +
                    str(league_model_run.full_n_samps) + '\n')

    out_fn = args.cohort + '.full_prevalence.tsv'
    with open(out_fn, 'w') as output:
        output.write('cohort\tevent\tn_occur\tn_samp\n')
        for eve in league_model_run.full_event_prev:
            output.write(args.cohort + '\t' + eve + '\t' + str(league_model_run.full_event_prev[eve]) + '\t' + str(
                league_model_run.full_n_samps) + '\n')

    out_fn = args.cohort + '.log_odds.tsv'
    with open(out_fn, 'w') as output:
        output.write('cohort\tevent\tevent_split\ttype\tperm_run\tlog_odds_early\n')
        for eve in league_model_run.log_odds_full_run.keys():
            for j, odds in enumerate(league_model_run.log_odds_full_run[eve]):
                if args.keep_samps_w_event is not None:
                    output.write(args.cohort + '\t' + eve + '\t' + str(j) + '\t' + str(
                        league_model_run.log_odds_full_run[eve][j]) +
                                 '\t' + args.keep_samps_w_event + '\twith\t' + '\n')
                elif args.remove_samps_w_event is not None:
                    output.write(args.cohort + '\t' + eve + '\t' + str(j) + '\t' + str(
                        league_model_run.log_odds_full_run[eve][j]) +
                                 '\t' + args.remove_samps_w_event + '\twithout\t' + '\n')
                else:
                    output.write(args.cohort + '\t' + eve + '\t' + str(j) + '\t' + str(
                        league_model_run.log_odds_full_run[eve][j]) +
                                 '\tNone\tna\t' + '\n')

    pkl.dump(league_model_run, open(args.cohort + '.fig_pkl.pkl', 'wb'))
    return odds_plot, pos_plot