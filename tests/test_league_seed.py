"""--seed makes LeagueModel reproducible: the main run, the permutations (serial and in workers) and the adaptive
refinement all draw from Generators spawned from the seed. Without --seed the global state is used, as before.

Run from the repository root: PYTHONPATH=. python -m unittest tests.test_league_seed
"""
import glob
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from LeagueModel.LeagueModelData import League

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVENTS = ['WGD', 'gain_1q', 'gain_7', 'gain_17q', 'loss_1p', 'loss_3p', 'loss_11q', 'loss_4p']


def comparisons(n_samples=12, seed=0):
    """A comparison table: every sample carries a random subset of EVENTS, every pair a random win probability."""
    state = np.random.RandomState(seed)
    rows = []
    for s in range(n_samples):
        events = [e for e in EVENTS if state.rand() < .7]
        for e1, e2 in itertools.combinations(events, 2):
            p1, p2 = state.dirichlet([1., 1., .5])[:2]
            rows.append({'sample': 'S%02d' % s, 'event1': e1, 'event2': e2,
                         'p_event1_win': p1, 'p_event2_win': p2, 'unknown': 1. - p1 - p2})
    return pd.DataFrame(rows)


class LeagueSeed(unittest.TestCase):
    def positions(self, seed, num_seasons=50):
        league = League(comparisons(), cohort=os.path.join(self.tmp, 'c'), final_event_list=list(EVENTS), seed=seed)
        league.run_permutation(num_seasons=num_seasons, samples=set(league.events_per_samp), final_event_list=list(EVENTS))
        return league.event_positions

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_same_seed_same_positions_whatever_the_global_state(self):
        np.random.seed(1)
        first = self.positions(7)
        np.random.seed(2)
        self.assertEqual(self.positions(7), first)

    def test_different_seeds_draw_differently(self):
        self.assertNotEqual(self.positions(7), self.positions(8))

    def test_seed_leaves_the_global_state_alone(self):
        np.random.seed(5)
        expected = np.random.rand()
        np.random.seed(5)
        self.positions(7)
        self.assertEqual(np.random.rand(), expected)

    def test_without_a_seed_the_global_state_is_used(self):
        np.random.seed(3)
        first = self.positions(None)
        np.random.seed(3)
        self.assertEqual(self.positions(None), first)
        np.random.seed(4)
        self.assertNotEqual(self.positions(None), first)


class LeagueModelCliSeed(unittest.TestCase):
    """The command line as Step 14d runs it (adaptive refinement), plus permutations in workers and serially."""
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.comp_fn = os.path.join(cls.tmp, 'comparisons.tsv')
        comparisons().to_csv(cls.comp_fn, sep='\t', index=False)
        cls.event_fn = os.path.join(cls.tmp, 'events.txt')
        with open(cls.event_fn, 'w') as fh:
            fh.write('\n'.join(EVENTS) + '\n')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def run_cli(self, name, *extra, hash_seed='0'):
        out = os.path.join(self.tmp, name)
        os.makedirs(out)
        cmd = [sys.executable, os.path.join(ROOT, 'PhylogicNDT.py'), 'LeagueModel', '--cohort', 'c',
               '--comparison_fn', self.comp_fn, '--final_event_list', self.event_fn,
               '--n_seasons', '64', '--num_games_against_each_opponent', '2',
               '--adaptive_refinement', '--adaptive_local_n_seasons', '64', '--adaptive_max_depth', '2',
               '--adaptive_support_threshold', '0.2', '--adaptive_min_window_events', '2'] + list(extra)
        env = dict(os.environ, PYTHONPATH=ROOT, PYTHONHASHSEED=hash_seed)
        result = subprocess.run(cmd, cwd=out, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, result.stdout.decode()[-2000:])
        tables = {}
        for fn in sorted(glob.glob(os.path.join(out, '*.tsv'))):
            with open(fn, 'rb') as fh:
                tables[os.path.basename(fn)] = fh.read()
        # The one table whose order is not the seed's: each sample's events come from a set, in PYTHONHASHSEED order.
        tables['c.all_events.tsv'] = sorted(tables['c.all_events.tsv'].splitlines())
        self.assertIn('c.adaptive_position_posterior.tsv', tables)
        return tables

    def test_same_seed_identical_tables_whatever_the_hash_seed(self):
        for name, extra in (('main', ['--n_perms', '0']),
                            ('workers', ['--n_perms', '6', '--n_jobs', '2', '--perm_chunk_size', '2']),
                            ('serial', ['--n_perms', '3', '--n_jobs', '1'])):
            with self.subTest(name):
                first = self.run_cli(name + '1', '--seed', '11', *extra)
                second = self.run_cli(name + '2', '--seed', '11', *extra, hash_seed='1')
                self.assertEqual(first, second)
                third = self.run_cli(name + '3', '--seed', '12', *extra)
                self.assertNotEqual(first['c.adaptive_position_posterior.tsv'],
                                    third['c.adaptive_position_posterior.tsv'])
                if name != 'main':
                    self.assertNotEqual(first['c.log_odds.tsv'], third['c.log_odds.tsv'])


if __name__ == '__main__':
    unittest.main()
