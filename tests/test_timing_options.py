"""Tests for the opt-in timing options, on synthetic segment tables and mutations run through the real readers.

Run from the repository root: python -m unittest tests.test_timing_options
"""
import gzip
import os
import random
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from data.Enums import Genome
from data.Patient import Patient
from LeagueModel.LeagueModelData import League
from output.PhylogicOutput import PhylogicOutput
from SinglePatientTiming import SinglePatientTiming as spt
from SinglePatientTiming.TimingEngine import TimingEngine, TimingOptions, Y_CALLABLE_BLOCKS_HG19

ALL_ON = dict(sex_aware=True, continuous_cn=True, clonal_unresolved=True, segment_homdels=True)
SEG_COLUMNS = ['Chromosome', 'Start.bp', 'End.bp', 'n_probes', 'seg_sigma', 'total_copy_ratio', 'rescaled_total_cn',
               'modal.a1', 'modal.a2', 'rescaled.cn.a1', 'rescaled.cn.a2', 'cancer.cell.frac.a1', 'ccf.ci95.low.a1',
               'ccf.ci95.high.a1', 'cancer.cell.frac.a2', 'ccf.ci95.low.a2', 'ccf.ci95.high.a2']
PURITY = 0.8
PLOIDY_FOR_COPY_RATIO = 2.0


def clonal_ccf_grid():
    grid = np.exp(-((np.linspace(0, 1, 101) - 1.) ** 2) / (2 * 0.04 ** 2))
    return grid / grid.sum()


def subclonal_ccf_grid(center=0.3):
    grid = np.exp(-((np.linspace(0, 1, 101) - center) ** 2) / (2 * 0.04 ** 2))
    return grid / grid.sum()


def seg_row(chrom, start, end, a1, a2, cont_a1=None, cont_a2=None, sigma=0.01, n_probes=500, total=None):
    cont_a1 = float(a1) if cont_a1 is None else cont_a1
    cont_a2 = float(a2) if cont_a2 is None else cont_a2
    total = cont_a1 + cont_a2 if total is None else total
    # Copy ratio on the line that links it to copy number for this purity, as in real segment tables.
    density = PURITY * PLOIDY_FOR_COPY_RATIO + 2 * (1 - PURITY)
    per_copy = (2. if str(chrom) in ('X', 'Y') and a1 == 0 and sigma < 0 else 1.)
    copy_ratio = (PURITY * total * per_copy + 2 * (1 - PURITY)) / density
    return [chrom, start, end, n_probes, sigma, copy_ratio, total, a1, a2, cont_a1, cont_a2,
            max(1 - np.nan_to_num(cont_a1, nan=1.), 0), 0, 1, max(1 - np.nan_to_num(cont_a2, nan=1.), 0), 0, 1]


def whole_genome(genome, male, baseline=(1, 1), skip=()):
    """Every autosome and X at the baseline state. For males X is 0/1 and Y is left for the test to add."""
    rows = []
    for chrom, size in genome.CHROM_DICT.items():
        if chrom == 'Y' or chrom in skip:
            continue
        a1, a2 = baseline
        if chrom == 'X' and male:
            a1, a2 = 0, baseline[1]
        rows.append(seg_row(chrom, 1000, size - 1000, a1, a2))
    return rows


def y_rows(cn, sigma=0.01):
    """chrY as rescued by map_to_absolute_copy_number: a1 is 0, and a2 holds the rounded value."""
    return [seg_row('Y', start, end, 0, int(round(cn)), cont_a1=0., cont_a2=cn, sigma=sigma, total=cn)
            for start, end in Y_CALLABLE_BLOCKS_HG19]


class TimingCase(unittest.TestCase):
    genome = Genome('hg19')

    def setUp(self):
        random.seed(1)
        np.random.seed(1)
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp)
        self.cwd = os.getcwd()

    def run_engine(self, seg_rows, muts=(), focal=None, options=None, artifact=None, sample='S1', purity=PURITY,
                   max_length=5000000):
        seg_path = os.path.join(self.tmp, sample + '.segtab.allelic.completed.txt.gz')
        with gzip.open(seg_path, 'wt') as handle:
            pd.DataFrame(seg_rows, columns=SEG_COLUMNS).to_csv(handle, sep='\t', index=False)
        maf_path = os.path.join(self.tmp, 'mut_ccfs.txt')
        grids = [f'preDP_ccf_{i / 100.0}' for i in range(101)]
        header = ['Patient_ID', 'Sample_ID', 'Hugo_Symbol', 'Chromosome', 'Start_position', 'Reference_Allele',
                  'Tumor_Seq_Allele', 't_ref_count', 't_alt_count', 'Protein_change', 'Variant_Classification',
                  'Variant_Type', 'Cluster_Assignment'] + grids
        # A clonal cluster must exist for clonality calls. A distant subclonal one makes CCF assignment meaningful.
        muts = list(muts) or [('1', 2000000 + i, 1, 'A') for i in range(5)] + [('2', 2000000 + i, 2, 'A') for i in range(5)]
        with open(maf_path, 'w') as handle:
            handle.write('\t'.join(header) + '\n')
            for index, (chrom, pos, cluster, gene) in enumerate(muts):
                grid = clonal_ccf_grid() if cluster == 1 else subclonal_ccf_grid()
                row = ['P1', sample, gene, chrom, pos, 'A', 'T', 40, 40 if cluster == 1 else 12, f'p.A{index}T',
                       'Missense_Mutation', 'SNP', cluster] + list(grid)
                handle.write('\t'.join(map(str, row)) + '\n')
        options = options or TimingOptions()
        patient = Patient(indiv_name='P1', genome_build='hg19')
        if options.male:
            patient.hemizygous_chroms = {'X', 'Y'}
        patient.addSample(maf_path, sample, input_type='post-clustering', timepoint_value=0, grid_size=101,
                          _additional_muts=None, seg_file=seg_path, purity=purity)
        patient.get_arm_level_cn_events()
        if focal is not None:
            focal_path = os.path.join(self.tmp, 'focal.txt')
            with open(focal_path, 'w') as handle:
                for chrom, start, end, kind in focal:
                    handle.write(f'{chrom}\t{start}\t{end}\t{kind}\t0.01\tregion\n')
            homdels = None
            if options.segment_homdels:
                homdels = dict(artifact_file=artifact, max_length=max_length, min_posterior=.95)
            patient.add_focal_cn_events(focal_regions=focal_path, segment_homdels=homdels)
        patient.preprocess_samples()
        engine = TimingEngine(patient, min_supporting_muts=3, options=options)
        engine.time_events()
        return engine

    def events(self, engine):
        return {name: events[0] for name, events in engine.all_cn_events.items() if events}


class MaleBaseline(TimingCase):
    def test_original_behavior_calls_normal_male_x_a_loss(self):
        rows = whole_genome(self.genome, male=True) + y_rows(1.)
        events = self.events(self.run_engine(rows))
        self.assertIn('loss_Xp', events)
        self.assertNotIn('loss_Y', events)

    def test_sex_aware_normal_male_has_no_sex_chromosome_events(self):
        rows = whole_genome(self.genome, male=True) + y_rows(1.)
        events = self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON)))
        self.assertEqual({name for name in events if 'X' in name or 'Y' in name}, set())

    def test_female_x_is_diploid_and_has_no_y_event(self):
        rows = whole_genome(self.genome, male=False)
        rows = [r for r in rows if r[0] != 'X'] + [seg_row('X', 1000, 60000000, 0, 1), seg_row('X', 62000000, 155000000, 0, 1)]
        events = self.events(self.run_engine(rows, options=TimingOptions(sex='XX', **ALL_ON)))
        self.assertIn('loss_Xp', events)  # a real loss in a female
        self.assertNotIn('loss_Y', events)

    def test_sex_aware_requires_sex(self):
        with self.assertRaises(ValueError):
            TimingOptions(sex_aware=True)


class LossOfY(TimingCase):
    def male(self, y_cn, baseline=(1, 1), **kwargs):
        rows = whole_genome(self.genome, male=True, baseline=baseline) + y_rows(y_cn, **kwargs)
        return self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON))

    def test_complete_loss_is_clonal_and_unresolved(self):
        event = self.events(self.male(0.))['loss_Y']
        self.assertEqual(event.Type, 'Chrom_loss')
        self.assertTrue(event.is_clonal)
        self.assertEqual(event.ccf_hat, 1.)
        self.assertEqual(event.timing_status, 'clonal_unresolved')
        self.assertIsNone(event.pi_dist)

    def test_partial_loss_is_subclonal_with_its_ccf(self):
        event = self.events(self.male(0.57))['loss_Y']
        self.assertFalse(event.is_clonal)
        self.assertAlmostEqual(event.ccf_hat, 0.43, places=6)
        self.assertEqual(event.timing_status, 'subclonal')
        self.assertEqual(int(np.argmax(event.pi_dist)), 100)

    def test_small_deviation_is_no_event(self):
        self.assertNotIn('loss_Y', self.events(self.male(0.9)))

    def test_too_little_callable_territory_gives_no_event(self):
        rows = whole_genome(self.genome, male=True) + y_rows(0.)[:1]
        self.assertNotIn('loss_Y', self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON))))

    def test_blocks_outside_callable_territory_are_ignored(self):
        rows = whole_genome(self.genome, male=True) + y_rows(1.)
        rows.append(seg_row('Y', 13130000, 13559000, 0, 0, cont_a1=0., cont_a2=0., total=0.))
        self.assertNotIn('loss_Y', self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON))))


class WgdRule(TimingCase):
    """A tetraploid male: enough 2/2 arms with timeable mutations for the WGD call."""

    def wgd_engine(self, y_cn=None, extra=()):
        rows = whole_genome(self.genome, male=True, baseline=(2, 2))
        rows += [] if y_cn is None else y_rows(y_cn)
        rows += list(extra)
        muts = []
        for chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']:
            muts += [(chrom, 2000000 + i * 1000, 1, 'A') for i in range(5)]
            muts += [(chrom, 140000000 + i * 1000, 1, 'A') for i in range(5)] if chrom in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12') else []
        return self.run_engine(rows, muts=muts, options=TimingOptions(sex='XY', **ALL_ON))

    def test_wgd_is_called_and_male_x_is_not_an_event(self):
        engine = self.wgd_engine(y_cn=2.)
        self.assertIsNotNone(engine.WGD)
        events = self.events(engine)
        self.assertEqual({name for name in events if 'X' in name or 'Y' in name}, set())

    def test_y_zero_precedes_wgd_and_y_one_follows_it(self):
        before = self.wgd_engine(y_cn=0.)
        after = self.wgd_engine(y_cn=1.)
        grid = np.linspace(0, 1, 101)
        mean_before = float((self.events(before)['loss_Y'].pi_dist * grid).sum())
        mean_after = float((self.events(after)['loss_Y'].pi_dist * grid).sum())
        wgd_before = float((before.WGD.pi_dist * grid).sum())
        wgd_after = float((after.WGD.pi_dist * grid).sum())
        self.assertLess(mean_before, wgd_before)
        self.assertGreater(mean_after, wgd_after)
        self.assertEqual(self.events(before)['loss_Y'].timing_status, 'estimated')
        self.assertEqual(self.events(before)['loss_Y'].allelic_cn, 0.)
        self.assertEqual(self.events(after)['loss_Y'].allelic_cn, 1.)


class ContinuousArmEvents(TimingCase):
    def arm_engine(self, cont_a1, modal_a1, continuous=True, cluster2=True):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, modal_a1, 1, cont_a1=cont_a1),
                 seg_row('17', 25000000, 81000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=continuous, clonal_unresolved=True)
        return self.run_engine(rows, options=options)

    def test_subclonal_arm_loss_rounded_to_normal_is_called(self):
        event = self.events(self.arm_engine(0.6, 1))['loss_17p']
        self.assertAlmostEqual(event.ccf_hat, 0.4)
        self.assertFalse(event.is_clonal)
        self.assertEqual(event.timing_status, 'subclonal')

    def test_integer_state_misses_it(self):
        self.assertNotIn('loss_17p', self.events(self.arm_engine(0.6, 1, continuous=False)))

    def test_rounded_down_loss_keeps_its_true_ccf(self):
        event = self.events(self.arm_engine(0.4, 0))['loss_17p']
        self.assertAlmostEqual(event.ccf_hat, 0.6)
        self.assertFalse(event.is_clonal)

    def test_clonal_loss_stays_clonal(self):
        event = self.events(self.arm_engine(0.03, 0))['loss_17p']
        self.assertAlmostEqual(event.ccf_hat, 0.97)
        self.assertTrue(event.is_clonal)
        self.assertEqual(event.timing_status, 'clonal_unresolved')

    def test_noise_near_normal_is_no_event(self):
        self.assertNotIn('loss_17p', self.events(self.arm_engine(0.9, 1)))

    def test_nan_continuous_value_falls_back_to_the_integer_state(self):
        event = self.events(self.arm_engine(float('nan'), 0))['loss_17p']
        self.assertEqual(event.ccf_hat, 1.)


class ClonalUnresolvedComparisons(TimingCase):
    @staticmethod
    def timed(name, pi=None, status=None):
        return spt.TimedEvent(name, pi, status)

    def test_rule(self):
        uniform = np.ones(101) / 101
        delta = np.zeros(101)
        delta[100] = 1.
        early = np.exp(-((np.linspace(0, 1, 101) - .2) ** 2) / .005)
        early /= early.sum()
        comps = spt.compare_event_list([
            self.timed('clonal_loss', None, 'clonal_unresolved'),
            self.timed('clonal_snv', uniform, 'clonal_unresolved'),
            self.timed('subclonal', delta, 'subclonal'),
            self.timed('timed', early, 'estimated'),
            self.timed('legacy', None, 'not_estimable'),
        ])
        self.assertEqual(comps[('clonal_loss', 'subclonal')], (1.0, 0.0, 0.0))
        self.assertEqual(comps[('clonal_snv', 'subclonal')], (1.0, 0.0, 0.0))
        self.assertEqual(comps[('clonal_loss', 'clonal_snv')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('clonal_loss', 'timed')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('clonal_loss', 'legacy')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('subclonal', 'legacy')], (0.0, 0.0, 1.0))
        self.assertGreater(comps[('subclonal', 'timed')][1], .99)

    def test_order_of_the_pair_does_not_matter(self):
        delta = np.zeros(101)
        delta[100] = 1.
        comps = spt.compare_event_list([self.timed('subclonal', delta, 'subclonal'),
                                        self.timed('clonal', None, 'clonal_unresolved')])
        self.assertEqual(comps[('subclonal', 'clonal')], (0.0, 1.0, 0.0))

    def test_without_statuses_nothing_changes(self):
        delta = np.zeros(101)
        delta[100] = 1.
        comps = spt.compare_event_list([self.timed('a', None), self.timed('b', delta)])
        self.assertEqual(comps[('a', 'b')], (0.0, 0.0, 1.0))

    def test_timing_table_round_trip(self):
        rows = whole_genome(self.genome, male=True) + y_rows(0.)
        rows = [r for r in rows if r[0] != '17'] + [seg_row('17', 1000, 22000000, 1, 1, cont_a1=.6),
                                                    seg_row('17', 25000000, 81000000, 1, 1)]
        engine = self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON))
        os.chdir(self.tmp)
        self.addCleanup(os.chdir, self.cwd)
        PhylogicOutput().write_timing_tsv(engine)
        table = pd.read_csv('P1.timing.tsv', sep='\t', low_memory=False)
        self.assertEqual(table.loc[table['Event Name'] == 'loss_Y', 'timing_status'].iloc[0], 'clonal_unresolved')
        comps = spt.compare_event_list(spt.read_timed_events('P1.timing.tsv', compare_all=True))
        pair = ('loss_Y', 'loss_17p') if ('loss_Y', 'loss_17p') in comps else ('loss_17p', 'loss_Y')
        expected = (1.0, 0.0, 0.0) if pair[0] == 'loss_Y' else (0.0, 1.0, 0.0)
        self.assertEqual(comps[pair], expected)

    def test_original_timing_table_keeps_its_columns(self):
        engine = self.run_engine(whole_genome(self.genome, male=False))
        os.chdir(self.tmp)
        self.addCleanup(os.chdir, self.cwd)
        PhylogicOutput().write_timing_tsv(engine)
        self.assertNotIn('timing_status', pd.read_csv('P1.timing.tsv', sep='\t', nrows=1).columns)


class SegmentHomdels(TimingCase):
    REGION = ('9', 21000000, 23000000, 'DEL')
    NAME = 'Focal_homdel_9:21000000-23000000'

    def rows(self, copies, sigma=0.01, start=21900000, end=22100000, a1=0, a2=0, n_probes=500, chrom='9', male=False):
        rows = [r for r in whole_genome(self.genome, male=male) if r[0] != chrom]
        size = self.genome.CHROM_DICT[chrom]
        flank = (0, 1) if male and chrom == 'X' else (1, 1)
        rows += [seg_row(chrom, 1000, start - 1, *flank),
                 seg_row(chrom, start, end, a1, a2, cont_a1=copies / 2., cont_a2=copies / 2., sigma=sigma,
                         n_probes=n_probes, total=copies),
                 seg_row(chrom, end + 1, size - 1000, *flank)]
        return rows

    def homdel(self, rows, focal=None, options=None, **kwargs):
        options = options or TimingOptions(sex='XX', **ALL_ON)
        return self.events(self.run_engine(rows, focal=focal or [self.REGION], options=options, **kwargs)).get(self.NAME)

    def test_clear_deletion_is_called_clonal_unresolved(self):
        event = self.homdel(self.rows(0.02))
        self.assertIsNotNone(event)
        self.assertEqual(event.timing_status, 'clonal_unresolved')
        self.assertTrue(event.is_clonal)

    def test_low_purity_rounding_is_rejected(self):
        # The original definition calls this modal 0/0 segment a homdel. Its copy number is closer to 1 than to 0.
        self.assertIsNone(self.homdel(self.rows(0.59, sigma=0.02)))

    def test_original_definition_accepts_it(self):
        options = TimingOptions(sex='XX', sex_aware=True, clonal_unresolved=True)
        self.assertIsNotNone(self.homdel(self.rows(0.59, sigma=0.02, start=21000000, end=23000000), options=options))

    def test_noisy_segment_is_rejected_by_the_posterior(self):
        self.assertIsNone(self.homdel(self.rows(0.3, sigma=0.5)))

    def test_too_long_segment_is_rejected(self):
        self.assertIsNone(self.homdel(self.rows(0.02, start=20000000, end=26000000)))

    def test_too_few_probes_is_rejected(self):
        self.assertIsNone(self.homdel(self.rows(0.02, n_probes=2)))

    def test_artifact_region_is_rejected(self):
        path = os.path.join(self.tmp, 'artifacts.txt')
        with open(path, 'w') as handle:
            handle.write('#chromosome\tstart\tend\trec\n9\t21800000\t22200000\t5\n')
        self.assertIsNone(self.homdel(self.rows(0.02), artifact=path))

    def test_amplification_region_is_ignored(self):
        self.assertIsNone(self.homdel(self.rows(0.02), focal=[('9', 21000000, 23000000, 'AMP')]))

    def test_modal_state_does_not_decide(self):
        self.assertIsNotNone(self.homdel(self.rows(0.02, a1=0, a2=1)))

    def test_subclonal_deletion_has_ccf_relative_to_flanks(self):
        engine = self.run_engine(self.rows(0.4), focal=[self.REGION], options=TimingOptions(sex='XX', **ALL_ON))
        event = self.events(engine)[self.NAME]
        cnv = [c for c in engine.sample_list[0].sample.low_coverage_mutations.values()
               if getattr(c, 'cn_category', '') == 'Focal_homdel'][0]
        self.assertAlmostEqual(cnv.ccf_hat, 0.8)
        self.assertFalse(event.is_clonal is None)

    def test_male_x_counts_the_single_copy(self):
        region = ('X', 76000000, 78000000, 'DEL')
        name = 'Focal_homdel_X:76000000-78000000'
        rows = self.rows(0.02, start=76900000, end=77100000, chrom='X', male=True)
        events = self.events(self.run_engine(rows, focal=[region], options=TimingOptions(sex='XY', **ALL_ON)))
        self.assertIn(name, events)

    def test_male_x_normal_copy_is_not_a_deletion(self):
        region = ('X', 76000000, 78000000, 'DEL')
        rows = [r for r in whole_genome(self.genome, male=True)]
        events = self.events(self.run_engine(rows, focal=[region], options=TimingOptions(sex='XY', **ALL_ON)))
        self.assertEqual([name for name in events if 'X' in name], [])

    def test_wgd_tumor_times_homdel_before_wgd(self):
        rows = [r for r in whole_genome(self.genome, male=False, baseline=(2, 2)) if r[0] != '9']
        size = self.genome.CHROM_DICT['9']
        rows += [seg_row('9', 1000, 21899999, 2, 2), seg_row('9', 21900000, 22100000, 0, 0, cont_a1=.01, cont_a2=.01, total=.02),
                 seg_row('9', 22100001, size - 1000, 2, 2)]
        muts = [(c, 2000000 + i * 1000, 1, 'A') for c in map(str, range(1, 13)) for i in range(5)]
        muts += [(c, 140000000 + i * 1000, 1, 'A') for c in map(str, range(1, 8)) for i in range(5)]
        engine = self.run_engine(rows, muts=muts, focal=[self.REGION], options=TimingOptions(sex='XX', **ALL_ON))
        self.assertIsNotNone(engine.WGD)
        event = self.events(engine)[self.NAME]
        grid = np.linspace(0, 1, 101)
        self.assertEqual(event.timing_status, 'estimated')
        self.assertLess(float((event.pi_dist * grid).sum()), float((engine.WGD.pi_dist * grid).sum()))


class LeagueInformation(unittest.TestCase):
    def league(self):
        league = object.__new__(League)
        league.final_event_list = ['a', 'b', 'c']
        league.event_prev = {'a': 5, 'b': 5, 'c': 1}
        league.mut_type = {}
        league.event_pairs_per_samp = {
            's1': {('a', 'b'): {('a', 'b'): .9, ('b', 'a'): .1, 'unknown': 0.}, ('a', 'c'): {('a', 'c'): 0., ('c', 'a'): 0., 'unknown': 1.}},
            's2': {('a', 'b'): {('a', 'b'): .8, ('b', 'a'): .1, 'unknown': .1}},
        }
        league.form_pairs_for_league_model()
        league.update_pairs_for_league_model()
        return league

    def test_counts_exclude_the_low_count_rescue(self):
        rows = self.league().get_event_information()
        self.assertAlmostEqual(rows['a']['n_informative_games'], 1.9)
        self.assertAlmostEqual(rows['a']['n_unknown_games'], 1.1)
        self.assertAlmostEqual(rows['c']['n_informative_games'], 0.)
        self.assertAlmostEqual(rows['c']['draw_fraction'], 1.)  # b and c never met, a and c only as unknown
        self.assertAlmostEqual(rows['b']['n_unknown_games'], .1)

    def test_scores_keep_the_rescue(self):
        league = self.league()
        pair = league.event_pairs[('b', 'c')]
        self.assertEqual(pair.observed_rates['b'] + pair.observed_rates['c'], 0)
        self.assertEqual(pair.win_rates['b'] + pair.win_rates['c'], 2)


if __name__ == '__main__':
    unittest.main()
