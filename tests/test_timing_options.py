"""Tests for the opt-in timing options, on synthetic segment tables and mutations run through the real readers.

Segment rows follow real ABSOLUTE tables: rescaled.cn sits on the modal integer and the signal is in hscr, except on male
X and Y, whose rows are rescued as (0, copy number) with no hscr.

Run from the repository root: python -m unittest tests.test_timing_options
"""
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from data.Enums import Genome
from data.Patient import Patient
from LeagueModel.LeagueModelData import League
from output.PhylogicOutput import PhylogicOutput
from SinglePatientTiming import SinglePatientTiming as spt
from SinglePatientTiming.TimingEngine import TimingEngine, TimingMut, TimingOptions, TERMINAL, Y_CALLABLE_BLOCKS_HG19

ALL_ON = dict(sex_aware=True, continuous_cn=True, clonal_unresolved=True, segment_homdels=True)
SEG_COLUMNS = ['Chromosome', 'Start.bp', 'End.bp', 'n_probes', 'seg_sigma', 'total_copy_ratio', 'rescaled_total_cn',
               'modal.a1', 'modal.a2', 'hscr.a1', 'hscr.a2', 'rescaled.cn.a1', 'rescaled.cn.a2', 'cancer.cell.frac.a1',
               'ccf.ci95.low.a1', 'ccf.ci95.high.a1', 'cancer.cell.frac.a2', 'ccf.ci95.low.a2', 'ccf.ci95.high.a2']
PURITY = 0.8
PLOIDY = 2.0
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def comb(purity=PURITY, ploidy=PLOIDY):
    density = purity * ploidy + 2 * (1 - purity)
    return 2 * (1 - purity) / density, purity / density


def seg_row(chrom, start, end, a1, a2, q1=None, q2=None, sigma=0.01, n_probes=500, purity=PURITY):
    """Autosomal or female X row. a1 and a2 are the modal state, q1 and q2 the true allelic copy numbers."""
    b, delta = comb(purity)
    q1 = float(a1) if q1 is None else q1
    q2 = float(a2) if q2 is None else q2
    hscr1 = '' if np.isnan(q1) else b + 2 * delta * q1
    hscr2 = '' if np.isnan(q2) else b + 2 * delta * q2
    known1, known2 = np.nan_to_num(q1, nan=float(a1)), np.nan_to_num(q2, nan=float(a2))
    total = known1 + known2
    # As in real tables, cancer.cell.frac is 1 - q for a loss, and for a gain only the distance above an integer.
    def fraction(q):
        return 1. - q if q < 1 else q - np.floor(q)
    ccf1, ccf2 = fraction(known1), fraction(known2)
    return [chrom, start, end, n_probes, sigma, b + delta * total, total, a1, a2, hscr1, hscr2, float(a1), float(a2),
            ccf1, 0, 1, ccf2, 0, 1]


def hemizygous_row(chrom, start, end, copies, sigma=0.01, n_probes=500, purity=PURITY):
    """Male X or Y row, as in real tables: one copy moves the copy ratio by 2 * delta, and a1 sits on the floor b."""
    b, delta = comb(purity)
    ratio = b + 2 * delta * copies
    return [chrom, start, end, n_probes, sigma, ratio, copies, 0, int(round(copies)), b, ratio, 0., copies,
            0, 0, 1, abs(copies - 1), 0, 1]


def anchor_rows(purity=PURITY):
    """Small segments of another state: real tables always vary, which is what reveals the purity and ploidy line."""
    return [seg_row('21', 15000000, 15020000, 1, 2, purity=purity)]


def whole_genome(genome, male, baseline=(1, 1), skip=(), purity=PURITY):
    """Every autosome and X at the baseline state. For males X is hemizygous and Y is left for the test to add."""
    rows = list(anchor_rows(purity))
    for chrom, size in genome.CHROM_DICT.items():
        if chrom == 'Y' or chrom in skip:
            continue
        if chrom == 'X' and male:
            rows.append(hemizygous_row('X', 1000, size - 1000, float(baseline[1]), purity=purity))
        elif chrom == '21':
            rows.append(seg_row(chrom, 15030000, size - 1000, *baseline, purity=purity))
        else:
            rows.append(seg_row(chrom, 1000, size - 1000, *baseline, purity=purity))
    return rows


def y_rows(copies, sigma=0.01, purity=PURITY):
    return [hemizygous_row('Y', start, end, copies, sigma=sigma, purity=purity) for start, end in Y_CALLABLE_BLOCKS_HG19]


def clonal_grid():
    grid = np.exp(-((np.linspace(0, 1, 101) - 1.) ** 2) / (2 * 0.04 ** 2))
    return grid / grid.sum()


def subclonal_grid(center=0.3):
    grid = np.exp(-((np.linspace(0, 1, 101) - center) ** 2) / (2 * 0.04 ** 2))
    return grid / grid.sum()


def clonal_muts(chromosomes, n=5, start=2000000):
    return [(c, start + i * 1000, 1, 'A') for c in chromosomes for i in range(n)]


class TimingCase(unittest.TestCase):
    genome = Genome('hg19')

    def setUp(self):
        random.seed(1)
        np.random.seed(1)
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp)
        self.cwd = os.getcwd()

    def write_inputs(self, sample, seg_rows, muts):
        seg_path = os.path.join(self.tmp, sample + '.segtab.allelic.completed.txt.gz')
        with gzip.open(seg_path, 'wt') as handle:
            pd.DataFrame(seg_rows, columns=SEG_COLUMNS).to_csv(handle, sep='\t', index=False)
        maf_path = os.path.join(self.tmp, sample + '.mut_ccfs.txt')
        grids = [f'preDP_ccf_{i / 100.0}' for i in range(101)]
        header = ['Patient_ID', 'Sample_ID', 'Hugo_Symbol', 'Chromosome', 'Start_position', 'Reference_Allele',
                  'Tumor_Seq_Allele', 't_ref_count', 't_alt_count', 'Protein_change', 'Variant_Classification',
                  'Variant_Type', 'Cluster_Assignment'] + grids
        # A clonal and a distant subclonal cluster, so CCF assignment means something. A mutation may give its counts.
        muts = list(muts) or clonal_muts(['1']) + [('2', 2000000 + i, 2, 'A') for i in range(5)]
        with open(maf_path, 'w') as handle:
            handle.write('\t'.join(header) + '\n')
            for index, mut in enumerate(muts):
                chrom, pos, cluster, gene = mut[:4]
                ref, alt = mut[4:6] if len(mut) > 4 else (40, 40 if cluster == 1 else 12)
                grid = clonal_grid() if cluster == 1 else subclonal_grid()
                row = ['P1', sample, gene, chrom, pos, 'A', 'T', ref, alt, f'p.A{index}T', 'Missense_Mutation', 'SNP', cluster] + list(grid)
                handle.write('\t'.join(map(str, row)) + '\n')
        return maf_path, seg_path

    def run_samples(self, samples, focal=None, options=None, artifact=None, max_length=5000000):
        """samples: list of (name, seg_rows, muts, purity)."""
        options = options or TimingOptions()
        patient = Patient(indiv_name='P1', genome_build='hg19')
        if options.male:
            patient.hemizygous_chroms = {'X', 'Y'}
        for name, seg_rows, muts, purity in samples:
            maf_path, seg_path = self.write_inputs(name, seg_rows, muts)
            patient.addSample(maf_path, name, input_type='post-clustering', timepoint_value=0, grid_size=101,
                              _additional_muts=None, seg_file=seg_path, purity=purity)
        patient.get_arm_level_cn_events()
        if focal is not None:
            focal_path = os.path.join(self.tmp, 'focal.txt')
            with open(focal_path, 'w') as handle:
                for chrom, start, end, kind in focal:
                    handle.write(f'{chrom}\t{start}\t{end}\t{kind}\t0.01\tregion\n')
            homdels = None
            if options.segment_homdels:
                homdels = dict(artifact_file=artifact, max_length=max_length, max_copies=.5, min_z=3.)
            patient.add_focal_cn_events(focal_regions=focal_path, segment_homdels=homdels, continuous=options.continuous_cn)
        patient.preprocess_samples()
        engine = TimingEngine(patient, min_supporting_muts=3, options=options)
        engine.time_events()
        return engine

    def run_engine(self, seg_rows, muts=(), focal=None, options=None, artifact=None, sample='S1', purity=PURITY, max_length=5000000):
        return self.run_samples([(sample, seg_rows, muts, purity)], focal, options, artifact, max_length)

    @staticmethod
    def events(engine):
        return {name: events[0] for name, events in engine.all_cn_events.items() if events}

    @staticmethod
    def pi_mean(dist):
        return float((np.asarray(dist) * np.linspace(0, 1, 101)).sum())


class MaleBaseline(TimingCase):
    def test_without_options_a_normal_male_x_is_called_lost(self):
        rows = whole_genome(self.genome, male=True) + y_rows(1.)
        events = self.events(self.run_engine(rows))
        self.assertIn('loss_Xp', events)
        self.assertNotIn('loss_Y', events)

    def test_normal_male_has_no_sex_chromosome_events(self):
        rows = whole_genome(self.genome, male=True) + y_rows(1.)
        events = self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON)))
        self.assertEqual({name for name in events if 'X' in name or 'Y' in name}, set())

    def test_female_x_is_diploid(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != 'X']
        rows += [seg_row('X', 1000, 60000000, 0, 1), seg_row('X', 62000000, 155000000, 0, 1)]
        events = self.events(self.run_engine(rows, options=TimingOptions(sex='XX', **ALL_ON)))
        self.assertIn('loss_Xp', events)  # a real loss in a female

    def test_sex_must_be_known(self):
        for sex in (None, '', 'Unknown', float('nan')):
            with self.assertRaises(ValueError):
                TimingOptions(sex=sex, sex_aware=True)
        for sex, expected in ((' XY', True), ('Male', True), ('xx', False), ('Female', False)):
            self.assertEqual(TimingOptions(sex=sex, sex_aware=True).male, expected)


class LossOfY(TimingCase):
    def male(self, y_cn, **kwargs):
        rows = whole_genome(self.genome, male=True) + y_rows(y_cn, **kwargs)
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

    def test_ccf_near_the_cutoff_stays_unknown(self):
        event = self.events(self.male(0.2))['loss_Y']
        self.assertEqual(event.timing_status, 'ambiguous_clonality')
        self.assertIsNone(event.pi_dist)
        self.assertFalse(event.is_clonal)

    def test_noisy_y_widens_the_ambiguous_band(self):
        self.assertEqual(self.events(self.male(0.45))['loss_Y'].timing_status, 'subclonal')
        self.assertEqual(self.events(self.male(0.45, sigma=2.))['loss_Y'].timing_status, 'ambiguous_clonality')

    def test_small_deviation_is_no_event(self):
        self.assertNotIn('loss_Y', self.events(self.male(0.9)))

    def test_no_y_limit_follows_the_copy_ratio_at_low_purity(self):
        # At purity 0.3 the copy-ratio floor of tumors without Y reads as 0.25 copies.
        purity = .3
        rows = whole_genome(self.genome, male=True, purity=purity) + y_rows(.25, purity=purity)
        event = self.events(self.run_engine(rows, purity=purity, options=TimingOptions(sex='XY', **ALL_ON)))['loss_Y']
        self.assertEqual(event.allelic_cn, 0.)
        self.assertEqual(event.ccf_hat, 1.)
        self.assertTrue(event.is_clonal)

    def test_too_little_callable_territory_gives_no_event(self):
        rows = whole_genome(self.genome, male=True) + y_rows(0.)[:1]
        self.assertNotIn('loss_Y', self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON))))

    def test_only_callable_territory_counts(self):
        # Long segments outside the blocks must not move the median, as cross-mapped signal would.
        rows = whole_genome(self.genome, male=True) + y_rows(1.)
        rows += [hemizygous_row('Y', 2650000, 6400000, 0.), hemizygous_row('Y', 24000000, 28790000, 0.)]
        self.assertNotIn('loss_Y', self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON))))

    def test_median_not_mean(self):
        # One small block far off must not move the value, though it moves a mean.
        rows = whole_genome(self.genome, male=True) + y_rows(1.)[:2] + [hemizygous_row('Y', 20800000, 23810000, 0.)]
        self.assertNotIn('loss_Y', self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON))))

    def test_wide_y_segment_gives_one_event(self):
        rows = whole_genome(self.genome, male=True) + [hemizygous_row('Y', 2650000, 28790000, 0.)]
        events = self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON)))
        self.assertEqual([name for name in events if 'Y' in name], ['loss_Y'])


class WgdRule(TimingCase):
    """A tetraploid male with timeable mutations on 2/2 arms."""
    CHROMS = [str(i) for i in range(1, 13)]

    def muts(self, major_fraction=1.):
        """Clonal SNVs on 2/2 arms, with this share at multiplicity 2."""
        depth, muts = 80, []
        # Expected VAF at purity 0.8 and total copy number 4: m=1 about 0.222, m=2 about 0.444.
        for c in self.CHROMS:
            for arm_start in (2000000, 140000000 if c in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12') else 2000000):
                for i in range(10):
                    alt = int(depth * (0.444 if i < 10 * major_fraction else 0.222))
                    muts.append((c, arm_start + i * 1000, 1, 'A', depth - alt, alt))
        return muts

    def engine(self, y_cn=None, extra=(), major_fraction=1., options=None, skip=()):
        rows = whole_genome(self.genome, male=True, baseline=(2, 2), skip=skip)
        rows += [] if y_cn is None else y_rows(y_cn)
        rows += list(extra)
        return self.run_engine(rows, muts=self.muts(major_fraction), options=options or TimingOptions(sex='XY', **ALL_ON))

    def test_wgd_is_called_and_male_x_is_not_an_event(self):
        engine = self.engine(y_cn=2.)
        self.assertIsNotNone(engine.WGD)
        self.assertEqual({name for name in self.events(engine) if 'X' in name or 'Y' in name}, set())

    def test_chromosome_y_does_not_support_wgd(self):
        engine = self.engine(y_cn=2.)
        for sample in engine.sample_list:
            self.assertFalse([s for s in sample.WGD.supporting_arm_states if s.cn_type_prefix == 'Chrom'])

    def test_each_arm_counts_once_for_wgd(self):
        engine = self.engine(y_cn=2.)
        self.assertLessEqual(engine.WGD.n_timed_gains, len(self.CHROMS) * 2 + 1)

    def test_subclonal_mutations_do_not_time_a_wgd(self):
        rows = whole_genome(self.genome, male=True, baseline=(2, 2)) + y_rows(2.)
        muts = [(c, 2000000 + i * 1000, 2, 'A') for c in self.CHROMS for i in range(10)] + [('1', 100, 1, 'A')] * 3
        engine = self.run_engine(rows, muts=muts, options=TimingOptions(sex='XY', **ALL_ON))
        self.assertEqual(engine.WGD.n_timed_gains, 0)
        self.assertEqual(engine.WGD.timing_status, 'clonal_unresolved')
        self.assertIsNone(engine.WGD.pi_dist)

    def test_ambiguous_y_after_wgd_stays_unknown(self):
        event = self.events(self.engine(y_cn=1.2, major_fraction=.5))['loss_Y']
        self.assertEqual(event.timing_status, 'ambiguous_clonality')
        self.assertIsNone(event.pi_dist)

    def test_two_two_arms_time_the_wgd(self):
        # Mutations on both copies arose before the doubling, so mostly multiplicity 2 means a late doubling.
        late = self.engine(y_cn=2., major_fraction=1.)
        early = self.engine(y_cn=2., major_fraction=.2)
        self.assertGreater(late.WGD.n_timed_gains, 10)
        self.assertEqual(early.WGD.timing_status, 'estimated')
        self.assertGreater(self.pi_mean(late.WGD.pi_dist), self.pi_mean(early.WGD.pi_dist) + 0.2)

    def test_y_zero_precedes_wgd(self):
        engine = self.engine(y_cn=0., major_fraction=.5)
        event = self.events(engine)['loss_Y']
        self.assertEqual(event.allelic_cn, 0.)
        self.assertEqual(event.timing_status, 'estimated')
        self.assertLess(self.pi_mean(event.pi_dist), self.pi_mean(engine.WGD.pi_dist))

    def test_y_one_follows_wgd(self):
        engine = self.engine(y_cn=1., major_fraction=.5)
        event = self.events(engine)['loss_Y']
        self.assertEqual(event.allelic_cn, 1.)
        self.assertTrue(event.is_clonal)
        self.assertGreater(self.pi_mean(event.pi_dist), self.pi_mean(engine.WGD.pi_dist))

    def test_y_between_zero_and_one_cannot_precede_wgd(self):
        # Cells still hold Y, so a clonal loss before doubling is impossible. Rounding must not decide this.
        for y_cn in (0.3, 0.49, 0.51, 0.8):
            engine = self.engine(y_cn=y_cn, major_fraction=.5)
            event = self.events(engine)['loss_Y']
            self.assertEqual(event.allelic_cn, 1., msg=str(y_cn))
            self.assertGreater(self.pi_mean(event.pi_dist), self.pi_mean(engine.WGD.pi_dist), msg=str(y_cn))

    def test_y_above_one_gives_the_ccf_of_one_copy(self):
        event = self.events(self.engine(y_cn=1.5, major_fraction=.5))['loss_Y']
        self.assertAlmostEqual(event.ccf_hat, .5)
        self.assertFalse(event.is_clonal)

    def test_wgd_without_evidence_is_not_presented_as_timed(self):
        # No mutations on gained arms: the doubling cannot be timed, so nothing can be timed against it.
        rows = whole_genome(self.genome, male=True, baseline=(2, 2)) + y_rows(0.)
        muts = [('21', 15000000 + i, 1, 'A') for i in range(4)]
        engine = self.run_engine(rows, muts=muts, options=TimingOptions(sex='XY', **ALL_ON))
        if engine.WGD is None:
            self.skipTest('No WGD call without mutations on the arms')
        self.assertEqual(engine.WGD.timing_status, 'clonal_unresolved')
        self.assertIsNone(engine.WGD.pi_dist)
        self.assertEqual(self.events(engine)['loss_Y'].timing_status, 'clonal_unresolved')


class ContinuousArmEvents(TimingCase):
    def arm_engine(self, q1, modal_a1, q2=None, modal_a2=1, continuous=True, muts=(), purity=PURITY):
        rows = [r for r in whole_genome(self.genome, male=False, purity=purity) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, modal_a1, modal_a2, q1=q1, q2=q2, purity=purity),
                 seg_row('17', 25000000, 81000000, 1, 1, purity=purity)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=continuous, clonal_unresolved=True)
        return self.run_engine(rows, muts=muts, options=options, purity=purity)

    def test_rescaled_columns_are_not_read_on_autosomes(self):
        # As in real NBE10 13q: modal 0/2, rescaled a1 0, while hscr places a1 at 0.30.
        event = self.events(self.arm_engine(0.3, 0, q2=2., modal_a2=2))['loss_17p']
        self.assertAlmostEqual(event.ccf_hat, 0.7)
        self.assertFalse(event.is_clonal)

    def test_subclonal_arm_loss_rounded_to_normal_is_called(self):
        event = self.events(self.arm_engine(0.6, 1))['loss_17p']
        self.assertAlmostEqual(event.ccf_hat, 0.4)
        self.assertEqual(event.timing_status, 'subclonal')

    def test_subclonal_gain_is_a_gain(self):
        event = self.events(self.arm_engine(1., 1, q2=1.4))['gain_17p']
        self.assertEqual(event.Type, 'Arm_gain')
        self.assertAlmostEqual(event.ccf_hat, 0.4)
        self.assertFalse(event.is_clonal)

    def test_integer_state_misses_it(self):
        self.assertNotIn('loss_17p', self.events(self.arm_engine(0.6, 1, continuous=False)))

    def test_clonal_loss_stays_clonal(self):
        event = self.events(self.arm_engine(0.03, 0))['loss_17p']
        self.assertAlmostEqual(event.ccf_hat, 0.97)
        self.assertTrue(event.is_clonal)
        self.assertEqual(event.timing_status, 'clonal_unresolved')

    def test_noise_near_normal_is_no_event(self):
        self.assertNotIn('loss_17p', self.events(self.arm_engine(0.9, 1)))

    def test_missing_hscr_falls_back_to_the_modal_state(self):
        event = self.events(self.arm_engine(float('nan'), 0))['loss_17p']
        self.assertEqual(event.ccf_hat, 1.)

    def test_ccf_cap_holds_for_high_gains(self):
        self.assertEqual(self.events(self.arm_engine(1., 1, q2=4., modal_a2=4))['gain_17p'].ccf_hat, 1.)

    def test_clonality_is_judged_per_allele(self):
        engine = self.arm_engine(1.3, 1, q2=2., modal_a2=2)
        events = {e.allele: e for e in engine.all_cn_events['gain_17p']}
        self.assertFalse(events['a1'].is_clonal)
        self.assertTrue(events['a2'].is_clonal)
        engine = self.arm_engine(0.05, 0, q2=0.6, modal_a2=1)
        events = {e.allele: e for e in engine.all_cn_events['loss_17p']}
        self.assertTrue(events['a1'].is_clonal)
        self.assertFalse(events['a2'].is_clonal)

    def test_segments_rounding_differently_are_one_state(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 11000000, 1, 1, q1=.52), seg_row('17', 11000001, 22000000, 0, 1, q1=.48),
                 seg_row('17', 25000000, 81000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
        event = self.events(self.run_engine(rows, options=options))['loss_17p']
        self.assertAlmostEqual(event.ccf_hat, .5, delta=.03)
        self.assertFalse(event.is_clonal)

    def test_balanced_segments_give_one_event(self):
        # ACEseq saw no imbalance, so both alleles carry half of a subclonal gain. This is one event with its full CCF.
        engine = self.arm_engine(1.25, 1, q2=1.25)
        self.assertEqual(len(engine.all_cn_events['gain_17p']), 1)
        self.assertAlmostEqual(engine.all_cn_events['gain_17p'][0].ccf_hat, .5)

    def test_plurality_state_is_called_when_the_median_falls_elsewhere(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 25000000, 50000000, 1, 3), seg_row('17', 50000001, 66000000, 1, 2), seg_row('17', 66000001, 81000000, 1, 1),
                 seg_row('17', 1000, 22000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
        event = self.events(self.run_engine(rows, options=options))['gain_17q']
        self.assertEqual(event.allele, 'a2')
        self.assertEqual(event.allelic_cn, 3.)

    def test_no_state_with_enough_of_the_arm_gives_no_call(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 25000000, 43000000, 1, 3), seg_row('17', 43000001, 62000000, 1, 2), seg_row('17', 62000001, 81000000, 1, 5),
                 seg_row('17', 1000, 22000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
        self.assertNotIn('gain_17q', self.events(self.run_engine(rows, options=options)))

    def test_values_apart_by_more_than_the_tolerance_are_different_states(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 12500000, 1, 1, q1=1.), seg_row('17', 12500001, 22000000, 1, 1, q1=.6), seg_row('17', 25000000, 81000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
        self.assertNotIn('loss_17p', self.events(self.run_engine(rows, options=options)))

    def test_noise_decides_between_subclonal_and_unknown(self):
        def status(sigma):
            rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
            rows += [seg_row('17', 1000, 22000000, 0, 1, q1=.45, sigma=sigma), seg_row('17', 25000000, 81000000, 1, 1)]
            options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
            return self.events(self.run_engine(rows, options=options))['loss_17p'].timing_status
        self.assertEqual(status(.01), 'subclonal')
        self.assertEqual(status(.5), 'ambiguous_clonality')

    def test_two_competing_states_give_no_call(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 11000000, 1, 1), seg_row('17', 11000001, 22000000, 0, 1),
                 seg_row('17', 25000000, 81000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
        self.assertNotIn('loss_17p', self.events(self.run_engine(rows, options=options)))

    def test_median_is_length_weighted(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 2000000, 1, 1), seg_row('17', 2000001, 3000000, 1, 1),
                 seg_row('17', 3000001, 22000000, 0, 1, q1=.3), seg_row('17', 25000000, 81000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
        self.assertAlmostEqual(self.events(self.run_engine(rows, options=options))['loss_17p'].ccf_hat, .7)

    def test_gain_timing_uses_the_modal_state(self):
        # Modal 1/2 with a1 at 0.6: the reads must still time the gain, by multiplicities on total copy number 3.
        depth = 100
        def muts(m):
            alt = int(depth * m * PURITY / (3 * PURITY + 2 * (1 - PURITY)))
            return [('17', 5000000 + i * 1000, 1, 'A', depth - alt, alt) for i in range(12)]
        early = self.arm_engine(0.97, 1, q2=2., modal_a2=2, muts=muts(2))
        late = self.arm_engine(0.97, 1, q2=2., modal_a2=2, muts=muts(1))
        off_early = self.arm_engine(0.6, 1, q2=2., modal_a2=2, muts=muts(2))
        off_late = self.arm_engine(0.6, 1, q2=2., modal_a2=2, muts=muts(1))
        self.assertGreater(self.pi_mean(self.events(early)['gain_17p'].pi_dist), self.pi_mean(self.events(late)['gain_17p'].pi_dist) + .3)
        self.assertGreater(self.pi_mean(self.events(off_early)['gain_17p'].pi_dist), self.pi_mean(self.events(off_late)['gain_17p'].pi_dist) + .3)

    def test_subclonal_supporters_do_not_time_a_gain(self):
        muts = [('17', 5000000 + i * 1000, 2, 'A') for i in range(8)] + clonal_muts(['1'])
        event = self.events(self.arm_engine(1., 1, q2=2., modal_a2=2, muts=muts))['gain_17p']
        self.assertEqual(event.timing_status, 'clonal_unresolved')
        self.assertIsNone(event.pi_dist)


class ClusterComparator(TimingCase):
    def test_large_clonal_cluster_does_not_overrule_a_clonal_ccf(self):
        # A narrow clonal cluster scores CCF 0.92 below a small subclonal cluster. Copy number must decide.
        rows = [r for r in whole_genome(self.genome, male=True) if r[0] != '17'] + y_rows(.08)
        rows += [seg_row('17', 1000, 22000000, 0, 1, q1=.08), seg_row('17', 25000000, 81000000, 1, 1)]
        muts = clonal_muts(map(str, range(1, 21)), n=100) + [('2', 3000000 + i * 1000, 2, 'A') for i in range(30)]
        events = self.events(self.run_engine(rows, muts=muts, options=TimingOptions(sex='XY', **ALL_ON)))
        self.assertAlmostEqual(events['loss_17p'].ccf_hat, .92)
        for name in ('loss_17p', 'loss_Y'):
            self.assertTrue(events[name].is_clonal)
            self.assertEqual(events[name].timing_status, 'clonal_unresolved')

    def test_large_clonal_cluster_with_a_homdel(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '9']
        size = self.genome.CHROM_DICT['9']
        rows += [seg_row('9', 1000, 21899999, 1, 1), seg_row('9', 21900000, 22100000, 0, 0, q1=.04, q2=.04), seg_row('9', 22100001, size - 1000, 1, 1)]
        muts = clonal_muts(map(str, range(1, 21)), n=100) + [('2', 3000000 + i * 1000, 2, 'A') for i in range(30)]
        engine = self.run_engine(rows, muts=muts, focal=[('9', 21000000, 23000000, 'DEL')], options=TimingOptions(sex='XX', **ALL_ON))
        self.assertTrue(self.events(engine)['Focal_homdel_9:21000000-23000000'].is_clonal)


class FocalCase(TimingCase):
    AMP = ('2', 16000000, 16200000, 'AMP')
    DEL = ('9', 8000000, 11000000, 'DEL')

    def with_segments(self, chrom, *middle, male=False, baseline=(1, 1), purity=PURITY):
        rows = [r for r in whole_genome(self.genome, male=male, baseline=baseline, purity=purity) if r[0] != chrom]
        size = self.genome.CHROM_DICT[chrom]
        start = middle[0][1]
        end = middle[-1][2]
        return rows + [seg_row(chrom, 1000, start - 1, *baseline, purity=purity)] + list(middle) + [seg_row(chrom, end + 1, size - 1000, *baseline, purity=purity)]

    def focal(self, rows, regions, options=None):
        return self.events(self.run_engine(rows, focal=regions, options=options or TimingOptions(sex='XX', **ALL_ON)))


class FocalEvents(FocalCase):
    def test_amplification_with_modal_one_one_is_seen(self):
        # As in real MYCN amplifications, the modal columns stay at 1/1 while total copy number is in the hundreds.
        rows = self.with_segments('2', seg_row('2', 15900000, 16300000, 1, 1, q1=1., q2=99.))
        event = self.focal(rows, [self.AMP])['Focal_gain_2:16000000-16200000']
        self.assertTrue(event.is_clonal)
        self.assertEqual(event.ccf_hat, 1.)

    def test_broad_arm_gain_is_not_focal(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '2']
        rows += [seg_row('2', 1000, 90000000, 1, 2), seg_row('2', 96000000, 243000000, 1, 1)]
        self.assertNotIn('Focal_gain_2:16000000-16200000', self.focal(rows, [self.AMP]))

    def test_gain_above_a_gained_arm_is_focal(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '2']
        rows += [seg_row('2', 1000, 15899999, 1, 2), seg_row('2', 15900000, 16300000, 1, 6), seg_row('2', 16300001, 90000000, 1, 2),
                 seg_row('2', 96000000, 243000000, 1, 1)]
        self.assertIn('Focal_gain_2:16000000-16200000', self.focal(rows, [self.AMP]))

    def test_normal_segment_beside_a_partial_arm_loss_is_not_a_gain(self):
        # As for CCND1 beside a distal 11q deletion: the arm median is low, but the window matches its near flank.
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '11']
        size = self.genome.CHROM_DICT['11']
        rows += [seg_row('11', 1000, 50000000, 1, 1), seg_row('11', 54000000, 69000000, 1, 1), seg_row('11', 69000001, 69500000, 1, 1),
                 seg_row('11', 69500001, 72000000, 1, 1), seg_row('11', 72000001, size - 1000, 0, 1)]
        self.assertNotIn('Focal_gain_11:69455855-69469242', self.focal(rows, [('11', 69455855, 69469242, 'AMP')]))

    def test_amplicon_of_several_segments_is_measured_from_the_baseline(self):
        rows = self.with_segments('2', seg_row('2', 15900000, 16000000, 1, 1, q1=1., q2=39.), seg_row('2', 16000001, 16200000, 1, 1, q1=1., q2=99.),
                                  seg_row('2', 16200001, 16300000, 1, 1, q1=1., q2=59.))
        self.assertEqual(self.focal(rows, [self.AMP])['Focal_gain_2:16000000-16200000'].ccf_hat, 1.)

    def test_amplicon_that_fills_its_flanks_is_still_seen(self):
        # The amplicon spans most of 4 Mb, so the flank medians are amplified too. It stands out from the arm instead.
        rows = self.with_segments('2', *[seg_row('2', 14000000 + i * 500000, 14000000 + (i + 1) * 500000 - 1, 1, 1, q1=1., q2=60. + i) for i in range(9)])
        region = ('2', 16000000, 16200000, 'AMP')
        self.assertEqual(self.focal(rows, [region])['Focal_gain_2:16000000-16200000'].ccf_hat, 1.)

    def test_high_amplification_survives_its_large_noise(self):
        # Segment noise scales with copy ratio. A real amplicon of 100 copies can have sigma near 15.
        rows = self.with_segments('2', seg_row('2', 15900000, 16300000, 1, 1, q1=1., q2=99., sigma=15.))
        self.assertIn('Focal_gain_2:16000000-16200000', self.focal(rows, [self.AMP]))

    def test_subclonal_focal_gain_has_its_ccf(self):
        rows = self.with_segments('2', seg_row('2', 15900000, 16300000, 1, 1, q1=1., q2=1.6))
        event = self.focal(rows, [self.AMP])['Focal_gain_2:16000000-16200000']
        self.assertAlmostEqual(event.ccf_hat, .6, places=2)
        self.assertEqual(event.timing_status, 'subclonal')

    def test_small_step_is_no_event(self):
        rows = self.with_segments('2', seg_row('2', 15900000, 16300000, 1, 1, q1=1., q2=1.3))
        self.assertNotIn('Focal_gain_2:16000000-16200000', self.focal(rows, [self.AMP]))

    def test_noisy_segment_is_no_event(self):
        rows = self.with_segments('2', seg_row('2', 15900000, 16300000, 1, 1, q1=1., q2=2., sigma=.5))
        self.assertNotIn('Focal_gain_2:16000000-16200000', self.focal(rows, [self.AMP]))

    def test_focal_loss_is_clonal_with_ccf_one(self):
        rows = self.with_segments('9', seg_row('9', 9000000, 9500000, 0, 1))
        event = self.focal(rows, [self.DEL])['Focal_loss_9:8000000-11000000']
        self.assertEqual(event.ccf_hat, 1.)
        self.assertTrue(event.is_clonal)

    def test_loss_from_a_doubled_arm(self):
        rows = self.with_segments('9', seg_row('9', 9000000, 9500000, 1, 2), baseline=(2, 2))
        depth, muts = 80, []
        for c in map(str, range(1, 9)):
            for i in range(10):
                alt = int(depth * (.444 if i % 2 else .222))
                muts += [(c, 2000000 + i * 1000, 1, 'A', depth - alt, alt), (c, 140000000 + i * 1000, 1, 'A', depth - alt, alt)]
        engine = self.run_engine(rows, muts=muts, focal=[self.DEL], options=TimingOptions(sex='XX', **ALL_ON))
        self.assertIsNotNone(engine.WGD)
        event = self.events(engine)['Focal_loss_9:8000000-11000000']
        self.assertEqual(event.allelic_cn, 1.)
        self.assertGreater(self.pi_mean(event.pi_dist), self.pi_mean(engine.WGD.pi_dist))

    def test_partial_deletion_on_male_x_is_a_loss(self):
        name = 'Focal_loss_X:76000000-78000000'
        size = self.genome.CHROM_DICT['X']
        rows = [r for r in whole_genome(self.genome, male=True) if r[0] != 'X'] + y_rows(1.)
        rows += [hemizygous_row('X', 1000, 76899999, 1.), hemizygous_row('X', 76900000, 77100000, .4), hemizygous_row('X', 77100001, size - 1000, 1.)]
        options = TimingOptions(sex='XY', **ALL_ON)
        events = self.focal(rows, [('X', 76000000, 78000000, 'DEL')], options)
        self.assertIn('Focal_homdel_X:76000000-78000000', events)
        rows[-2] = hemizygous_row('X', 76900000, 77100000, .6)
        event = self.focal(rows, [('X', 76000000, 78000000, 'DEL')], options)[name]
        self.assertAlmostEqual(event.ccf_hat, .4, places=2)
        self.assertEqual(event.timing_status, 'subclonal')


class Hemizygous(TimingCase):
    def test_multiplicity_uses_one_normal_copy(self):
        # Low purity: with two assumed normal copies, m=1 on 0/2 reads as far earlier than it is.
        purity, depth = .3, 200
        expected = purity / (2 * purity + (1 - purity))
        alt = int(depth * expected)
        rows = [r for r in whole_genome(self.genome, male=True, purity=purity) if r[0] != 'X']
        rows += [hemizygous_row('X', 1000, 58000000, 2., purity=purity), hemizygous_row('X', 62000000, 155000000, 1., purity=purity)] + y_rows(1., purity=purity)
        muts = [('X', 5000000 + i * 1000, 1, 'A', depth - alt, alt) for i in range(15)]
        engine = self.run_engine(rows, muts=muts, purity=purity, options=TimingOptions(sex='XY', **ALL_ON))
        event = self.events(engine)['gain_Xp']
        self.assertLess(self.pi_mean(event.pi_dist), .35)

    def test_male_x_normal_copy_is_not_a_deletion(self):
        region = ('X', 76000000, 78000000, 'DEL')
        rows = whole_genome(self.genome, male=True)
        events = self.events(self.run_engine(rows, focal=[region], options=TimingOptions(sex='XY', **ALL_ON)))
        self.assertEqual([name for name in events if 'X' in name], [])


class Multisample(TimingCase):
    def test_shared_subclonal_arm_loss_keeps_ccf_per_sample(self):
        def rows(q1):
            base = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
            return base + [seg_row('17', 1000, 22000000, 1, 1, q1=q1), seg_row('17', 25000000, 81000000, 1, 1)]
        muts = clonal_muts(['1'])
        engine = self.run_samples([('A', rows(.6), muts, PURITY), ('B', rows(.65), muts, PURITY)],
                                  options=TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True))
        event = self.events(engine)['loss_17p']
        self.assertFalse(event.is_clonal)
        self.assertAlmostEqual(event.ccf_hat, .35)

    def test_event_private_to_one_sample_is_dropped(self):
        def rows(q1, q2):
            base = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
            return base + [seg_row('17', 1000, 22000000, 1, 1, q1=q1, q2=q2), seg_row('17', 25000000, 81000000, 1, 1)]
        muts = clonal_muts(['1'])
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
        engine = self.run_samples([('A', rows(1., 1.), muts, PURITY), ('B', rows(1., 1.5), muts, PURITY)], options=options)
        self.assertNotIn('gain_17p', self.events(engine))
        # Each sample has an event on its own allele. Neither is shared, so neither reaches the patient level.
        engine = self.run_samples([('A', rows(.6, 1.), muts, PURITY), ('B', rows(1., 1.6), muts, PURITY)], options=options)
        self.assertEqual({n for n in self.events(engine) if '17p' in n}, set())

    def test_wgd_in_one_sample_only_creates_no_false_gains(self):
        muts = clonal_muts(map(str, range(1, 9)), n=10) + clonal_muts(map(str, range(1, 9)), n=10, start=140000000)
        doubled = whole_genome(self.genome, male=False, baseline=(2, 2))
        plain = whole_genome(self.genome, male=False)
        engine = self.run_samples([('A', doubled, muts, PURITY), ('B', plain, muts, PURITY)],
                                  options=TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True))
        self.assertIsNone(engine.WGD)
        self.assertEqual({n for n in self.events(engine) if n.startswith('gain_')}, set())

    def test_homdel_private_to_one_sample_is_subclonal(self):
        def rows(copies, a):
            base = [r for r in whole_genome(self.genome, male=False) if r[0] != '9']
            size = self.genome.CHROM_DICT['9']
            return base + [seg_row('9', 1000, 21899999, 1, 1), seg_row('9', 21900000, 22100000, a, a, q1=copies / 2, q2=copies / 2),
                           seg_row('9', 22100001, size - 1000, 1, 1)]
        muts = clonal_muts(['1'])
        engine = self.run_samples([('A', rows(.02, 0), muts, PURITY), ('B', rows(2., 1), muts, PURITY)],
                                  focal=[('9', 21000000, 23000000, 'DEL')], options=TimingOptions(sex='XX', **ALL_ON))
        event = self.events(engine)['Focal_homdel_9:21000000-23000000']
        self.assertFalse(event.is_clonal)
        self.assertEqual(event.timing_status, 'subclonal')


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
            self.timed('ambiguous', None, 'ambiguous_clonality'),
        ])
        self.assertEqual(comps[('clonal_loss', 'subclonal')], (1.0, 0.0, 0.0))
        self.assertEqual(comps[('clonal_snv', 'subclonal')], (1.0, 0.0, 0.0))
        self.assertEqual(comps[('clonal_loss', 'clonal_snv')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('clonal_loss', 'timed')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('clonal_loss', 'legacy')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('clonal_loss', 'ambiguous')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('subclonal', 'ambiguous')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('subclonal', 'legacy')], (0.0, 0.0, 1.0))
        self.assertGreater(comps[('subclonal', 'timed')][1], .99)

    def test_order_of_the_pair_does_not_matter(self):
        delta = np.zeros(101)
        delta[100] = 1.
        comps = spt.compare_event_list([self.timed('subclonal', delta, 'subclonal'), self.timed('clonal', None, 'clonal_unresolved')])
        self.assertEqual(comps[('subclonal', 'clonal')], (0.0, 1.0, 0.0))

    def test_without_statuses_nothing_changes(self):
        delta = np.zeros(101)
        delta[100] = 1.
        self.assertEqual(spt.compare_event_list([self.timed('a'), self.timed('b', delta)])[('a', 'b')], (0.0, 0.0, 1.0))

    def table(self, rows, options):
        engine = self.run_engine(rows, options=options)
        os.chdir(self.tmp)
        self.addCleanup(os.chdir, self.cwd)
        PhylogicOutput().write_timing_tsv(engine)
        return engine, pd.read_csv('P1.timing.tsv', sep='\t', low_memory=False)

    def test_timing_table_round_trip(self):
        rows = whole_genome(self.genome, male=True) + y_rows(0.)
        rows = [r for r in rows if r[0] != '17'] + [seg_row('17', 1000, 22000000, 1, 1, q1=.6), seg_row('17', 25000000, 81000000, 1, 1)]
        _, table = self.table(rows, TimingOptions(sex='XY', **ALL_ON))
        y = table[table['Event Name'] == 'loss_Y'].iloc[0]
        self.assertEqual(y['timing_status'], 'clonal_unresolved')
        self.assertEqual(str(y['is_clonal']), 'True')
        snvs = table[table['Start_position'].notna()]
        self.assertEqual(set(snvs['timing_status']), {'clonal_unresolved', 'subclonal'})
        for events in (spt.read_timed_events('P1.timing.tsv', compare_all=True), self.selected(['loss_Y', 'loss_17p'])):
            comps = spt.compare_event_list(events)
            pair = ('loss_Y', 'loss_17p') if ('loss_Y', 'loss_17p') in comps else ('loss_17p', 'loss_Y')
            self.assertEqual(comps[pair], (1.0, 0.0, 0.0) if pair[0] == 'loss_Y' else (0.0, 1.0, 0.0))

    def selected(self, names):
        with open('events.txt', 'w') as handle:
            handle.write('\n'.join(names) + '\n')
        return spt.read_timed_events('P1.timing.tsv', event_file='events.txt')

    def test_earliest_allele_represents_a_shared_name(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, 0, 1, q1=.03, q2=.6), seg_row('17', 25000000, 81000000, 1, 1)]
        _, table = self.table(rows, TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True))
        self.assertEqual(sorted(table[table['Event Name'] == 'loss_17p']['timing_status']), ['clonal_unresolved', 'subclonal'])
        events = {e.event_name: e for e in spt.read_timed_events('P1.timing.tsv', compare_all=True)}
        self.assertEqual(events['loss_17p'].timing_status, 'clonal_unresolved')

    def test_earlier_of_two_timed_gains_represents_a_shared_name(self):
        grid = np.linspace(0, 1, 101)
        def row(center):
            pi = np.exp(-((grid - center) ** 2) / .002)
            return list(pi / pi.sum())
        columns = ['Event Name'] + [f'pi_{float(x) / 100}' for x in range(101)] + ['timing_status']
        pd.DataFrame([['gain_17q'] + row(.8) + ['estimated'], ['gain_17q'] + row(.2) + ['estimated'], ['gain_1q'] + row(.5) + ['estimated']],
                     columns=columns).to_csv(os.path.join(self.tmp, 'dup.tsv'), sep='\t', index=False)
        events = {e.event_name: e for e in spt.read_timed_events(os.path.join(self.tmp, 'dup.tsv'), compare_all=True)}
        self.assertLess(self.pi_mean(events['gain_17q'].pi_dist), .3)

    def test_original_timing_table_keeps_its_columns(self):
        _, table = self.table(whole_genome(self.genome, male=False), TimingOptions())
        self.assertNotIn('timing_status', table.columns)


class SegmentHomdels(TimingCase):
    REGION = ('9', 21000000, 23000000, 'DEL')
    NAME = 'Focal_homdel_9:21000000-23000000'

    def rows(self, copies, sigma=0.01, start=21900000, end=22100000, a1=0, a2=0, n_probes=500):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '9']
        size = self.genome.CHROM_DICT['9']
        rows += [seg_row('9', 1000, start - 1, 1, 1),
                 seg_row('9', start, end, a1, a2, q1=copies / 2., q2=copies / 2., sigma=sigma, n_probes=n_probes),
                 seg_row('9', end + 1, size - 1000, 1, 1)]
        return rows

    def homdel(self, rows, focal=None, options=None, **kwargs):
        options = options or TimingOptions(sex='XX', **ALL_ON)
        return self.events(self.run_engine(rows, focal=focal or [self.REGION], options=options, **kwargs))

    def test_clear_deletion_is_called_clonal_unresolved(self):
        event = self.homdel(self.rows(0.02))[self.NAME]
        self.assertEqual(event.timing_status, 'clonal_unresolved')
        self.assertTrue(event.is_clonal)

    def test_rounding_to_zero_is_not_a_homdel_but_stays_a_loss(self):
        # The whole window is modal 0/0, so the window summary reports a loss on both alleles.
        events = self.homdel(self.rows(0.59, sigma=0.02, start=21000000, end=23000000))
        self.assertNotIn(self.NAME, events)
        self.assertIn('Focal_loss_9:21000000-23000000', events)

    def test_original_definition_accepts_it(self):
        options = TimingOptions(sex='XX', sex_aware=True, clonal_unresolved=True)
        self.assertIn(self.NAME, self.homdel(self.rows(0.59, sigma=0.02, start=21000000, end=23000000), options=options))

    def test_copies_above_half_are_rejected(self):
        self.assertIn(self.NAME, self.homdel(self.rows(0.45, sigma=0.001)))
        self.assertNotIn(self.NAME, self.homdel(self.rows(0.55, sigma=0.001)))

    def test_noise_comes_from_the_segment_on_the_copy_number_scale(self):
        # At slope about 2.5 copies per copy ratio, sigma 0.15 is about 0.4 copies: 0.2 copies cannot be told from 1.
        self.assertIn(self.NAME, self.homdel(self.rows(0.2, sigma=0.01)))
        self.assertNotIn(self.NAME, self.homdel(self.rows(0.2, sigma=0.15)))

    def test_too_long_segment_is_rejected(self):
        self.assertNotIn(self.NAME, self.homdel(self.rows(0.02, start=20000000, end=26000000)))

    def test_too_few_probes_is_rejected(self):
        self.assertNotIn(self.NAME, self.homdel(self.rows(0.02, n_probes=2)))

    def artifact(self, text):
        path = os.path.join(self.tmp, 'artifacts.txt')
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_artifact_region_is_rejected_entirely(self):
        events = self.homdel(self.rows(0.02), artifact=self.artifact('#chromosome\tstart\tend\trec\n9\t21800000\t22200000\t5\n'))
        self.assertNotIn(self.NAME, events)
        self.assertNotIn('Focal_loss_9:21000000-23000000', events)

    def test_homdel_does_not_also_emit_a_loss(self):
        self.assertNotIn('Focal_loss_9:21000000-23000000', self.homdel(self.rows(0.02)))

    def test_small_artifact_overlap_is_kept(self):
        self.assertIn(self.NAME, self.homdel(self.rows(0.02), artifact=self.artifact('9\t21900000\t21910000\t5\n')))

    def test_artifact_regions_add_up(self):
        text = 'chromosome\tstart\tend\n9\t21900000\t21940000\t1\n9\t21940000\t21980000\t1\n9\t21980000\t22020000\t1\n'
        self.assertNotIn(self.NAME, self.homdel(self.rows(0.02), artifact=self.artifact(text)))

    def test_artifact_given_as_one_position_counts(self):
        self.assertNotIn(self.NAME, self.homdel(self.rows(0.02), artifact=self.artifact('9\t21950000\t21950000\t1\n')))

    def test_amplification_region_is_ignored(self):
        self.assertNotIn(self.NAME, self.homdel(self.rows(0.02), focal=[('9', 21000000, 23000000, 'AMP')]))

    def test_modal_state_does_not_decide(self):
        self.assertIn(self.NAME, self.homdel(self.rows(0.02, a1=0, a2=1)))

    def test_ccf_is_relative_to_flanks_and_reaches_the_event(self):
        engine = self.run_engine(self.rows(0.1), focal=[self.REGION], options=TimingOptions(sex='XX', **ALL_ON))
        event = self.events(engine)[self.NAME]
        self.assertAlmostEqual(event.ccf_hat, 0.95, places=2)
        self.assertTrue(event.is_clonal)
        self.assertEqual(event.timing_status, 'clonal_unresolved')
        engine = self.run_engine(self.rows(0.28), focal=[self.REGION], options=TimingOptions(sex='XX', **ALL_ON))
        event = self.events(engine)[self.NAME]
        self.assertAlmostEqual(event.ccf_hat, 0.86, places=2)
        self.assertEqual(event.timing_status, 'ambiguous_clonality')
        engine = self.run_engine(self.rows(0.45), focal=[self.REGION], options=TimingOptions(sex='XX', cn_clonal_min_ccf=.95, **ALL_ON))
        event = self.events(engine)[self.NAME]
        self.assertAlmostEqual(event.ccf_hat, 0.775, places=2)
        self.assertEqual(event.timing_status, 'subclonal')

    def test_flank_lengths_are_clipped_to_the_window(self):
        # A long segment reaches only 200 kb into the 1 Mb flank, so it must not outweigh what lies next to the deletion.
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '9']
        size = self.genome.CHROM_DICT['9']
        rows += [seg_row('9', 1000, 20900000, 4, 4), seg_row('9', 20900001, 21699999, 1, 1),
                 seg_row('9', 21700000, 21900000, 0, 0, q1=.05, q2=.05), seg_row('9', 21900001, 22900000, 1, 1),
                 seg_row('9', 22900001, size - 1000, 1, 1)]
        event = self.events(self.run_engine(rows, focal=[self.REGION], options=TimingOptions(sex='XX', **ALL_ON)))[self.NAME]
        self.assertAlmostEqual(event.ccf_hat, 0.95, places=2)

    def test_male_x_counts_the_single_copy(self):
        name = 'Focal_homdel_X:76000000-78000000'
        size = self.genome.CHROM_DICT['X']
        rows = [r for r in whole_genome(self.genome, male=True) if r[0] != 'X'] + y_rows(1.)
        rows += [hemizygous_row('X', 1000, 76899999, 1.), hemizygous_row('X', 76900000, 77100000, .02), hemizygous_row('X', 77100001, size - 1000, 1.)]
        options = TimingOptions(sex='XY', **ALL_ON)
        self.assertIn(name, self.events(self.run_engine(rows, focal=[('X', 76000000, 78000000, 'DEL')], options=options)))

    def test_hemizygous_noise_uses_half_the_slope(self):
        # Doubling the noise would reject this value, and removing the halving would widen it enough to reject too.
        name = 'Focal_homdel_X:76000000-78000000'
        size = self.genome.CHROM_DICT['X']
        rows = [r for r in whole_genome(self.genome, male=True) if r[0] != 'X'] + y_rows(1.)
        rows += [hemizygous_row('X', 1000, 76899999, 1.), hemizygous_row('X', 76900000, 77100000, .2, sigma=.15), hemizygous_row('X', 77100001, size - 1000, 1.)]
        options = TimingOptions(sex='XY', **ALL_ON)
        self.assertIn(name, self.events(self.run_engine(rows, focal=[('X', 76000000, 78000000, 'DEL')], options=options)))

    def test_wgd_tumor_times_homdel_before_wgd(self):
        rows = [r for r in whole_genome(self.genome, male=False, baseline=(2, 2)) if r[0] != '9']
        size = self.genome.CHROM_DICT['9']
        rows += [seg_row('9', 1000, 21899999, 2, 2), seg_row('9', 21900000, 22100000, 0, 0, q1=.01, q2=.01), seg_row('9', 22100001, size - 1000, 2, 2)]
        depth = 80
        muts = []
        for c in map(str, range(1, 9)):
            for i in range(10):
                alt = int(depth * (.444 if i % 2 else .222))
                muts.append((c, 2000000 + i * 1000, 1, 'A', depth - alt, alt))
                muts.append((c, 140000000 + i * 1000, 1, 'A', depth - alt, alt))
        engine = self.run_engine(rows, muts=muts, focal=[self.REGION], options=TimingOptions(sex='XX', **ALL_ON))
        self.assertIsNotNone(engine.WGD)
        event = self.events(engine)[self.NAME]
        self.assertEqual(event.timing_status, 'estimated')
        self.assertLess(self.pi_mean(event.pi_dist), self.pi_mean(engine.WGD.pi_dist))


class Reader(TimingCase):
    def test_duplicate_rows_with_missing_columns_collapse(self):
        # NaN in a segment's data breaks Interval equality. Original behavior collapses identical rows.
        rows = whole_genome(self.genome, male=False)
        duplicate = seg_row('22', 1000, 2000, 1, 1, q1=float('nan'), q2=float('nan'))
        duplicate[6] = ''
        rows += [duplicate, list(duplicate)]
        maf_path, seg_path = self.write_inputs('S1', rows, ())
        patient = Patient(indiv_name='P1', genome_build='hg19')
        patient.addSample(maf_path, 'S1', input_type='post-clustering', timepoint_value=0, grid_size=101,
                          _additional_muts=None, seg_file=seg_path, purity=PURITY)
        tree = patient.sample_list[0].CnProfile['22']
        self.assertEqual(len([iv for iv in tree if iv.begin == 1000 and iv.end == 2000]), 1)


class CommandLine(TimingCase):
    def run_cli(self, rows, sex, *flags, focal=None, extra=None, muts=()):
        maf_path, seg_path = self.write_inputs('S1', rows, muts)
        sif = os.path.join(self.tmp, 'input.sif')
        table = {'sample_id': ['S1'], 'maf_fn': [maf_path], 'seg_fn': [seg_path], 'purity': [PURITY], 'timepoint': [0], 'tmb': [1]}
        if sex is not None:
            table['sex'] = [sex]
        for column, value in (extra or {}).items():
            table[column] = [value]
        pd.DataFrame(table).to_csv(sif, sep='\t', index=False)
        command = [sys.executable, os.path.join(REPO, 'PhylogicNDT.py'), 'Timing', '-i', 'P1', '-sif', sif, *flags]
        if focal:
            path = os.path.join(self.tmp, 'focal.txt')
            with open(path, 'w') as handle:
                handle.write(focal)
            command += ['--cn_peaks', path]
        return subprocess.run(command, cwd=self.tmp, capture_output=True, text=True,
                              env=dict(os.environ, PYTHONPATH=REPO, PYTHONDONTWRITEBYTECODE='1'))

    def table(self):
        return pd.read_csv(os.path.join(self.tmp, 'P1.timing.tsv'), sep='\t', low_memory=False)

    def test_flags_reach_the_engine(self):
        rows = whole_genome(self.genome, male=True) + y_rows(.57)
        rows = [r for r in rows if r[0] != '17'] + [seg_row('17', 1000, 22000000, 1, 1, q1=.6), seg_row('17', 25000000, 81000000, 1, 1)]
        done = self.run_cli(rows, 'XY', '--sex_aware', '--continuous_cn', '--clonal_unresolved')
        self.assertEqual(done.returncode, 0, done.stderr[-800:])
        table = self.table().set_index('Event Name', drop=False)
        self.assertEqual(table.loc['loss_Y', 'timing_status'], 'subclonal')
        self.assertAlmostEqual(float(table.loc['loss_17p', 'ccf_hat']), .4)
        self.assertNotIn('loss_Xp', table.index)

    def test_homdel_flags_reach_the_call(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '9']
        size = self.genome.CHROM_DICT['9']
        rows += [seg_row('9', 1000, 21899999, 1, 1), seg_row('9', 21900000, 22100000, 0, 0, q1=.01, q2=.01), seg_row('9', 22100001, size - 1000, 1, 1)]
        artifact = os.path.join(self.tmp, 'artifacts.txt')
        with open(artifact, 'w') as handle:
            handle.write('9\t21800000\t22200000\t5\n')
        focal = '9\t21000000\t23000000\tDEL\t0.01\tregion\n'
        flags = ['--sex_aware', '--continuous_cn', '--clonal_unresolved', '--segment_homdels']
        self.assertEqual(self.run_cli(rows, 'XX', *flags, focal=focal).returncode, 0)
        self.assertIn('Focal_homdel_9:21000000-23000000', set(self.table()['Event Name']))
        self.assertEqual(self.run_cli(rows, 'XX', *flags, '--homdel_artifact_file', artifact, focal=focal).returncode, 0)
        self.assertNotIn('Focal_homdel_9:21000000-23000000', set(self.table()['Event Name']))

    def test_thresholds_reach_the_calls(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '9']
        size = self.genome.CHROM_DICT['9']
        rows += [seg_row('9', 1000, 21899999, 1, 1), seg_row('9', 21900000, 22100000, 0, 0, q1=.1, q2=.1), seg_row('9', 22100001, size - 1000, 1, 1)]
        focal = '9\t21000000\t23000000\tDEL\t0.01\tregion\n'
        flags = ['--sex_aware', '--continuous_cn', '--clonal_unresolved', '--segment_homdels']
        name = 'Focal_homdel_9:21000000-23000000'
        self.assertEqual(self.run_cli(rows, 'XX', *flags, focal=focal).returncode, 0)
        self.assertIn(name, set(self.table()['Event Name']))
        self.assertEqual(self.run_cli(rows, 'XX', *flags, '--homdel_max_copies', '0.1', focal=focal).returncode, 0)
        self.assertNotIn(name, set(self.table()['Event Name']))
        self.assertEqual(self.run_cli(rows, 'XX', *flags, '--homdel_max_length', '1000', focal=focal).returncode, 0)
        self.assertNotIn(name, set(self.table()['Event Name']))
        self.assertEqual(self.run_cli(rows, 'XX', *flags, '--homdel_min_z', '50', focal=focal).returncode, 0)
        self.assertNotIn(name, set(self.table()['Event Name']))
        self.assertEqual(self.run_cli(rows, 'XX', *flags, '--cn_clonal_min_ccf', '0.5', focal=focal).returncode, 0)
        table = self.table().set_index('Event Name')
        self.assertEqual(str(table.loc[name, 'is_clonal']), 'True')
        self.assertEqual(self.run_cli(rows, 'XX', *flags, '--cn_clonal_min_ccf', '1.0', focal=focal).returncode, 0)
        self.assertEqual(self.table().set_index('Event Name').loc[name, 'timing_status'], 'subclonal')

    def test_sif_builder_writes_sex(self):
        command = [sys.executable, os.path.join(REPO, 'create_patient_sif.py'), '--patient_id', 'P1.REP0', '--sample_names', 'S1',
                   '--absolute_mafs', 'a.maf', '--absolute_segtabs', 'a.seg', '--absolute_purities', '.8', '--sex', 'XY',
                   '--outfile', os.path.join(self.tmp, 'out.sif')]
        self.assertEqual(subprocess.run(command, capture_output=True).returncode, 0)
        self.assertEqual(pd.read_csv(os.path.join(self.tmp, 'out.sif'), sep='\t')['sex'].iloc[0], 'XY')
        timing = pd.read_csv(os.path.join(self.tmp, 'out.sif.timing.txt.gz'), sep='\t')
        self.assertEqual(timing['sex'].iloc[0], 'XY')
        self.assertEqual(timing['maf_fn'].iloc[0], 'P1.REP0.mut_ccfs.txt.gz')

    def test_bad_input_stops_the_run(self):
        rows = whole_genome(self.genome, male=True) + y_rows(1.)
        for sex in ('Unknown', None):
            self.assertNotEqual(self.run_cli(rows, sex, '--sex_aware').returncode, 0)
        self.assertNotEqual(self.run_cli(rows, 'XY', '--segment_homdels', '--clonal_unresolved', '--continuous_cn').returncode, 0)
        self.assertNotEqual(self.run_cli(rows, 'XY', '--continuous_cn', '--intersect_cn_trees').returncode, 0)
        self.assertNotEqual(self.run_cli(rows, 'XY', '--sex_aware', '--genome_build', 'hg38').returncode, 0)
        self.assertNotEqual(self.run_cli(rows, 'XY', '--homdel_artifact_file', 'x').returncode, 0)


def doubled_arm_muts(chromosomes=tuple(map(str, range(1, 9))), depth=80):
    """Clonal SNVs at multiplicity 1 and 2 on both arms of 2/2 chromosomes, which time a WGD."""
    muts = []
    for c in chromosomes:
        for i in range(10):
            alt = int(depth * (.444 if i % 2 else .222))
            muts += [(c, 2000000 + i * 1000, 1, 'A', depth - alt, alt), (c, 140000000 + i * 1000, 1, 'A', depth - alt, alt)]
    return muts


class WgdOverride(TimingCase):
    """Fix 2: a doubled status given in the SIF replaces Timing's arm rule for that sample."""
    def partly_doubled(self):
        # 16 of 41 called arms at 2/2: short of the half the arm rule asks for.
        doubled = [str(c) for c in range(1, 9)]
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] not in doubled]
        rows += [seg_row(c, 1000, self.genome.CHROM_DICT[c] - 1000, 2, 2) for c in doubled]
        return rows

    def options(self, override=None):
        return TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True, wgd_override=override)

    def test_true_calls_a_wgd_the_rule_misses(self):
        rows = self.partly_doubled()
        self.assertIsNone(self.run_engine(rows, muts=doubled_arm_muts(), options=self.options()).WGD)
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=self.options({'S1': True}))
        self.assertIsNotNone(engine.WGD)
        self.assertEqual(engine.WGD.timing_status, 'estimated')
        # Arms at 1/1 now read as losses after the doubling, the 2/2 arms as no event.
        self.assertIn('loss_9p', self.events(engine))
        self.assertNotIn('gain_1p', self.events(engine))

    def test_false_removes_a_wgd_the_rule_calls(self):
        rows = whole_genome(self.genome, male=False, baseline=(2, 2))
        self.assertIsNotNone(self.run_engine(rows, muts=doubled_arm_muts(), options=self.options()).WGD)
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=self.options({'S1': False}))
        self.assertIsNone(engine.WGD)
        self.assertIn('gain_1p', self.events(engine))

    def test_other_samples_keep_the_rule(self):
        # Only S2 is given. It doubles; S1 keeps its rule, which says no here, so the patient has no WGD.
        muts = doubled_arm_muts()
        engine = self.run_samples([('S1', self.partly_doubled(), muts, PURITY), ('S2', self.partly_doubled(), muts, PURITY)],
                                  options=self.options({'S2': True}))
        by_name = {sample.sample_name: sample for sample in engine.sample_list}
        self.assertIsNone(by_name['S1'].WGD)
        self.assertIsNotNone(by_name['S2'].WGD)
        self.assertIsNone(engine.WGD)
        # Here S1's rule says yes, and with S2 given as doubled the patient has one.
        doubled = whole_genome(self.genome, male=False, baseline=(2, 2))
        engine = self.run_samples([('S1', doubled, muts, PURITY), ('S2', self.partly_doubled(), muts, PURITY)],
                                  options=self.options({'S2': True}))
        self.assertTrue(all(sample.WGD is not None for sample in engine.sample_list))
        self.assertIsNotNone(engine.WGD)

    def test_sif_column_reaches_the_engine(self):
        cli = CommandLine('run_cli')
        cli.tmp, cli.genome = self.tmp, self.genome
        rows = self.partly_doubled()
        flags = ['--sex_aware', '--continuous_cn', '--clonal_unresolved']
        for value, expected in (('True', True), ('', False), ('False', False)):
            done = cli.run_cli(rows, 'XX', *flags, extra={'doubled': value}, muts=doubled_arm_muts())
            self.assertEqual(done.returncode, 0, done.stderr[-800:])
            self.assertEqual('WGD' in set(cli.table()['Event Name']), expected, msg=value)
        self.assertNotEqual(cli.run_cli(rows, 'XX', *flags, extra={'doubled': 'maybe'}).returncode, 0)


class UncappedExcess(TimingCase):
    """Fix 3: a gain far above the next copy is clonal however noisy, because its excess is tested before the cap."""
    def gain(self, q2, sigma, continuous=True):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, 1, int(round(q2)), q1=1., q2=q2, sigma=sigma), seg_row('17', 25000000, 81000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=continuous, clonal_unresolved=True)
        return self.events(self.run_engine(rows, options=options))['gain_17p']

    def test_high_noisy_gain_is_clonal(self):
        # As NBE7 17q at 6.1 copies: CCF capped at 1 with sd 0.16 lies within the margin of the 0.8 cutoff.
        event = self.gain(6.1, .09)
        self.assertGreater(event.ccf_sd, .15)
        self.assertTrue(event.is_clonal)
        self.assertEqual(event.timing_status, 'clonal_unresolved')
        self.assertEqual(event.ccf_hat, 1.)

    def test_one_copy_gain_is_judged_as_before(self):
        self.assertEqual(self.gain(1.9, .09).timing_status, 'ambiguous_clonality')
        self.assertEqual(self.gain(1.4, .01).timing_status, 'subclonal')

    def test_loss_keeps_the_capped_test(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, 0, 1, q1=.12, sigma=.09), seg_row('17', 25000000, 81000000, 1, 1)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True)
        self.assertEqual(self.events(self.run_engine(rows, options=options))['loss_17p'].timing_status, 'ambiguous_clonality')


class FocalGainShape(FocalCase):
    """Fix 4: a focal gain needs two real flanks and a short segment, and high level is relative to ploidy."""
    def test_gain_against_the_centromere_has_one_flank(self):
        # A short gained segment starting at the 11q centromere edge, measured against its distal side only.
        size = self.genome.CHROM_DICT['11']
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '11']
        rows += [seg_row('11', 1000, 51000000, 1, 1), seg_row('11', 54000000, 56000000, 1, 2), seg_row('11', 56000001, size - 1000, 1, 1)]
        self.assertNotIn('Focal_gain_11:55000000-55100000', self.focal(rows, [('11', 55000000, 55100000, 'AMP')]))

    def test_gain_against_the_chromosome_end_has_one_flank(self):
        rows = self.with_segments('5', seg_row('5', 10000, 1300000, 1, 2))
        rows = [r for r in rows if not (r[0] == '5' and r[2] == 9999)]
        self.assertNotIn('Focal_gain_5:1253262-1295184', self.focal(rows, [('5', 1253262, 1295184, 'AMP')]))

    def test_focal_gain_near_the_chromosome_end_keeps_a_short_flank(self):
        # TERT lies 1.25 Mb from the 5p end: a flank of about 1 Mb is real.
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '5']
        rows += [seg_row('5', 10000, 1099999, 1, 1), seg_row('5', 1100000, 1300000, 1, 3),
                 seg_row('5', 1300001, self.genome.CHROM_DICT['5'] - 1000, 1, 1)]
        self.assertIn('Focal_gain_5:1253262-1295184', self.focal(rows, [('5', 1253262, 1295184, 'AMP')]))

    def test_amplicon_at_an_arm_end_fills_its_own_flank(self):
        # A high-level amplification needs no second flank: here it starts at the 5p end.
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '5']
        rows += [seg_row('5', 10000, 1300000, 1, 1, q1=1., q2=40.), seg_row('5', 1300001, self.genome.CHROM_DICT['5'] - 1000, 1, 1)]
        event = self.focal(rows, [('5', 1253262, 1295184, 'AMP')])['Focal_gain_5:1253262-1295184']
        self.assertTrue(event.high_level)
        # Against the centromere too.
        size = self.genome.CHROM_DICT['11']
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '11']
        rows += [seg_row('11', 1000, 51000000, 1, 1), seg_row('11', 54000000, 56000000, 1, 1, q1=1., q2=40.),
                 seg_row('11', 56000001, size - 1000, 1, 1)]
        self.assertIn('Focal_gain_11:55000000-55100000', self.focal(rows, [('11', 55000000, 55100000, 'AMP')]))

    def test_deletions_keep_the_four_copy_noise_exemption(self):
        # A noisy loss of 5 copies at ploidy about 4: the exemption from the noise test starts at 4 copies for losses,
        # not at twice the ploidy, which applies to gains only.
        size = self.genome.CHROM_DICT['9']
        rows = [r for r in whole_genome(self.genome, male=False, baseline=(2, 2)) if r[0] != '9']
        rows += [seg_row('9', 1000, 8999999, 3, 3), seg_row('9', 9000000, 9500000, 0, 1, sigma=1.), seg_row('9', 9500001, size - 1000, 3, 3)]
        self.assertIn('Focal_loss_9:8000000-11000000', self.focal(rows, [self.DEL]))

    def test_long_gain_is_not_focal(self):
        self.assertNotIn('Focal_gain_2:16000000-16200000', self.focal(self.with_segments('2', seg_row('2', 14000000, 17500000, 1, 2)), [self.AMP]))
        self.assertIn('Focal_gain_2:16000000-16200000', self.focal(self.with_segments('2', seg_row('2', 14500000, 17000000, 1, 2)), [self.AMP]))

    def test_long_amplicon_segment_is_not_focal_either(self):
        self.assertNotIn('Focal_gain_2:16000000-16200000',
                         self.focal(self.with_segments('2', seg_row('2', 14000000, 17500000, 1, 1, q1=1., q2=60.)), [self.AMP]))

    def plateau(self, baseline, step):
        """Nine 0.5 Mb segments all `step` copies above the arm, so the flanks are as high as the window."""
        a2 = baseline[1] + step
        middle = [seg_row('2', 14000000 + i * 500000, 14000000 + (i + 1) * 500000 - 1, *baseline, q1=float(baseline[0]), q2=float(a2))
                  for i in range(9)]
        return self.with_segments('2', *middle, baseline=baseline)

    def test_high_level_is_relative_to_ploidy(self):
        # Five copies above the arm: high level at ploidy 2 (from 4 copies), not at ploidy 4 (from 8).
        self.assertIn('Focal_gain_2:16000000-16200000', self.focal(self.plateau((1, 1), 5), [self.AMP]))
        self.assertNotIn('Focal_gain_2:16000000-16200000', self.focal(self.plateau((2, 2), 5), [self.AMP]))
        self.assertIn('Focal_gain_2:16000000-16200000', self.focal(self.plateau((2, 2), 9), [self.AMP]))

    def test_cli_parameters_reach_the_call(self):
        cli = CommandLine('run_cli')
        cli.tmp, cli.genome = self.tmp, self.genome
        rows = self.with_segments('2', seg_row('2', 14500000, 17000000, 1, 2))
        focal = '2\t16000000\t16200000\tAMP\t0.01\tregion\n'
        flags = ['--sex_aware', '--continuous_cn', '--clonal_unresolved']
        self.assertEqual(cli.run_cli(rows, 'XX', *flags, focal=focal).returncode, 0)
        self.assertIn('Focal_gain_2:16000000-16200000', set(cli.table()['Event Name']))
        self.assertEqual(cli.run_cli(rows, 'XX', *flags, '--focal_gain_max_length', '2000000', focal=focal).returncode, 0)
        self.assertNotIn('Focal_gain_2:16000000-16200000', set(cli.table()['Event Name']))
        rows = self.plateau((1, 1), 5)
        self.assertEqual(cli.run_cli(rows, 'XX', *flags, '--focal_amp_ploidy', '3', focal=focal).returncode, 0)
        self.assertNotIn('Focal_gain_2:16000000-16200000', set(cli.table()['Event Name']))


class AmpliconsAndWgd(FocalCase):
    """Fix 5: in a doubled tumor, no focal gain is ordered against the WGD by its copy number, whatever its level."""
    def doubled(self, *middle, muts=(), baseline=(2, 2), focal=None):
        return self.run_engine(self.with_segments('2', *middle, baseline=baseline), muts=doubled_arm_muts() + list(muts),
                               focal=[focal or self.AMP], options=TimingOptions(sex='XX', **ALL_ON))

    def test_mutations_in_a_focal_gain_are_unresolved_too(self):
        # The window names the event; the amplified segment, 2 Mb here, reaches well beyond it.
        window = ('2', 16080686, 16087129, 'AMP')
        in_window = [('2', 16083000 + i * 1000, 1, 'A', 44, 36) for i in range(2)]
        in_segment = [('2', 15500000 + i * 1000, 1, 'A', 44, 36) for i in range(2)]
        flank = [('2', 17500000 + i * 1000, 1, 'A', 44, 36) for i in range(2)]
        for middle in (seg_row('2', 15000000, 17000000, 2, 2, q1=2., q2=99.), seg_row('2', 15000000, 17000000, 2, 4)):
            engine = self.doubled(middle, muts=in_window + in_segment + flank, focal=window)
            event = self.events(engine)['Focal_gain_2:16080686-16087129']
            self.assertEqual((event.start, event.end), (16080686, 16087129))
            self.assertEqual(event.segments, ((15000000, 17000000),))
            muts = {m.pos: m for m in engine.mutations.values() if m.chrN == '2'}
            for _, pos, *_ in in_window + in_segment:
                self.assertIsNone(muts[pos].p_before_wgd)
                self.assertEqual((muts[pos].timing_status, muts[pos].timing_reason), ('clonal_unresolved', 'amplicon_vs_wgd'))
                self.assertTrue(np.allclose(muts[pos].pi_dist, 1. / 101))
            for _, pos, *_ in flank:
                self.assertIsNotNone(muts[pos].p_before_wgd)
                self.assertIsNone(getattr(muts[pos], 'timing_reason', None))
        os.chdir(self.tmp)
        self.addCleanup(os.chdir, self.cwd)
        PhylogicOutput().write_timing_tsv(engine)
        table = pd.read_csv('P1.timing.tsv', sep='\t', low_memory=False, dtype={'Chromosome': str})
        rows = table[(table['Chromosome'] == '2') & table['Start_position'].between(15000000, 16999999)]
        self.assertEqual(len(rows), 4)
        self.assertEqual(set(rows['timing_reason']), {'amplicon_vs_wgd'})
        self.assertTrue(rows['p_before_wgd'].isna().all())

    def test_the_span_covers_an_amplicon_of_several_segments(self):
        # The window lies in the higher segment; the lower one, still far above the arm, is part of the same amplicon.
        window = ('2', 16300000, 16310000, 'AMP')
        snvs = [('2', 15500000, 1, 'A', 44, 36), ('2', 16300500, 1, 'A', 44, 36), ('2', 17800000, 1, 'A', 44, 36)]
        engine = self.doubled(seg_row('2', 15000000, 16000000, 2, 2, q1=2., q2=36.), seg_row('2', 16000001, 17000000, 2, 2, q1=2., q2=38.),
                              muts=snvs, focal=window)
        self.assertEqual(self.events(engine)['Focal_gain_2:16300000-16310000'].segments, ((15000000, 17000000),))
        muts = {m.pos: m for m in engine.mutations.values() if m.chrN == '2'}
        self.assertEqual([muts[pos].p_before_wgd is None for pos in (15500000, 16300500, 17800000)], [True, True, False])

    def test_the_span_runs_through_a_dip_that_is_still_high_level(self):
        # Peak at 62 copies on a doubled arm of 4: the midpoint is 33, the high-level line 4 + 2 x ploidy, about 12.
        # The dip at 22 copies is below the midpoint but high-level itself, so the amplicon goes on beyond it.
        window = ('2', 15200000, 15210000, 'AMP')
        snvs = [('2', 16200000, 1, 'A', 44, 36), ('2', 17500000, 1, 'A', 44, 36)]
        engine = self.doubled(seg_row('2', 15000000, 15500000, 2, 2, q1=2., q2=60.), seg_row('2', 15500001, 15800000, 2, 2, q1=2., q2=20.),
                              seg_row('2', 15800001, 16500000, 2, 2, q1=2., q2=50.), muts=snvs, focal=window)
        self.assertEqual(self.events(engine)['Focal_gain_2:15200000-15210000'].segments, ((15000000, 16500000),))
        muts = {m.pos: m for m in engine.mutations.values() if m.chrN == '2'}
        self.assertEqual([muts[pos].p_before_wgd is None for pos in (16200000, 17500000)], [True, False])

    def test_a_low_level_gain_keeps_the_midpoint(self):
        # A gain of 2 copies with a neighbour 0.8 above the arm: below the midpoint, so not part of the gain.
        window = ('2', 15200000, 15210000, 'AMP')
        engine = self.doubled(seg_row('2', 15000000, 15500000, 2, 4), seg_row('2', 15500001, 16500000, 2, 2, q1=2., q2=2.8), focal=window)
        self.assertEqual(self.events(engine)['Focal_gain_2:15200000-15210000'].segments, ((15000000, 15500000),))

    def test_segments_of_every_carrier_count(self):
        # Two samples measure the amplicon on different segments: an SNV inside either one is not timed.
        def rows(start, end):
            return self.with_segments('2', seg_row('2', start, end, 2, 2, q1=2., q2=99.), baseline=(2, 2))
        window = ('2', 16080686, 16087129, 'AMP')
        snvs = [('2', 15200000, 1, 'A', 44, 36), ('2', 16800000, 1, 'A', 44, 36), ('2', 17500000, 1, 'A', 44, 36)]
        muts = doubled_arm_muts() + snvs
        engine = self.run_samples([('S1', rows(15000000, 16500000), muts, PURITY), ('S2', rows(15600000, 17000000), muts, PURITY)],
                                  focal=[window], options=TimingOptions(sex='XX', **ALL_ON))
        self.assertEqual(self.events(engine)['Focal_gain_2:16080686-16087129'].segments, ((15000000, 16500000), (15600000, 17000000)))
        muts = {m.pos: m for m in engine.mutations.values() if m.chrN == '2'}
        self.assertEqual([muts[pos].p_before_wgd is None for pos in (15200000, 16800000, 17500000)], [True, True, False])

    def test_a_carrier_outside_the_consensus_state_counts(self):
        # The carriers reach different states, so only S1 holds the consensus. S2's longer amplicon still counts.
        def rows(end, copies):
            return self.with_segments('2', seg_row('2', 15000000, end, 2, 2, q1=2., q2=copies), baseline=(2, 2))
        window = ('2', 15200000, 15210000, 'AMP')
        muts = doubled_arm_muts() + [('2', 16500000, 1, 'A', 44, 36)]
        engine = self.run_samples([('S1', rows(16000000, 60.), muts, PURITY), ('S2', rows(17000000, 30.), muts, PURITY)],
                                  focal=[window], options=TimingOptions(sex='XX', **ALL_ON))
        event = self.events(engine)['Focal_gain_2:15200000-15210000']
        self.assertEqual(event.carrier_sample_idx, (0,))
        self.assertEqual(event.segments, ((15000000, 16000000), (15000000, 17000000)))
        mut = next(m for m in engine.mutations.values() if m.chrN == '2' and m.pos == 16500000)
        self.assertIsNone(mut.p_before_wgd)
        self.assertEqual(mut.timing_reason, 'amplicon_vs_wgd')

    def test_mutations_in_a_focal_gain_without_wgd_are_untouched(self):
        inside = [('2', 16050000 + i * 1000, 1, 'A', 44, 36) for i in range(3)]
        engine = self.doubled(seg_row('2', 15900000, 16300000, 1, 1, q1=1., q2=99.), muts=inside, baseline=(1, 1))
        self.assertIsNone(engine.WGD)
        muts = [m for m in engine.mutations.values() if m.chrN == '2' and 16000000 <= m.pos < 16200000]
        self.assertTrue(muts)
        self.assertTrue(all(getattr(m, 'timing_reason', None) is None for m in muts))

    def test_amplicon_is_unresolved_against_the_wgd(self):
        engine = self.doubled(seg_row('2', 15900000, 16300000, 2, 2, q1=2., q2=99.))
        self.assertIsNotNone(engine.WGD)
        event = self.events(engine)['Focal_gain_2:16000000-16200000']
        self.assertTrue(event.high_level)
        self.assertEqual(event.timing_status, 'clonal_unresolved')
        self.assertEqual(event.timing_reason, 'amplicon_vs_wgd')
        self.assertIsNone(event.pi_dist)

    def test_low_level_focal_gain_is_unresolved_too(self):
        # Two copies above a doubled arm: the parity rule would place it before the WGD.
        engine = self.doubled(seg_row('2', 15900000, 16300000, 2, 4))
        event = self.events(engine)['Focal_gain_2:16000000-16200000']
        self.assertFalse(event.high_level)
        self.assertEqual(event.timing_status, 'clonal_unresolved')
        self.assertEqual(event.timing_reason, 'amplicon_vs_wgd')
        self.assertIsNone(event.pi_dist)
        self.assertIsNone(event.wgd_side)

    def test_focal_losses_keep_the_rule(self):
        size = self.genome.CHROM_DICT['9']
        rows = [r for r in whole_genome(self.genome, male=False, baseline=(2, 2)) if r[0] != '9']
        rows += [seg_row('9', 1000, 8999999, 2, 2), seg_row('9', 9000000, 9500000, 1, 2), seg_row('9', 9500001, size - 1000, 2, 2)]
        engine = self.run_engine(rows, muts=doubled_arm_muts(), focal=[self.DEL], options=TimingOptions(sex='XX', **ALL_ON))
        event = self.events(engine)['Focal_loss_9:8000000-11000000']
        self.assertEqual((event.timing_status, event.timing_source, event.wgd_side), ('estimated', 'wgd_rule', 'after'))


class SameSideOfWgd(TimingCase):
    """Fix 6: two events the WGD rule places on one side of the WGD cannot be ordered against each other."""
    @staticmethod
    def dists():
        grid = np.linspace(0, 1, 101)
        wgd = np.exp(-((grid - .4) ** 2) / .002)
        wgd /= wgd.sum()
        before = 1 - np.cumsum(wgd)
        after = np.cumsum(wgd)
        return wgd, before / before.sum(), after / after.sum()

    def test_rule(self):
        wgd, before, after = self.dists()
        comps = spt.compare_event_list([
            spt.TimedEvent('WGD', wgd, 'estimated'),
            spt.TimedEvent('loss_Y', before, 'estimated', 'before'),
            spt.TimedEvent('loss_11q', before, 'estimated', 'before'),
            spt.TimedEvent('loss_3p', after, 'estimated', 'after'),
            spt.TimedEvent('gain_7p', after, 'estimated', 'after'),
        ])
        self.assertEqual(comps[('loss_Y', 'loss_11q')], (0.0, 0.0, 1.0))
        self.assertEqual(comps[('loss_3p', 'gain_7p')], (0.0, 0.0, 1.0))
        self.assertGreater(comps[('loss_Y', 'loss_3p')][0], .9)
        self.assertGreater(comps[('WGD', 'loss_Y')][1], .8)
        # Without sides, identical distributions still score as a game, as before.
        self.assertGreater(spt.compare_event_list([spt.TimedEvent('a', before, 'estimated'),
                                                   spt.TimedEvent('b', before, 'estimated')])[('a', 'b')][0], .2)

    def test_table_round_trip(self):
        rows = [r for r in whole_genome(self.genome, male=True, baseline=(2, 2)) if r[0] != '17'] + y_rows(0.)
        rows += [seg_row('17', 1000, 22000000, 0, 2), seg_row('17', 25000000, 81000000, 2, 2)]
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=TimingOptions(sex='XY', **ALL_ON))
        events = self.events(engine)
        self.assertEqual((events['loss_Y'].wgd_side, events['loss_17p'].wgd_side), ('before', 'before'))
        os.chdir(self.tmp)
        self.addCleanup(os.chdir, self.cwd)
        PhylogicOutput().write_timing_tsv(engine)
        table = pd.read_csv('P1.timing.tsv', sep='\t', low_memory=False).set_index('Event Name')
        self.assertEqual(table.loc['loss_Y', 'wgd_side'], 'before')
        self.assertTrue(pd.isna(table.loc['WGD', 'wgd_side']))
        comps = spt.compare_event_list(spt.read_timed_events('P1.timing.tsv', compare_all=True))
        pair = ('loss_Y', 'loss_17p') if ('loss_Y', 'loss_17p') in comps else ('loss_17p', 'loss_Y')
        self.assertEqual(comps[pair], (0.0, 0.0, 1.0))
        self.assertGreater(max(comps[('WGD', 'loss_Y')][:2]), .8)

    def test_engine_comparisons_use_the_sides(self):
        rows = [r for r in whole_genome(self.genome, male=True, baseline=(2, 2)) if r[0] != '17'] + y_rows(0.)
        rows += [seg_row('17', 1000, 22000000, 0, 2), seg_row('17', 25000000, 81000000, 2, 2)]
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=TimingOptions(sex='XY', **ALL_ON))
        comps = spt.compare_events(engine)
        pair = ('loss_Y', 'loss_17p') if ('loss_Y', 'loss_17p') in comps else ('loss_17p', 'loss_Y')
        self.assertEqual(comps[pair], (0.0, 0.0, 1.0))


class WgdAnchoredPairs(TimingCase):
    """Extension of fix 6: SNVs timed against the WGD carry their probability of preceding it, as rule events do."""
    def doubled(self):
        rows = [r for r in whole_genome(self.genome, male=True, baseline=(2, 2)) if r[0] != '17'] + y_rows(0.)
        rows += [seg_row('17', 1000, 22000000, 0, 2), seg_row('17', 25000000, 81000000, 1, 2)]
        return self.run_engine(rows, muts=doubled_arm_muts(), options=TimingOptions(sex='XY', **ALL_ON))

    def test_rule(self):
        grid = np.linspace(0, 1, 101)
        wgd = np.exp(-((grid - .4) ** 2) / .002)
        wgd /= wgd.sum()
        before, after = 1 - np.cumsum(wgd), np.cumsum(wgd)
        before, after = before / before.sum(), after / after.sum()
        snv = before * .9 + after * .1
        comps = spt.compare_event_list([
            spt.TimedEvent('WGD', wgd, 'estimated'),
            spt.TimedEvent('loss_Y', before, 'estimated', 'before', 1.),
            spt.TimedEvent('loss_17q', after, 'estimated', 'after', 0.),
            spt.TimedEvent('snv', snv, 'estimated', None, .9),
            spt.TimedEvent('old_table', before, 'estimated', 'before'),
            spt.TimedEvent('gain_timed', after, 'estimated'),
        ])
        self.assertEqual(comps[('loss_Y', 'loss_17q')], (1.0, 0.0, 0.0))
        self.assertEqual(comps[('loss_Y', 'old_table')], (0.0, 0.0, 1.0))
        for got, expected in zip(comps[('loss_Y', 'snv')], (.1, 0., .9)):
            self.assertAlmostEqual(got, expected)
        for got, expected in zip(comps[('loss_17q', 'snv')], (0., .9, .1)):
            self.assertAlmostEqual(got, expected)
        # Against the WGD itself, and against an event without an anchor, the distributions decide as before.
        self.assertGreater(comps[('WGD', 'loss_Y')][1], .8)
        self.assertEqual(comps[('WGD', 'snv')], spt.compare_event_list([spt.TimedEvent('WGD', wgd), spt.TimedEvent('snv', snv)])[('WGD', 'snv')])
        self.assertEqual(comps[('snv', 'gain_timed')], spt.compare_event_list([spt.TimedEvent('snv', snv), spt.TimedEvent('gain_timed', after)])[('snv', 'gain_timed')])

    def test_uniform_distribution_stays_unknown(self):
        # A multiplicity-1 SNV on 1/2 has p_before = p_after: its distribution is uniform, its q only the WGD's weight.
        grid = np.linspace(0, 1, 101)
        wgd = np.exp(-((grid - .4) ** 2) / .002)
        wgd /= wgd.sum()
        before = 1 - np.cumsum(wgd)
        comps = spt.compare_event_list([spt.TimedEvent('loss_Y', before / before.sum(), 'estimated', 'before', 1.),
                                        spt.TimedEvent('snv', np.ones(101) / 101, 'estimated', None, .41)])
        self.assertEqual(comps[('loss_Y', 'snv')], (0.0, 0.0, 1.0))

    def test_multiplicity_one_on_one_two_is_unknown_in_the_engine(self):
        rows = [r for r in whole_genome(self.genome, male=True, baseline=(2, 2)) if r[0] != '17'] + y_rows(0.)
        rows += [seg_row('17', 1000, 22000000, 0, 2), seg_row('17', 25000000, 81000000, 1, 2)]
        depth = 100
        alt = int(depth * PURITY / (3 * PURITY + 2 * (1 - PURITY)))
        muts = doubled_arm_muts() + [('17', 40000000 + i * 1000, 1, 'A', depth - alt, alt) for i in range(3)]
        engine = self.run_engine(rows, muts=muts, options=TimingOptions(sex='XY', **ALL_ON))
        snvs = [m for m in engine.mutations.values() if m.chrN == '17' and m.pos >= 40000000]
        self.assertTrue(all(m.p_before_wgd is not None for m in snvs))
        comps = spt.compare_events(engine, drivers=[snvs[0].gene])
        pair = next(key for key in comps if 'loss_Y' in key and snvs[0].event_name in key)
        self.assertEqual(comps[pair], (0.0, 0.0, 1.0))

    def test_snvs_get_their_probability_from_multiplicity(self):
        engine = self.doubled()
        timed = [mut for mut in engine.mutations.values() if mut.pi_dist is not None and mut.is_clonal]
        self.assertTrue(timed)
        # On a 2/2 arm, multiplicity 2 (alt fraction 0.444) means before the doubling, multiplicity 1 after it.
        early = [m.p_before_wgd for m in timed if m.alt_cnt[0] / (m.alt_cnt[0] + m.ref_cnt[0]) > .4]
        late = [m.p_before_wgd for m in timed if m.alt_cnt[0] / (m.alt_cnt[0] + m.ref_cnt[0]) < .3]
        self.assertTrue(early and late)
        self.assertGreater(min(early), .9)
        self.assertLess(max(late), .1)
        events = self.events(engine)
        self.assertEqual((events['loss_Y'].p_before_wgd, events['loss_17p'].p_before_wgd), (1., 1.))

    @staticmethod
    def timed_mut(n_samples, alt, ref=80):
        """A clonal SNV on 2/2 in n identical samples, timed against an early WGD, whose two sides weigh unequally."""
        class Sample(object):
            def __init__(self):
                self.purity, self.options, self.sample_name = PURITY, TimingOptions(), 'S'

        class Gain(object):
            grid = np.linspace(0, 1, 101)
            pi_dist = np.exp(-((grid - .3) ** 2) / .005)
            pi_dist /= pi_dist.sum()

        mut = TimingMut([Sample() for _ in range(n_samples)], 'G', '1', 100, 'T', 'A', np.full(n_samples, alt),
                        np.full(n_samples, ref - alt), np.full(n_samples, 2.), np.full(n_samples, 2.),
                        np.tile(clonal_grid(), (n_samples, 1)))
        mut.get_pi_dist(Gain())
        cdf = np.cumsum(Gain.pi_dist)
        return mut, np.sum(1 - cdf), np.sum(cdf)

    def test_q_is_the_before_weight_of_its_distribution(self):
        # The distribution is q (1 - cdf) / A + (1 - q) cdf / B: q is read back from it by least squares.
        grid = np.linspace(0, 1, 101)
        wgd = np.exp(-((grid - .3) ** 2) / .005)
        cdf = np.cumsum(wgd / wgd.sum())
        for alt in (12, 20, 26, 30, 36):
            mut, a, b = self.timed_mut(1, alt)
            u, v = (1 - cdf) / a, cdf / b
            weight = float(np.dot(u - v, mut.pi_dist - v) / np.dot(u - v, u - v))
            self.assertAlmostEqual(mut.p_before_gain, weight, places=6, msg=str(alt))

    def test_samples_combine_as_their_distributions_do(self):
        one, a, b = self.timed_mut(1, 26)
        two = self.timed_mut(2, 26)[0]
        q = one.p_before_gain
        ratio = q / (1 - q) * b / a  # p_before / p_after of one sample
        self.assertTrue(.05 < q < .95)
        self.assertAlmostEqual(two.p_before_gain, a * ratio ** 2 / (a * ratio ** 2 + b))

    def test_table_round_trip(self):
        engine = self.doubled()
        os.chdir(self.tmp)
        self.addCleanup(os.chdir, self.cwd)
        PhylogicOutput().write_timing_tsv(engine)
        table = pd.read_csv('P1.timing.tsv', sep='\t', low_memory=False)
        snvs = table[table['Start_position'].notna()]
        timed = snvs[snvs['timing_status'] == 'estimated']
        self.assertTrue(timed['p_before_wgd'].notna().all())
        self.assertTrue(snvs[snvs['timing_status'] != 'estimated']['p_before_wgd'].isna().all())
        self.assertTrue(pd.isna(table.set_index('Event Name').loc['WGD', 'p_before_wgd']))
        events = {e.event_name: e for e in spt.read_timed_events('P1.timing.tsv', compare_all=True)}
        snv = next(e for name, e in events.items() if not name.startswith(('gain', 'loss', 'WGD')) and e.p_before_wgd is not None)
        self.assertEqual(events['loss_Y'].p_before_wgd, 1.)
        comps = spt.compare_event_list([events['loss_Y'], snv])
        q = snv.p_before_wgd
        for got, expected in zip(comps[('loss_Y', snv.event_name)], (1. - q, 0., q)):
            self.assertAlmostEqual(got, expected)

    def test_the_winning_row_keeps_its_probability(self):
        grid = np.linspace(0, 1, 101)
        def row(center):
            pi = np.exp(-((grid - center) ** 2) / .002)
            return list(pi / pi.sum())
        columns = ['Event Name'] + [f'pi_{float(x) / 100}' for x in range(101)] + ['timing_status', 'wgd_side', 'p_before_wgd']
        pd.DataFrame([['loss_17p'] + row(.8) + ['estimated', 'after', 0.], ['loss_17p'] + row(.2) + ['estimated', 'before', 1.]],
                     columns=columns).to_csv(os.path.join(self.tmp, 'dup.tsv'), sep='\t', index=False)
        event = spt.read_timed_events(os.path.join(self.tmp, 'dup.tsv'), compare_all=True)[0]
        self.assertEqual((event.wgd_side, event.p_before_wgd), ('before', 1.))

    def test_without_clonal_unresolved_nothing_is_recorded(self):
        rows = whole_genome(self.genome, male=False, baseline=(2, 2))
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=TimingOptions(sex='XX', sex_aware=True, continuous_cn=True))
        self.assertTrue(any(mut.pi_dist is not None for mut in engine.mutations.values()))
        self.assertTrue(all(mut.p_before_wgd is None for mut in engine.mutations.values()))


class LoyOdds(TimingCase):
    """Fix 7: a doubled tumor without Y is placed before its WGD only where Tau's odds reach the minimum."""
    def loss_y(self, y_cn, odds, minimum=3.):
        rows = whole_genome(self.genome, male=True, baseline=(2, 2)) + y_rows(y_cn)
        options = TimingOptions(sex='XY', loy_wgd_min_odds=minimum, loy_odds=odds, **ALL_ON)
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=options)
        self.assertIsNotNone(engine.WGD)
        return engine, self.events(engine)['loss_Y']

    def test_strong_odds_keep_the_rule(self):
        engine, event = self.loss_y(0., {'S1': 3.06})
        self.assertEqual(event.timing_status, 'estimated')
        self.assertLess(self.pi_mean(event.pi_dist), self.pi_mean(engine.WGD.pi_dist))

    def test_weak_or_missing_odds_leave_it_unresolved(self):
        for odds in ({'S1': .73}, {}):
            _, event = self.loss_y(0., odds)
            self.assertEqual(event.timing_status, 'clonal_unresolved', msg=str(odds))
            self.assertEqual(event.timing_reason, 'loy_wgd_odds_below_min')
            self.assertIsNone(event.pi_dist)

    def test_one_copy_left_is_unaffected(self):
        engine, event = self.loss_y(1., {})
        self.assertEqual(event.timing_status, 'estimated')
        self.assertGreater(self.pi_mean(event.pi_dist), self.pi_mean(engine.WGD.pi_dist))

    def test_without_the_option_nothing_changes(self):
        rows = whole_genome(self.genome, male=True, baseline=(2, 2)) + y_rows(0.)
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=TimingOptions(sex='XY', loy_odds={'S1': .1}, **ALL_ON))
        self.assertEqual(self.events(engine)['loss_Y'].timing_status, 'estimated')

    def test_cli(self):
        cli = CommandLine('run_cli')
        cli.tmp, cli.genome = self.tmp, self.genome
        rows = whole_genome(self.genome, male=True, baseline=(2, 2)) + y_rows(0.)
        flags = ['--sex_aware', '--continuous_cn', '--clonal_unresolved', '--loy_wgd_min_odds', '3']
        for odds, expected in (('0.73', 'clonal_unresolved'), ('', 'clonal_unresolved'), ('48.2', 'estimated')):
            done = cli.run_cli(rows, 'XY', *flags, extra={'odds_pre_wgd': odds}, muts=doubled_arm_muts())
            self.assertEqual(done.returncode, 0, done.stderr[-800:])
            self.assertEqual(cli.table().set_index('Event Name').loc['loss_Y', 'timing_status'], expected, msg=odds)
        self.assertNotEqual(cli.run_cli(rows, 'XY', *flags, muts=doubled_arm_muts()).returncode, 0)
        self.assertNotEqual(cli.run_cli(rows, 'XY', '--loy_wgd_min_odds', '3', extra={'odds_pre_wgd': '5'}).returncode, 0)


class Seed(TimingCase):
    """Fix 8: --seed makes Timing's own simulation of detection power reproducible."""
    def gained(self, seed):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, 1, 2), seg_row('17', 25000000, 81000000, 1, 1)]
        depth = 100
        alt = int(depth * 2 * PURITY / (3 * PURITY + 2 * (1 - PURITY)))
        muts = [('17', 5000000 + i * 1000, 1, 'A', depth - alt, alt) for i in range(12)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True, seed=seed)
        return self.events(self.run_engine(rows, muts=muts, options=options))['gain_17p']

    def two_gains(self, seed):
        """An engine with gains on 17p and 5p, each timed by its own clonal mutations."""
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] not in ('17', '5')]
        rows += [seg_row('17', 1000, 22000000, 1, 2), seg_row('17', 25000000, 81000000, 1, 1),
                 seg_row('5', 1000, 46000000, 1, 2), seg_row('5', 50000000, 180000000, 1, 1)]
        depth = 100
        def alt(m):
            return int(depth * m * PURITY / (3 * PURITY + 2 * (1 - PURITY)))
        muts = [('17', 5000000 + i * 1000, 1, 'A', depth - alt(2), alt(2)) for i in range(12)]
        muts += [('5', 5000000 + i * 1000, 1, 'A', depth - alt(1 + i % 2), alt(1 + i % 2)) for i in range(12)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True, seed=seed)
        events = self.events(self.run_engine(rows, muts=muts, options=options))
        return events['gain_17p'], events['gain_5p']

    @staticmethod
    def time_in_order(first, second):
        """Time two events one after the other, from the same global random state, and return their distributions."""
        np.random.seed(99)
        random.seed(99)
        first.get_pi_dist_for_gain()
        second.get_pi_dist_for_gain()
        return first.pi_dist.copy(), second.pi_dist.copy()

    def test_same_seed_same_result_in_any_order(self):
        seventeen, five = self.two_gains(5)
        a17, a5 = self.time_in_order(seventeen, five)
        b5, b17 = self.time_in_order(five, seventeen)
        self.assertTrue(np.array_equal(a17, b17))
        self.assertTrue(np.array_equal(a5, b5))
        self.assertTrue(np.array_equal(a17, self.two_gains(5)[0].pi_dist))
        # Without a seed both draw from the global random state, so the order matters.
        seventeen, five = self.two_gains(None)
        a17, a5 = self.time_in_order(seventeen, five)
        b5, b17 = self.time_in_order(five, seventeen)
        self.assertFalse(np.array_equal(a17, b17))

    def test_without_a_seed_the_global_state_is_used(self):
        event = self.gained(None)
        np.random.seed(3)
        random.seed(3)
        first = event._simulate_p2()
        np.random.seed(3)
        random.seed(3)
        self.assertTrue(np.array_equal(first, event._simulate_p2()))
        self.assertFalse(np.array_equal(first, event._simulate_p2()))

    def test_cli_runs_with_one_seed_are_identical(self):
        cli = CommandLine('run_cli')
        cli.tmp, cli.genome = self.tmp, self.genome
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, 1, 2), seg_row('17', 25000000, 81000000, 1, 1)]
        muts = [('17', 5000000 + i * 1000, 1, 'A', 60, 40) for i in range(12)]
        tables = []
        for _ in range(2):
            self.assertEqual(cli.run_cli(rows, 'XX', '--seed', '3', muts=muts).returncode, 0)
            with open(os.path.join(self.tmp, 'P1.timing.tsv')) as handle:
                tables.append(handle.read())
        self.assertEqual(tables[0], tables[1])


class TerminalArmEvents(TimingCase):
    """Fix 8: terminal partial-arm events such as 1p36 loss, called where the whole-arm call misses them."""
    def chrom(self, chrom, *segments, fraction=.2, male=False, baseline=(1, 1), muts=()):
        rows = [r for r in whole_genome(self.genome, male=male, baseline=baseline) if r[0] != chrom] + list(segments)
        options = TimingOptions(sex='XY' if male else 'XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True,
                                terminal_arm_min_fraction=fraction)
        return self.run_engine(rows, muts=muts, options=options)

    def one_p(self, deleted_to, fraction=.2):
        size = self.genome.CHROM_DICT['1']
        return self.chrom('1', seg_row('1', 1000, deleted_to, 0, 1), seg_row('1', deleted_to + 1, 124000000, 1, 1),
                          seg_row('1', 126000000, size - 1000, 1, 1), fraction=fraction)

    def test_terminal_loss_is_called_under_the_arm_name(self):
        engine = self.one_p(30000000)
        event = self.events(engine)['loss_1p']
        self.assertEqual(event.Type, 'Arm_loss')
        self.assertTrue(event.is_clonal)
        self.assertIsNotNone(event.state.terminal)
        self.assertNotIn('loss_1p', self.events(self.one_p(30000000, fraction=None)))

    def test_short_terminal_run_is_not_called(self):
        self.assertNotIn('loss_1p', self.events(self.one_p(20000000)))
        self.assertIn('loss_1p', self.events(self.one_p(20000000, fraction=.15)))

    def test_interstitial_loss_is_not_terminal(self):
        size = self.genome.CHROM_DICT['1']
        engine = self.chrom('1', seg_row('1', 1000, 40000000, 1, 1), seg_row('1', 40000001, 70000000, 0, 1),
                            seg_row('1', 70000001, 124000000, 1, 1), seg_row('1', 126000000, size - 1000, 1, 1))
        self.assertNotIn('loss_1p', self.events(engine))

    def test_short_noisy_segment_does_not_break_the_run(self):
        size = self.genome.CHROM_DICT['1']
        engine = self.chrom('1', seg_row('1', 1000, 15000000, 0, 1), seg_row('1', 15000001, 16000000, 1, 1),
                            seg_row('1', 16000001, 32000000, 0, 1), seg_row('1', 32000001, 124000000, 1, 1),
                            seg_row('1', 126000000, size - 1000, 1, 1))
        self.assertAlmostEqual(self.events(engine)['loss_1p'].state.terminal[1], 32000001, delta=1)

    def test_half_arm_loss_that_splits_the_whole_arm_call(self):
        size = self.genome.CHROM_DICT['11']
        segments = [seg_row('11', 1000, 51000000, 1, 1), seg_row('11', 54000000, 94000000, 1, 1), seg_row('11', 94000001, size - 1000, 0, 1)]
        self.assertIn('loss_11q', self.events(self.chrom('11', *segments)))
        self.assertNotIn('loss_11q', self.events(self.chrom('11', *segments, fraction=None)))

    def test_no_second_event_when_the_whole_arm_call_has_it(self):
        size = self.genome.CHROM_DICT['11']
        engine = self.chrom('11', seg_row('11', 1000, 51000000, 1, 1), seg_row('11', 54000000, 70000000, 1, 1),
                            seg_row('11', 70000001, size - 1000, 0, 1))
        self.assertEqual(len(engine.all_cn_events['loss_11q']), 1)
        self.assertIsNone(engine.all_cn_events['loss_11q'][0].state.terminal)

    def test_terminal_gain(self):
        size = self.genome.CHROM_DICT['17']
        engine = self.chrom('17', seg_row('17', 1000, 22000000, 1, 1), seg_row('17', 25000000, 55000000, 1, 1),
                            seg_row('17', 55000001, size - 1000, 1, 2))
        self.assertEqual(self.events(engine)['gain_17q'].Type, 'Arm_gain')

    def test_terminal_states_do_not_vote_for_wgd(self):
        size = self.genome.CHROM_DICT['1']
        segments = [seg_row('1', 1000, 30000000, 1, 2), seg_row('1', 30000001, 124000000, 2, 2), seg_row('1', 126000000, size - 1000, 2, 2)]
        engine = self.chrom('1', *segments, baseline=(2, 2), muts=doubled_arm_muts())
        self.assertIsNotNone(engine.WGD)
        self.assertIn('1p' + TERMINAL, engine.concordant_cn_states)
        for wgd in [engine.WGD] + [sample.WGD for sample in engine.sample_list]:
            arms = [(state.chrN, state.arm) for state in wgd.supporting_arm_states]
            self.assertEqual(len(arms), len(set(arms)))
            self.assertEqual([state for state in wgd.supporting_arm_states if (state.chrN, state.arm) == ('1', 'p')][0].cn_a1, 2.)
        event = self.events(engine)['loss_1p']
        self.assertEqual(event.allelic_cn, 1.)
        self.assertGreater(self.pi_mean(event.pi_dist), self.pi_mean(engine.WGD.pi_dist))

    def test_needs_continuous_cn(self):
        with self.assertRaises(ValueError):
            TimingOptions(terminal_arm_min_fraction=.2)

    def test_short_arm_needs_the_minimum_length(self):
        # 18p is 15 Mb: a quarter of it is 3.75 Mb, shorter than the 5 Mb minimum.
        size = self.genome.CHROM_DICT['18']
        segments = [seg_row('18', 1000, 3750000, 0, 1), seg_row('18', 3750001, 15000000, 1, 1), seg_row('18', 17500000, size - 1000, 1, 1)]
        self.assertNotIn('loss_18p', self.events(self.chrom('18', *segments)))
        segments[0], segments[1] = seg_row('18', 1000, 6000000, 0, 1), seg_row('18', 6000001, 15000000, 1, 1)
        self.assertIn('loss_18p', self.events(self.chrom('18', *segments)))

    def test_cli(self):
        cli = CommandLine('run_cli')
        cli.tmp, cli.genome = self.tmp, self.genome
        size = self.genome.CHROM_DICT['1']
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '1']
        rows += [seg_row('1', 1000, 30000000, 0, 1), seg_row('1', 30000001, 124000000, 1, 1), seg_row('1', 126000000, size - 1000, 1, 1)]
        flags = ['--sex_aware', '--continuous_cn', '--clonal_unresolved']
        self.assertEqual(cli.run_cli(rows, 'XX', *flags, '--terminal_arm_events', '0.1').returncode, 0)
        self.assertIn('loss_1p', set(cli.table()['Event Name']))
        self.assertEqual(cli.run_cli(rows, 'XX', *flags, '--terminal_arm_events', '0.1', '--terminal_arm_min_length', '40000000').returncode, 0)
        self.assertNotIn('loss_1p', set(cli.table()['Event Name']))
        self.assertNotEqual(cli.run_cli(rows, 'XX', '--sex_aware', '--terminal_arm_events', '0.1').returncode, 0)


class WholeChromosomeEvents(TimingCase):
    """Fix 8: a change of both arms in one state is one whole-chromosome event."""
    def engine(self, rows, whole=True, male=False, muts=(), focal=None, segment_homdels=False):
        options = TimingOptions(sex='XY' if male else 'XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True,
                                segment_homdels=segment_homdels, whole_chromosome_events=whole)
        return self.run_engine(rows, muts=muts, options=options, focal=focal)

    def seventeen(self, p, q):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        return rows + [seg_row('17', 1000, 22000000, *p), seg_row('17', 25000000, 81000000, *q)]

    def test_whole_chromosome_gain(self):
        events = self.events(self.engine(self.seventeen((1, 2), (1, 2))))
        self.assertEqual(events['gain_17'].Type, 'Chrom_gain')
        self.assertTrue(events['gain_17'].is_clonal)
        self.assertNotIn('gain_17q', events)
        self.assertNotIn('gain_17p', events)
        events = self.events(self.engine(self.seventeen((1, 2), (1, 2)), whole=False))
        self.assertIn('gain_17q', events)
        self.assertNotIn('gain_17', events)

    def test_segmental_gain_keeps_its_arm_name(self):
        events = self.events(self.engine(self.seventeen((1, 1), (1, 2))))
        self.assertIn('gain_17q', events)
        self.assertNotIn('gain_17', events)

    def test_different_states_keep_arm_names(self):
        events = self.events(self.engine(self.seventeen((1, 2), (1, 3))))
        self.assertIn('gain_17p', events)
        self.assertIn('gain_17q', events)
        self.assertNotIn('gain_17', events)

    def test_clonality_is_the_lower_of_the_arms(self):
        events = self.events(self.engine(self.seventeen((1, 1), (1, 1))))
        self.assertNotIn('gain_17', events)
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, 1, 1, q2=1.45), seg_row('17', 25000000, 81000000, 1, 1, q2=1.55)]
        event = self.events(self.engine(rows))['gain_17']
        self.assertAlmostEqual(event.ccf_hat, .45)
        self.assertEqual(event.timing_status, 'subclonal')

    def test_mutations_of_both_arms_time_it(self):
        depth = 100
        def muts(m):
            alt = int(depth * m * PURITY / (3 * PURITY + 2 * (1 - PURITY)))
            return [(c, start + i * 1000, 1, 'A', depth - alt, alt) for c, start in (('17', 5000000), ('17', 40000000)) for i in range(6)]
        early = self.events(self.engine(self.seventeen((1, 2), (1, 2)), muts=muts(2)))['gain_17']
        late = self.events(self.engine(self.seventeen((1, 2), (1, 2)), muts=muts(1)))['gain_17']
        self.assertEqual(len(early.supporting_muts), 12)
        self.assertEqual(early.timing_status, 'estimated')
        self.assertGreater(self.pi_mean(early.pi_dist), self.pi_mean(late.pi_dist) + .3)

    def test_every_sample_must_hold_the_state(self):
        # S2's 17q is subclonal. The patient-level states take each value from the sample closest to the integer state.
        def rows(q2):
            base = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
            return base + [seg_row('17', 1000, 22000000, 1, 2, q2=2.), seg_row('17', 25000000, 81000000, 1, 2, q2=q2)]
        options = TimingOptions(sex='XX', sex_aware=True, continuous_cn=True, clonal_unresolved=True, whole_chromosome_events=True)
        engine = self.run_samples([('S1', rows(2.), (), PURITY), ('S2', rows(1.5), (), PURITY)], options=options)
        events = self.events(engine)
        self.assertNotIn('gain_17', events)
        self.assertTrue(events['gain_17p'].is_clonal)
        self.assertEqual(events['gain_17q'].timing_status, 'subclonal')
        engine = self.run_samples([('S1', rows(2.), (), PURITY), ('S2', rows(2.), (), PURITY)], options=options)
        self.assertTrue(self.events(engine)['gain_17'].is_clonal)

    def test_a_subclonal_arm_makes_it_subclonal(self):
        def merged(p_q2, q_q2):
            rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
            rows += [seg_row('17', 1000, 22000000, 1, 2, q2=p_q2), seg_row('17', 25000000, 81000000, 1, 2, q2=q_q2)]
            return self.events(self.engine(rows))['gain_17']
        # 17p subclonal at 0.55, 17q too close to the cutoff to tell at 0.72.
        event = merged(1.55, 1.72)
        self.assertEqual(event.timing_status, 'subclonal')
        self.assertFalse(event.clonality_ambiguous)
        self.assertAlmostEqual(event.ccf_hat, .55)
        # 17p clonal, 17q unknown at 0.85: unknown.
        event = merged(2., 1.85)
        self.assertEqual(event.timing_status, 'ambiguous_clonality')
        self.assertTrue(event.clonality_ambiguous)
        self.assertFalse(event.is_clonal)

    def test_values_must_be_finite(self):
        engine = self.engine(self.seventeen((1, 2), (1, 2)), whole=False)
        p = engine.concordant_cn_states['17p']
        p.continuous = (float('nan'), p.continuous[1])
        engine.merge_whole_chromosome_events()
        self.assertNotIn('gain_17', self.events(engine))
        p.continuous = (1., p.continuous[1])
        engine.merge_whole_chromosome_events()
        self.assertIn('gain_17', self.events(engine))

    def test_focal_event_finds_the_prior_of_its_merged_arm(self):
        # A focal gain reaching the state of its arm borrows the arm's timing. With both arms merged, that is the
        # timing of the whole-chromosome event.
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '17']
        rows += [seg_row('17', 1000, 22000000, 1, 2), seg_row('17', 25000000, 55000000, 1, 2), seg_row('17', 55000001, 57000000, 1, 1),
                 seg_row('17', 57000001, 57500000, 1, 2), seg_row('17', 57500001, 59500000, 1, 1), seg_row('17', 59500001, 81000000, 1, 2)]
        depth = 100
        alt = int(depth * 2 * PURITY / (3 * PURITY + 2 * (1 - PURITY)))
        muts = [('17', start + i * 1000, 1, 'A', depth - alt, alt) for start in (5000000, 40000000) for i in range(6)]
        events = self.events(self.engine(rows, muts=muts, focal=[('17', 57100000, 57200000, 'AMP')]))
        focal = events['Focal_gain_17:57100000-57200000']
        self.assertEqual((focal.timing_status, focal.timing_source), ('prior_only', 'arm_prior'))
        self.assertTrue(np.allclose(focal.pi_dist, events['gain_17'].pi_dist))

    def test_integer_states_merge(self):
        # Without continuous_cn the events have no allele, and a 2/2 arm holds two gains of one type.
        options = TimingOptions(sex='XX', sex_aware=True, clonal_unresolved=True, whole_chromosome_events=True)
        engine = self.run_engine(self.seventeen((2, 2), (2, 2)), options=options)
        self.assertEqual([e.Type for e in engine.all_cn_events['gain_17']], ['Chrom_gain', 'Chrom_gain'])
        self.assertEqual(engine.all_cn_events['gain_17p'] + engine.all_cn_events['gain_17q'], [])
        engine = self.run_engine(self.seventeen((1, 2), (1, 2)), options=options)
        self.assertEqual(len(engine.all_cn_events['gain_17']), 1)

    def test_male_x_and_y(self):
        rows = [r for r in whole_genome(self.genome, male=True) if r[0] != 'X'] + y_rows(0.)
        size = self.genome.CHROM_DICT['X']
        rows += [hemizygous_row('X', 1000, 58000000, 2.), hemizygous_row('X', 62000000, size - 1000, 2.)]
        events = self.events(self.engine(rows, male=True))
        self.assertEqual(events['gain_X'].Type, 'Chrom_gain')
        self.assertEqual(events['gain_X'].allele, 'a2')
        self.assertNotIn('gain_Xp', events)
        self.assertEqual(events['loss_Y'].Type, 'Chrom_loss')
        self.assertTrue(events['loss_Y'].is_clonal)

    def test_homdel_is_unaffected(self):
        rows = [r for r in whole_genome(self.genome, male=False) if r[0] != '9']
        size = self.genome.CHROM_DICT['9']
        rows += [seg_row('9', 1000, 21899999, 0, 1), seg_row('9', 21900000, 22100000, 0, 0, q1=.01, q2=.01),
                 seg_row('9', 22100001, 38000000, 0, 1), seg_row('9', 71000000, size - 1000, 0, 1)]
        events = self.events(self.engine(rows, focal=[('9', 21000000, 23000000, 'DEL')], segment_homdels=True))
        self.assertIn('Focal_homdel_9:21000000-23000000', events)
        self.assertIn('loss_9', events)
        self.assertNotIn('loss_9p', events)


class GainOfY(TimingCase):
    """Fix 8: chrY gains are emitted as gain_Y, relative to one copy, or two after a WGD."""
    def test_two_copies_without_wgd(self):
        rows = whole_genome(self.genome, male=True) + y_rows(1.95)
        event = self.events(self.run_engine(rows, options=TimingOptions(sex='XY', **ALL_ON)))['gain_Y']
        self.assertEqual(event.Type, 'Chrom_gain')
        self.assertEqual(event.allelic_cn, 2.)
        self.assertTrue(event.is_clonal)
        self.assertEqual(event.timing_status, 'clonal_unresolved')

    def test_after_wgd(self):
        rows = whole_genome(self.genome, male=True, baseline=(2, 2)) + y_rows(3.)
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=TimingOptions(sex='XY', **ALL_ON))
        event = self.events(engine)['gain_Y']
        self.assertEqual(event.allelic_cn, 3.)
        self.assertEqual(event.wgd_side, 'after')
        rows = whole_genome(self.genome, male=True, baseline=(2, 2)) + y_rows(2.)
        engine = self.run_engine(rows, muts=doubled_arm_muts(), options=TimingOptions(sex='XY', **ALL_ON))
        self.assertEqual([name for name in self.events(engine) if 'Y' in name], [])


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
        self.assertAlmostEqual(rows['c']['draw_fraction'], 1.)
        self.assertAlmostEqual(rows['b']['n_unknown_games'], .1)

    def test_scores_keep_the_rescue(self):
        pair = self.league().event_pairs[('b', 'c')]
        self.assertEqual(pair.observed_rates['b'] + pair.observed_rates['c'], 0)
        self.assertEqual(pair.win_rates['b'] + pair.win_rates['c'], 2)


if __name__ == '__main__':
    unittest.main()
