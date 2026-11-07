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
from SinglePatientTiming.TimingEngine import TimingEngine, TimingOptions, Y_CALLABLE_BLOCKS_HG19

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


class FocalEvents(TimingCase):
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
    def run_cli(self, rows, sex, *flags, focal=None):
        maf_path, seg_path = self.write_inputs('S1', rows, ())
        sif = os.path.join(self.tmp, 'input.sif')
        table = {'sample_id': ['S1'], 'maf_fn': [maf_path], 'seg_fn': [seg_path], 'purity': [PURITY], 'timepoint': [0], 'tmb': [1]}
        if sex is not None:
            table['sex'] = [sex]
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
