"""Check that Timing and TimingComparison give identical tables with all new flags off, on master and on this branch.

Usage: python tests/regress_timing_flags_off.py <checkout> <output_dir> <patient_dir> [<patient_dir> ...]
Run it once per checkout with PYTHONHASHSEED=0, then compare the two output folders. Patients come from an existing
somatic_workflow output folder with a REP0 timing SIF and mut_ccfs. FOCAL_REGIONS names the --cn_peaks file.
"""
import os
import random
import runpy
import sys

import numpy as np
import pandas as pd


def run_cli(checkout, arguments):
    random.seed(7)
    np.random.seed(7)
    sys.argv = [os.path.join(checkout, 'PhylogicNDT.py')] + arguments
    try:
        runpy.run_path(sys.argv[0], run_name='__main__')
    except SystemExit as stop:
        if stop.code not in (None, 0):
            raise


def main():
    checkout, output_dir = map(os.path.abspath, sys.argv[1:3])
    sys.path.insert(0, checkout)
    os.chdir(checkout)
    for patient_dir in map(os.path.abspath, sys.argv[3:]):
        patient = os.path.basename(patient_dir)
        out = os.path.join(output_dir, patient)
        os.makedirs(out, exist_ok=True)
        sif = pd.read_csv(os.path.join(patient_dir, f'{patient}.REP0.sif.timing.txt.gz'), sep='\t')
        sif['maf_fn'] = [os.path.join(patient_dir, os.path.basename(p)) for p in sif['maf_fn']]
        sif['seg_fn'] = [os.path.join(patient_dir, os.path.basename(p)) for p in sif['seg_fn']]
        sif_path = os.path.join(out, 'input.sif')
        sif.to_csv(sif_path, sep='\t', index=False)
        os.chdir(out)
        run_cli(checkout, ['Timing', '-i', patient, '-sif', sif_path, '--cn_peaks', os.environ['FOCAL_REGIONS']])
        names = pd.read_csv(f'{patient}.timing.tsv', sep='\t', usecols=['Event Name'], low_memory=False)['Event Name']
        names.drop_duplicates().head(300).to_csv('events.txt', index=False, header=False)
        run_cli(checkout, ['TimingComparison', '-i', patient, '--timing_tsv', f'{patient}.timing.tsv', '--event_file', 'events.txt'])
        os.chdir(checkout)
        print(patient, 'done', flush=True)


if __name__ == '__main__':
    main()
