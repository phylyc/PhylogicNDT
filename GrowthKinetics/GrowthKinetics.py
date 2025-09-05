import logging
import gzip
from collections import defaultdict


def run_tool(args):
    logging.debug('Arguments {}'.format(args))
    import data.Patient as Patient
    from .GrowthKineticsEngine import GrowthKineticsEngine

    patient_data = Patient.Patient(indiv_name=args.indiv_id)

    if args.sif:  # if sif file is specified
        with open(args.sif, 'r') as sif_file:
            header = sif_file.readline().strip('\n').split('\t')
            for line in sif_file:
                row = dict(zip(header, line.strip('\n').split('\t')))
                sample_id = row['sample_id']
                maf_fn = row['maf_fn']
                seg_fn = row['seg_fn']
                purity = float(row['purity'])
                timepoint = float(row['timepoint'])
                tmb = float(row.get("tmb", 1))
                patient_data.addSample(maf_fn, sample_id, timepoint_value=timepoint,
                                       _additional_muts=None, seg_file=seg_fn,
                                       purity=purity, tmb=tmb)

    mcmc_trace_cell_abundance, num_itertaions = load_mcmc_trace_abundances(args.abundance_mcmc_trace)
    gk_engine = GrowthKineticsEngine(patient_data, args.wbc)
    gk_engine.estimate_growth_rate(mcmc_trace_cell_abundance, times=args.time, n_iter=min(num_itertaions, args.n_iter))

    # Output and visualization
    import output.PhylogicOutput
    phylogicoutput = output.PhylogicOutput.PhylogicOutput()
    phylogicoutput.write_growth_rate_tsv(gk_engine.growth_rates, args.indiv_id)
    phylogicoutput.plot_growth_rates(gk_engine.growth_rates, args.indiv_id)


def load_mcmc_trace_abundances(in_file):
    iterations = set()
    cell_abundance_mcmc_trace = defaultdict(lambda: defaultdict(list))
    open_func = gzip.open if in_file.endswith(".gz") else open
    with open_func(in_file, 'rt') as reader:
        for line in reader:
            values = line.strip('\n').split('\t')
            if line.startswith('Patient_ID'):
                header = {k: v for v, k in enumerate(values)}
            else:
                sample_id = values[header['Sample_ID']]
                cluster_id = int(values[header['Cluster_ID']])
                abundance = int(values[header['Abundance']])
                iterations.add(values[header['Iteration']])
                cell_abundance_mcmc_trace[sample_id][cluster_id].append(abundance)
    return cell_abundance_mcmc_trace, len(iterations)