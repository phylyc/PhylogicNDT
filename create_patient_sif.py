import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient_id", required=True)
    parser.add_argument("--sample_names", nargs="*", required=False)
    parser.add_argument("--absolute_mafs", nargs="+", required=True)
    parser.add_argument("--absolute_segtabs", nargs="*", required=False)
    parser.add_argument("--absolute_purities", nargs="+", required=True, type=float)
    parser.add_argument("--timepoints", nargs='*', type=int, required=False, help="Optional explicit timepoint ordering (one per sample)")
    parser.add_argument("--tumor_mutation_burdens", nargs='*', type=float, required=False, help="Optional tumor mutation burdens (one per sample)")
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    # If args.sample_names not provided, assume sample_ids is just the MAF filename minus ".ABS_MAF.txt"
    sample_ids = [os.path.basename(x).replace(".ABS_MAF.txt", "") for x in args.absolute_mafs] if not args.sample_names else args.sample_names
    n = len(sample_ids)

    maf_fns = args.absolute_mafs
    seg_fns = args.absolute_segtabs if args.absolute_segtabs else [""] * n
    purities = args.absolute_purities
    timepoints = args.timepoints if args.timepoints else [0] * n # Default to timepoint to all 0 if not provided
    tmb = args.tumor_mutation_burdens if args.tumor_mutation_burdens else [1] * n

    if len(maf_fns) != n:
        raise ValueError(f"Number of maf_fns ({len(maf_fns)}) does not match number of samples ({n}).")
    if len(seg_fns) != n:
        raise ValueError(f"Number of seg_fns ({len(seg_fns)}) does not match number of samples ({n}).")
    if len(purities) != n:
        raise ValueError(f"Number of purities ({len(purities)}) does not match number of samples ({n}).")
    if len(timepoints) != n:
        raise ValueError(f"Number of timepoints ({len(timepoints)}) does not match number of samples ({n}).")
    if len(tmb) != n:
        raise ValueError(f"Number of tumor_mutation_burdens ({len(tmb)}) does not match number of samples ({n}).")

    df = pd.DataFrame({
        "sample_id": sample_ids,
        "maf_fn": maf_fns,
        "seg_fn": seg_fns,
        "purity": purities,
        "timepoint": timepoints,
        "tmb": tmb
    })
    df.to_csv(args.outfile, sep="\t", index=False)

    # SIF for timing of mutations in patient based on clustering results
    df = pd.DataFrame({
        "sample_id": sample_ids,
        "maf_fn": [f"{args.patient_id}.mut_ccfs.txt"] * n,
        "seg_fn": seg_fns,
        "purity": purities,
        "timepoint": timepoints,
        "tmb": tmb
    })
    df.to_csv(args.outfile + ".timing.txt", sep="\t", index=False)


if __name__ == "__main__":
    main()
