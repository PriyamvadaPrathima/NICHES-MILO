#!/usr/bin/env python
"""
prep_niches_input.py -- Step 1 of the generalized NICHES+Milo pipeline.

Subsets an already QC'd / normalized / cell-type-annotated h5ad down to the
cell types you want in the communication network, restricts to a strict
two-group comparison, and adds three bookkeeping columns the rest of the
pipeline depends on:

    _condition_binary   0/1 encoding of the two-group comparison
    _niches_split       the unit NICHES imputes+networks per sample; ENCODES
                         the condition (see PIPELINE_STEPS.md, "Why the split
                         unit encodes the condition") so step 3 can recover
                         cond_binary without a separate lookup table
    _grouping           coarse cell-type bucket passed to NICHES as its
                         `cell_types` argument (sender/receiver identity)

Run with `uv run` from the repository root, e.g.:

    uv run python prep_niches_input.py \\
        --input data/my_dataset.h5ad \\
        --output data/my_dataset_LIM.h5ad \\
        --sample-col sample_id \\
        --condition-col treatment \\
        --reference-value control \\
        --treatment-value drug \\
        --celltype-col cell_type \\
        --cell-types Macrophage DC "T cell"

See PIPELINE_STEPS.md for the full step-by-step reference (inputs/outputs
for every stage of the pipeline) and README.md for setup instructions.
"""
import argparse
import json
import sys

import scanpy as sc


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True,
                    help="Input h5ad. Expected to already be QC'd, normalized, "
                         "and cell-type annotated -- this script does not do QC.")
    p.add_argument("--output", required=True,
                    help="Output path for the subsetted h5ad that "
                         "run_niches_generic.R will read.")
    p.add_argument("--sample-col", required=True,
                    help="obs column with the sample/batch ID -- the biological "
                         "replicate unit NICHES imputes separately and Milo "
                         "later counts nhood cells per.")
    p.add_argument("--condition-col", required=True,
                    help="obs column with the treatment/condition label.")
    p.add_argument("--reference-value", required=True,
                    help="Value in --condition-col treated as the reference/"
                         "control group (encoded 0).")
    p.add_argument("--treatment-value", required=True,
                    help="Value in --condition-col treated as the treatment "
                         "group (encoded 1).")
    p.add_argument("--celltype-col", default="celltype",
                    help="obs column with fine-grained per-cell type labels "
                         "(default: 'celltype').")
    p.add_argument("--cell-types", nargs="+", required=True,
                    help="Which values of --celltype-col to keep for the "
                         "communication analysis (the sender/receiver "
                         "populations -- e.g. Macrophage DC 'T cell').")
    p.add_argument("--grouping-col", default=None,
                    help="obs column with a coarse cell-type grouping you've "
                         "already computed. If omitted, --grouping-map is used "
                         "if given, else --celltype-col values are used "
                         "directly (no collapsing).")
    p.add_argument("--grouping-map", default=None,
                    help="Path to a JSON file mapping {celltype_value: "
                         "group_name}, used to build the coarse grouping when "
                         "--grouping-col is not supplied. See "
                         "example_grouping_map.json for the format.")
    p.add_argument("--replicate-col", default=None,
                    help="Optional obs column with a finer replicate/"
                         "multiplexing ID (e.g. a hashtag oligo) to combine "
                         "with --sample-col before NICHES imputation. Omit if "
                         "--sample-col alone already identifies each "
                         "imputation/replicate unit.")
    p.add_argument("--min-cells-per-split", type=int, default=20,
                    help="Warn if any (sample[, replicate]) split unit has "
                         "fewer than this many cells after subsetting "
                         "(default: 20). ALRA imputation is unreliable on "
                         "very small groups.")
    return p.parse_args(argv)


def build_grouping(adata, celltype_col, grouping_col, grouping_map_path):
    if grouping_col is not None:
        if grouping_col not in adata.obs.columns:
            sys.exit(f"--grouping-col '{grouping_col}' not found in adata.obs. "
                      f"Available columns: {list(adata.obs.columns)}")
        return adata.obs[grouping_col].astype(str)

    if grouping_map_path is not None:
        with open(grouping_map_path) as fh:
            mapping = json.load(fh)
        present = set(adata.obs[celltype_col].astype(str).unique())
        missing = present - set(mapping)
        if missing:
            sys.exit(f"--grouping-map ({grouping_map_path}) is missing entries "
                      f"for celltype value(s) present in the data: {sorted(missing)}")
        return adata.obs[celltype_col].astype(str).map(mapping)

    # default: no collapsing, grouping == fine celltype
    return adata.obs[celltype_col].astype(str)


def main(argv=None):
    args = parse_args(argv)
    adata = sc.read_h5ad(args.input)

    for col in (args.sample_col, args.condition_col, args.celltype_col):
        if col not in adata.obs.columns:
            sys.exit(f"Column '{col}' not found in adata.obs. "
                      f"Available columns: {list(adata.obs.columns)}")
    if args.replicate_col is not None and args.replicate_col not in adata.obs.columns:
        sys.exit(f"--replicate-col '{args.replicate_col}' not found in adata.obs.")

    cond_vals = set(adata.obs[args.condition_col].astype(str).unique())
    for v in (args.reference_value, args.treatment_value):
        if v not in cond_vals:
            sys.exit(f"--condition-col '{args.condition_col}' does not contain "
                      f"value '{v}'. Found values: {sorted(cond_vals)}")
    if args.reference_value == args.treatment_value:
        sys.exit("--reference-value and --treatment-value must differ.")

    # restrict to exactly the two groups being compared
    keep = adata.obs[args.condition_col].astype(str).isin(
        [args.reference_value, args.treatment_value]
    )
    n_dropped = int((~keep).sum())
    if n_dropped:
        print(f"Dropping {n_dropped} cells outside the two-group comparison "
              f"({args.reference_value!r} vs {args.treatment_value!r}).")
    adata = adata[keep].copy()

    # restrict to cell types of interest
    present_types = set(adata.obs[args.celltype_col].astype(str).unique())
    ct_missing = set(args.cell_types) - present_types
    if ct_missing:
        sys.exit(f"--cell-types value(s) not found in '{args.celltype_col}' "
                  f"(within the two-group subset): {sorted(ct_missing)}. "
                  f"Available: {sorted(present_types)}")
    adata = adata[adata.obs[args.celltype_col].astype(str).isin(args.cell_types)].copy()
    if adata.n_obs == 0:
        sys.exit("No cells remain after filtering to --cell-types. "
                  "Check --celltype-col / --cell-types.")

    # 0/1 condition encoding
    binary_map = {args.reference_value: 0, args.treatment_value: 1}
    adata.obs["_condition_binary"] = (
        adata.obs[args.condition_col].astype(str).map(binary_map).astype(int)
    )

    # coarse grouping for NICHES cell_types=
    adata.obs["_grouping"] = build_grouping(
        adata, args.celltype_col, args.grouping_col, args.grouping_map
    ).values

    # split unit for per-sample ALRA imputation; embeds the condition so
    # step 3 can recover cond_binary from the 'Condition' tag NICHES writes.
    if args.replicate_col is not None:
        split_id = (
            adata.obs[args.sample_col].astype(str) + "-"
            + adata.obs[args.replicate_col].astype(str)
        )
    else:
        split_id = adata.obs[args.sample_col].astype(str)
    adata.obs["_niches_split"] = (
        adata.obs["_condition_binary"].astype(str) + "__" + split_id
    )

    counts = adata.obs["_niches_split"].value_counts()
    small = counts[counts < args.min_cells_per_split]
    if len(small):
        print(f"WARNING: the following split units have fewer than "
              f"{args.min_cells_per_split} cells; ALRA imputation may be "
              f"unreliable for them:")
        print(small.to_string())

    print("\nCells per split unit ('condition__sample[-replicate]'):")
    print(counts.sort_index().to_string())
    print(f"\nTotal cells retained: {adata.n_obs}")
    print(f"Cell types retained ({args.celltype_col}): {args.cell_types}")
    print(f"Grouping values passed to NICHES (_grouping): "
          f"{sorted(adata.obs['_grouping'].unique())}")

    adata.write(args.output)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
