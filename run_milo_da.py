#!/usr/bin/env python
"""
run_milo_da.py -- Step 3 of the generalized NICHES+Milo pipeline.

Reads the NICHES cell-to-cell network h5ad from run_niches_generic.R, runs
Milo differential-abundance testing for the strict two-group comparison set
up in prep_niches_input.py, Louvain-clusters the significant neighborhoods,
and writes the final annotated h5ad (with `sc_louvain` cluster assignments)
that downstream enrichment/figure scripts read.

Run with `uv run` from the repository root so that `milo_helpers` (the
shared library module at the repo root) is importable:

    uv run python run_milo_da.py \\
        --input data/my_dataset_NICHES.h5ad \\
        --output data/my_dataset_MILO.h5ad \\
        --wilcox-out data/my_dataset_wilcoxon.xlsx

See PIPELINE_STEPS.md for the full step-by-step reference (inputs/outputs
for every stage of the pipeline) and README.md for setup instructions.
"""
import argparse
import os
import pickle
import sys

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import milopy
import milopy.core as milo

# Newer anndata added a strict guard against writing pandas' newer nullable
# string arrays (pd.StringDtype) to h5ad, off by default, because older
# anndata (<0.11) can't read them back -- see the write_nullable error this
# opts out of. Somewhere in this pipeline's R/Python round-trips an obs
# index ends up as that newer dtype; this is anndata's own documented
# opt-in for writing it anyway, not a workaround for a bug.
anndata.settings.allow_write_nullable_strings = True

from milo_helpers import plot_nhood_clusters
from _lib import (
    write_wilcoxon_results, make_nhoods_fixed, DA_nhoods_fixed,
    count_nhoods_fixed, cluster_nhoods_fixed, get_sc_louvain_fixed,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True,
                    help="NICHES h5ad from run_niches_generic.R.")
    p.add_argument("--output", required=True,
                    help="Output path for the final Milo-annotated h5ad.")
    p.add_argument("--split-col", default="Condition",
                    help="obs column NICHES tagged with the split-unit "
                         "string (default: 'Condition').")
    p.add_argument("--condition-delim", default="__",
                    help="Delimiter used in --split-col to separate the 0/1 "
                         "condition prefix from the sample id -- MUST match "
                         "what prep_niches_input.py used (default: '__').")
    p.add_argument("--annotation-col", default="celltype.Joint",
                    help="obs column used to label each neighborhood by its "
                         "dominant sender/receiver cell types (default: "
                         "'celltype.Joint', created automatically by NICHES "
                         "when a matching celltype field is in "
                         "meta.data.to.map). Use 'VectorType' if you mapped "
                         "your celltype column under a different name.")
    p.add_argument("--min-connect", type=int, default=2,
                    help="Minimum shared-cell connections to merge two DA "
                         "neighborhoods into the same cluster (default: 2, "
                         "matches the published analyses).")
    p.add_argument("--max-difflfc", type=float, default=2.2,
                    help="Maximum logFC difference to still merge two DA "
                         "neighborhoods (default: 2.2, matches the "
                         "published analyses).")
    p.add_argument("--alpha", type=float, default=0.1,
                    help="SpatialFDR threshold for the neighborhood plot "
                         "(default: 0.1).")
    p.add_argument("--min-size", type=int, default=5,
                    help="Point-size scaling factor for the neighborhood "
                         "plot (default: 5).")
    p.add_argument("--plot-out", default=None,
                    help="If given, save the neighborhood cluster plot to "
                         "this path (e.g. clusters.pdf). Skipped if omitted.")
    p.add_argument("--wilcox-out", default=None,
                    help="If given, run Wilcoxon rank-sum DE across "
                         "sc_louvain clusters and write results to this "
                         ".xlsx path (one sheet per cluster).")
    p.add_argument("--adata-checkpoint", default=None,
                    help="Path to cache the AnnData right after Milo's "
                         "differential-abundance test completes (the step "
                         "that depends on the rpy2/R/edgeR bridge) -- saved "
                         "as a pickle, not h5ad, since adata.uns['nhood_adata'] "
                         "is itself a nested AnnData that h5ad can't round-trip "
                         "natively. If this file already exists, it's loaded "
                         "instead of recomputing neighbors/nhoods/DA, so a "
                         "rerun after a later failure (clustering, Wilcoxon, "
                         "the final write) doesn't need the R bridge again. "
                         "Default: '<output>_da.pkl' next to --output. Pass "
                         "an empty string ('') to disable.")
    return p.parse_args(argv)


def run_milo_generic(adata, split_col, condition_delim, annotation_col):
    """Two-group Milo DA, generalized from milo_helpers.run_milo().

    Recovers the 0/1 condition encoding from `split_col` (values look like
    "<0-or-1><delim><sample id>", written by prep_niches_input.py and
    propagated onto the NICHES output's 'Condition' field by
    run_niches_generic.R), counts cells per split-unit ("pseudo-replicate")
    in each neighborhood, and tests for differential abundance with design
    `~ cond_binary`.
    """
    if split_col not in adata.obs.columns:
        sys.exit(f"--split-col '{split_col}' not found in adata.obs. "
                  f"Available columns: {list(adata.obs.columns)}")

    cond_prefix = adata.obs[split_col].astype(str).str.split(condition_delim).str[0]
    try:
        adata.obs["cond_binary"] = cond_prefix.astype(int)
    except ValueError:
        sys.exit(
            f"Could not parse a 0/1 condition prefix out of '{split_col}' "
            f"using delimiter {condition_delim!r}. Example values seen: "
            f"{adata.obs[split_col].astype(str).unique()[:5].tolist()}. "
            "Does --condition-delim match what prep_niches_input.py used?"
        )

    uniq = sorted(adata.obs["cond_binary"].unique())
    if uniq != [0, 1]:
        sys.exit(
            f"Expected exactly two condition groups (0 and 1) after parsing "
            f"'{split_col}', found: {uniq}. This pipeline is strict "
            "two-group only -- see README.md if you need >2 groups."
        )

    adata.obs[split_col] = adata.obs[split_col].astype("category")
    adata.obs["rep_code"] = adata.obs[split_col].cat.codes

    sc.pp.neighbors(adata)
    # milo.make_nhoods() itself is broken against current pandas (see
    # _lib.make_nhoods_fixed's docstring) -- using the patched copy instead.
    make_nhoods_fixed(adata)
    # milo.count_nhoods() itself is also broken: pd.get_dummies() now
    # defaults to bool dtype, which doesn't survive as numeric once it
    # reaches edgeR via rpy2 (see _lib.count_nhoods_fixed's docstring) --
    # using the patched copy instead.
    count_nhoods_fixed(adata, sample_col="rep_code")
    # milo.DA_nhoods() itself is also broken: it opens with rpy2's
    # numpy2ri.activate()/pandas2ri.activate(), which newer rpy2 turns into
    # a hard `raise DeprecationWarning(...)` instead of activating anything
    # (see _lib.DA_nhoods_fixed's docstring) -- using the patched copy.
    DA_nhoods_fixed(adata, design="~ cond_binary")
    milopy.utils.build_nhood_graph(adata)

    if annotation_col not in adata.obs.columns:
        sys.exit(
            f"--annotation-col '{annotation_col}' not found in adata.obs. "
            "NICHES creates this automatically from a mapped celltype field "
            "when its Sending and Receiving values agree -- check that "
            "meta.data.to.map in run_niches_generic.R included your "
            "celltype column, or pass --annotation-col to point at the "
            "right field (e.g. 'VectorType')."
        )
    milopy.utils.annotate_nhoods(adata, anno_col=annotation_col)
    nhood_obs = adata.uns["nhood_adata"].obs
    # milopy.utils.annotate_nhoods() leaves "nhood_annotation" as a pandas
    # Categorical whose categories are exactly the values seen in
    # --annotation-col. Modern pandas refuses to assign a value that isn't
    # already a declared category (older pandas added it silently) --
    # "Mixed" needs to be added as a category first.
    if isinstance(nhood_obs["nhood_annotation"].dtype, pd.CategoricalDtype) and \
            "Mixed" not in nhood_obs["nhood_annotation"].cat.categories:
        nhood_obs["nhood_annotation"] = nhood_obs["nhood_annotation"].cat.add_categories(["Mixed"])
    nhood_obs.loc[nhood_obs["nhood_annotation_frac"] < 0.5, "nhood_annotation"] = "Mixed"

    return adata


def main(argv=None):
    args = parse_args(argv)

    adata_checkpoint = args.adata_checkpoint
    if adata_checkpoint is None:
        adata_checkpoint = os.path.splitext(args.output)[0] + "_da.pkl"
    if adata_checkpoint == "":
        adata_checkpoint = None

    if adata_checkpoint and os.path.exists(adata_checkpoint):
        print(f"Found post-DA checkpoint, loading instead of recomputing: {adata_checkpoint}")
        with open(adata_checkpoint, "rb") as fh:
            adata = pickle.load(fh)
    else:
        adata = sc.read_h5ad(args.input)
        adata = run_milo_generic(adata, args.split_col, args.condition_delim, args.annotation_col)
        if adata_checkpoint:
            with open(adata_checkpoint, "wb") as fh:
                pickle.dump(adata, fh)
            print(f"Saved post-DA checkpoint: {adata_checkpoint}")

    print("Clustering differentially abundant neighborhoods...")
    # milo_helpers.cluster_nhoods() itself is broken against current
    # networkx (nx.from_numpy_matrix was removed -- see
    # _lib.cluster_nhoods_fixed's docstring) -- using the patched copy.
    partition = cluster_nhoods_fixed(adata, args.min_connect, args.max_difflfc)

    # milo_helpers.get_sc_louvain() (used below) expects a 'louvain' column
    # on adata.uns['nhood_adata'].obs: the per-nhood cluster label from
    # `partition`, with non-significant neighborhoods (SpatialFDR > alpha)
    # marked as float('inf') so they're excluded from the cluster mapping.
    # The published analysis scripts call get_sc_louvain() without ever
    # setting this column explicitly -- plot_nhood_clusters() only sets an
    # equivalent column on a throwaway *copy* of nhood_adata, so it never
    # reaches the object get_sc_louvain() reads from. This looks like a step
    # that dropped out of the reference scripts; we set it explicitly here
    # so get_sc_louvain() has something real to read.
    nhood_adata = adata.uns["nhood_adata"]
    nhood_adata.obs["louvain"] = np.array(list(partition.values()), dtype=float)
    nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"] > args.alpha, "louvain"] = float("inf")

    plot_nhood_clusters(
        adata, list(partition.values()), "Louvain cluster",
        alpha=args.alpha, min_size=args.min_size,
    )
    if args.plot_out:
        import matplotlib.pyplot as plt
        plt.savefig(args.plot_out, bbox_inches="tight")
        print(f"Saved cluster plot to {args.plot_out}")

    # milo_helpers.get_sc_louvain() itself is broken for a real edge case
    # (zero significant neighborhoods -- see _lib.get_sc_louvain_fixed's
    # docstring) -- using the patched copy.
    adata.obs["sc_louvain"] = get_sc_louvain_fixed(adata)
    adata.obs["sc_louvain"] = adata.obs["sc_louvain"].astype("category")
    adata.obs["louvain_str"] = list(map(str, adata.obs["sc_louvain"].values))

    n_clusters = len(adata.obs["sc_louvain"].unique())
    print(f"Found {n_clusters} sc_louvain cluster label(s) (including -1 = unassigned).")

    if args.wilcox_out:
        print("Running Wilcoxon DE across clusters...")
        sc.tl.rank_genes_groups(adata, "louvain_str", method="wilcoxon", key_added="louvain-wilc")
        write_wilcoxon_results(args.wilcox_out, adata, np.unique(adata.obs["sc_louvain"]), "louvain-wilc")
        print(f"Wrote {args.wilcox_out}")

    adata.write(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
