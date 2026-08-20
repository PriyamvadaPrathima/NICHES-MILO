#!/usr/bin/env python
"""
plot_lr_heatmap.py -- L-R enrichment heatmap for NICHES+Milo results,
generalized from Bridges et al.'s figures/Fig3_heatmap.py (see
milo_helpers.py's header for attribution).

Uses sc.pl.matrixplot the same way the original figure script does:
displayed L-R pairs (var_names) as rows, sc_louvain clusters as columns,
standard_scale='var' so each row is min-max scaled to [0, 1] across the
displayed clusters (a value of 1 marks that pair's "home" cluster).

Unlike the original -- which hand-picks a curated LR_PAIRS list and
HIGHLIGHT_CLUSTERS after visually inspecting one specific comparison --
this version can select both automatically:

  - Clusters: all significant sc_louvain clusters by default (everything
    except '-1' = unassigned), or pass --clusters to pick specific ones,
    same as the original's HIGHLIGHT_CLUSTERS.

  - L-R pairs: top --top-n-per-cluster pairs by Wilcoxon score per
    displayed cluster, read straight out of adata.uns[--wilcox-key] (the
    rank_genes_groups results run_milo_da.py stores there when called with
    --wilcox-out), or pass --lr-pairs for an explicit curated list, same
    as the original's LR_PAIRS.

Two more things this version adds on top of the original figure script:

  - Columns are labeled '<cluster number>: <dominant VectorType>' by
    default (e.g. '30: Macrophage-Treg') instead of a bare Louvain cluster
    number, since that number is an arbitrary community-detection ID with
    no inherent meaning. Pass --no-celltype-labels to see raw numbers.

  - Before plotting, prints each displayed cluster's own Milo
    differential-abundance stats (mean logFC, mean/max SpatialFDR across
    its neighborhoods) -- this is the actual reference-vs-treatment significance,
    which is NOT the same thing as the heatmap's coloring (a per-L-R-pair
    cluster-vs-rest Wilcoxon test, describing what defines each cluster's
    signaling identity, not whether that cluster differs between
    conditions).

Run with `uv run` from the repository root:

    uv run python plot_lr_heatmap.py \\
        --input data/my_dataset_MILO.h5ad \\
        --output data/my_dataset_lr_heatmap.png

    # Or, to reproduce the original figure's exact curated-list style:
    uv run python plot_lr_heatmap.py \\
        --input data/my_dataset_MILO.h5ad \\
        --output data/my_dataset_lr_heatmap.png \\
        --clusters 25 20 23 28 -1 \\
        --lr-pairs "Ctla4—Cd86" "Cd28—Cd86" "Il10—Il10ra"

See PIPELINE_STEPS.md / README.md for the rest of the pipeline.
"""
import argparse
import sys

import scanpy as sc

from _lib import (
    select_clusters, dominant_celltype_labels, print_cluster_significance,
    auto_select_lr_pairs,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, help="NICHES+Milo h5ad from run_milo_da.py.")
    p.add_argument("--output", required=True, help="Output image path (e.g. heatmap.png).")
    p.add_argument("--groupby", default="louvain_str",
                    help="obs column to group by (default: 'louvain_str', "
                         "written by run_milo_da.py).")
    p.add_argument("--clusters", nargs="+", default=None,
                    help="Which cluster labels to display, in display order "
                         "(default: all clusters present except '-1' = "
                         "unassigned, sorted numerically where possible).")
    p.add_argument("--include-unassigned", action="store_true",
                    help="Include the '-1' (unassigned) group as a reference "
                         "column even when --clusters isn't given.")
    p.add_argument("--top-n-clusters", type=int, default=None,
                    help="If given (and --clusters isn't), narrow the "
                         "default cluster list down to this many, ranked by "
                         "--rank-clusters-by. Ignored if --clusters is set.")
    p.add_argument("--rank-clusters-by", choices=["size", "spatialfdr"], default="size",
                    help="How to rank clusters for --top-n-clusters: 'size' "
                         "= most cell pairs assigned (default, always "
                         "available); 'spatialfdr' = most significant on "
                         "average (needs adata.uns['nhood_adata'] with a "
                         "'louvain' column, i.e. this must be the direct "
                         "output of run_milo_da.py).")
    p.add_argument("--celltype-filter", default=None,
                    help="Restrict to clusters whose DOMINANT --celltype-col "
                         "value contains this text (case-insensitive), e.g. "
                         "--celltype-filter Macrophage for every cluster "
                         "that's mostly macrophage-involving cell pairs. "
                         "Applied before --top-n-clusters, so combining the "
                         "two gives the top N *matching* clusters.")
    p.add_argument("--min-cluster-size", type=int, default=None,
                    help="Drop clusters with fewer than this many cell "
                         "pairs before display -- a quality floor against "
                         "clusters Milo called significant on a handful of "
                         "cell pairs. Applied before --top-n-clusters, same "
                         "as --celltype-filter, and even against clusters "
                         "named explicitly via --clusters. Off by default.")
    p.add_argument("--lr-pairs", nargs="+", default=None,
                    help="Explicit list of L-R pair names (var_names) to "
                         "display -- same curated-list style as the "
                         "original Fig3_heatmap.py. Overrides "
                         "--top-n-per-cluster.")
    p.add_argument("--top-n-per-cluster", type=int, default=5,
                    help="If --lr-pairs isn't given, auto-select this many "
                         "top-scoring L-R pairs per displayed cluster from "
                         "the Wilcoxon results in .uns (default: 5).")
    p.add_argument("--wilcox-key", default="louvain-wilc",
                    help="adata.uns key holding rank_genes_groups results "
                         "(default: 'louvain-wilc', written by "
                         "run_milo_da.py's --wilcox-out step).")
    p.add_argument("--cmap", default="Reds", help="Matplotlib colormap (default: 'Reds').")
    p.add_argument("--celltype-col", default="VectorType",
                    help="obs column giving each cell pair's sender-receiver "
                         "cell type identity (default: 'VectorType'). Used "
                         "to label columns by cell type instead of a bare "
                         "cluster number -- pass --no-celltype-labels to "
                         "turn this off and show raw cluster numbers.")
    p.add_argument("--no-celltype-labels", action="store_true",
                    help="Label columns with the raw cluster number/name "
                         "from --groupby instead of '<number>: <dominant "
                         "cell type pair>'.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    adata = sc.read_h5ad(args.input)

    if args.groupby not in adata.obs.columns:
        sys.exit(
            f"--groupby '{args.groupby}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    clusters = select_clusters(
        adata, args.groupby, args.clusters, args.top_n_clusters,
        args.rank_clusters_by, args.include_unassigned,
        args.celltype_filter, args.celltype_col,
        min_cluster_size=args.min_cluster_size,
    )

    idx = adata.obs[args.groupby].astype(str).isin(clusters)
    adata_highlight = adata[idx]
    print(f"Showing {len(clusters)} cluster(s): {clusters} ({adata_highlight.n_obs} cell pairs)")

    print_cluster_significance(adata, clusters, args.groupby)

    if args.lr_pairs:
        lr_pairs = args.lr_pairs
    else:
        lr_pairs = auto_select_lr_pairs(adata, clusters, args.wilcox_key, args.top_n_per_cluster)
        print(
            f"Auto-selected {len(lr_pairs)} L-R pair(s) from the top "
            f"{args.top_n_per_cluster} per cluster in "
            f"adata.uns[{args.wilcox_key!r}]: {lr_pairs}"
        )

    missing_lr = [lr for lr in lr_pairs if lr not in adata_highlight.var_names]
    if missing_lr:
        print(f"Warning: L-R pair(s) not found in var_names, dropping: {missing_lr}")
        lr_pairs = [lr for lr in lr_pairs if lr not in missing_lr]
    if not lr_pairs:
        sys.exit("No valid L-R pairs left to plot.")

    import matplotlib
    matplotlib.use("Agg")

    plot_groupby = args.groupby
    plot_categories_order = clusters
    if not args.no_celltype_labels:
        label_map = dominant_celltype_labels(adata_highlight, args.groupby, clusters, args.celltype_col)
        if any(label_map[c] != c for c in clusters):
            adata_highlight = adata_highlight.copy()
            adata_highlight.obs["_cluster_label"] = (
                adata_highlight.obs[args.groupby].astype(str).map(label_map).astype("category")
            )
            plot_groupby = "_cluster_label"
            plot_categories_order = [label_map[c] for c in clusters]
            print(f"Labeling columns by dominant '{args.celltype_col}' per cluster: {plot_categories_order}")

    fig = sc.pl.matrixplot(
        adata_highlight,
        lr_pairs,
        groupby=plot_groupby,
        dendrogram=False,
        swap_axes=True,
        categories_order=plot_categories_order,
        standard_scale="var",
        cmap=args.cmap,
        return_fig=True,
        show=False,
    )
    fig.savefig(args.output, bbox_inches="tight", dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
