#!/usr/bin/env python
"""
plot_lr_comparison_heatmap.py -- reference- vs treatment-condition differential
expression of L-R mechanism scores, per cell-pair cluster.

This answers a different question than plot_lr_heatmap.py. That script's
heatmap colors show a *cluster-vs-rest* Wilcoxon test on the raw NICHES
scores -- "which L-R pairs define this cluster's signaling identity?" It
says nothing about whether a given L-R pair's score actually changed
between conditions.

This script runs its own comparison instead, independently within each
displayed cluster: for every (L-R pair, cluster) cell, it compares NICHES
scores between --condition-col groups 1 (treatment) and 0 (reference),
among only the cell pairs assigned to that cluster. Two numbers come out of
that comparison, per cell of the heatmap:

  - Effect size (color): Cohen's d = (mean_treatment - mean_reference) /
    pooled_std. Positive = higher in the treatment condition; negative =
    higher in the reference condition. This is deliberately NOT scanpy's built-in
    "logfoldchanges" -- that assumes the input is already log1p-transformed
    (it un-logs with expm1 before taking a ratio), which NICHES scores
    aren't guaranteed to be. Cohen's d makes no such assumption and is
    scale-free, so L-R pairs with very different score magnitudes are still
    comparable on one shared color scale.

  - Significance (asterisk): a Mann-Whitney U / Wilcoxon rank-sum test
    (same test family run_milo_da.py's --wilcox-out already uses elsewhere
    in this pipeline), Benjamini-Hochberg-corrected across every (L-R pair,
    cluster) combination tested together. Cells with adjusted p < --alpha
    get a "*". Gray cells mean the comparison couldn't be run at all (too
    few cell pairs in one condition, within that cluster).

Which L-R pairs to test/show, by default, reuses plot_lr_heatmap.py's own
selection: the top --top-n-per-cluster pairs (by cluster-vs-rest Wilcoxon
score, from adata.uns[--wilcox-key]) per displayed cluster -- i.e. "for the
L-R pairs that define each cell-pair cluster's identity, how does their
signaling shift between the reference and treatment conditions within that
cluster?" Pass --lr-pairs for an
explicit list instead.

Run with `uv run` from the repository root:

    uv run python plot_lr_comparison_heatmap.py \\
        --input data/my_dataset_MILO.h5ad \\
        --output data/my_dataset_lr_comparison_heatmap.png

See PIPELINE_STEPS.md / README.md for the rest of the pipeline.
"""
import argparse
import sys

import numpy as np
import scanpy as sc

from _lib import (
    select_clusters, dominant_celltype_labels, print_cluster_significance,
    auto_select_lr_pairs, compute_comparison_stats,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, help="NICHES+Milo h5ad from run_milo_da.py.")
    p.add_argument("--output", required=True, help="Output image path (e.g. heatmap.png).")
    p.add_argument("--groupby", default="louvain_str",
                    help="obs column identifying cell-pair clusters (default: "
                         "'louvain_str', written by run_milo_da.py).")
    p.add_argument("--condition-col", default="cond_binary",
                    help="obs column with the 0 (reference)/1 (treatment) "
                         "condition encoding (default: 'cond_binary', "
                         "written by run_milo_da.py from prep_niches_input.py's "
                         "--reference-value/--treatment-value).")
    p.add_argument("--clusters", nargs="+", default=None,
                    help="Which cluster labels to display (default: all "
                         "except '-1' = unassigned, sorted numerically where "
                         "possible).")
    p.add_argument("--include-unassigned", action="store_true",
                    help="Include the '-1' (unassigned) group as a column "
                         "even when --clusters isn't given.")
    p.add_argument("--top-n-clusters", type=int, default=None,
                    help="Narrow the default cluster list to this many, "
                         "ranked by --rank-clusters-by. Ignored if --clusters is set.")
    p.add_argument("--rank-clusters-by", choices=["size", "spatialfdr"], default="size",
                    help="How to rank clusters for --top-n-clusters (default: 'size').")
    p.add_argument("--celltype-filter", default=None,
                    help="Restrict to clusters whose DOMINANT --celltype-col "
                         "value contains this text (case-insensitive), e.g. "
                         "--celltype-filter Macrophage. Applied before "
                         "--top-n-clusters, so combining the two gives the "
                         "top N *matching* clusters.")
    p.add_argument("--min-cluster-size", type=int, default=None,
                    help="Drop clusters with fewer than this many total "
                         "cell pairs before display -- a quality floor "
                         "against clusters Milo called significant on a "
                         "handful of cell pairs. Distinct from "
                         "--min-cells-per-group below, which checks each "
                         "condition group's size within an already-shown "
                         "cluster. Applied before --top-n-clusters, same as "
                         "--celltype-filter. Off by default.")
    p.add_argument("--lr-pairs", nargs="+", default=None,
                    help="Explicit list of L-R pair names (var_names) to test/display.")
    p.add_argument("--top-n-per-cluster", type=int, default=5,
                    help="If --lr-pairs isn't given, auto-select this many "
                         "top cluster-identity L-R pairs per displayed "
                         "cluster from the Wilcoxon results in .uns "
                         "(default: 5) -- same selection plot_lr_heatmap.py uses.")
    p.add_argument("--wilcox-key", default="louvain-wilc",
                    help="adata.uns key holding cluster-identity rank_genes_groups "
                         "results, used only for --lr-pairs auto-selection "
                         "(default: 'louvain-wilc').")
    p.add_argument("--alpha", type=float, default=0.05,
                    help="Significance threshold (BH-adjusted p-value) for "
                         "marking cells with an asterisk (default: 0.05). "
                         "This is this script's OWN reference-vs-treatment test threshold "
                         "-- unrelated to Step 3's Milo --alpha (SpatialFDR).")
    p.add_argument("--min-cells-per-group", type=int, default=3,
                    help="Minimum cell pairs required in each condition "
                         "group, per cluster, to run the comparison "
                         "(default: 3). Below this, that cell is left blank/gray.")
    p.add_argument("--cmap", default="RdBu_r", help="Diverging matplotlib colormap (default: 'RdBu_r').")
    p.add_argument("--celltype-col", default="VectorType",
                    help="obs column for dominant-cell-type column labels "
                         "(default: 'VectorType'); pass --no-celltype-labels to disable.")
    p.add_argument("--no-celltype-labels", action="store_true",
                    help="Label columns with raw cluster numbers instead of "
                         "'<number>: <dominant cell type pair>'.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    adata = sc.read_h5ad(args.input)

    if args.groupby not in adata.obs.columns:
        sys.exit(f"--groupby '{args.groupby}' not found in adata.obs. "
                  f"Available: {list(adata.obs.columns)}")
    if args.condition_col not in adata.obs.columns:
        sys.exit(f"--condition-col '{args.condition_col}' not found in adata.obs. "
                  f"Available: {list(adata.obs.columns)}")

    clusters = select_clusters(
        adata, args.groupby, args.clusters, args.top_n_clusters,
        args.rank_clusters_by, args.include_unassigned,
        args.celltype_filter, args.celltype_col,
        min_cluster_size=args.min_cluster_size,
    )
    print(f"Showing {len(clusters)} cluster(s): {clusters}")
    print_cluster_significance(adata, clusters, args.groupby)

    if args.lr_pairs:
        lr_pairs = args.lr_pairs
    else:
        lr_pairs = auto_select_lr_pairs(adata, clusters, args.wilcox_key, args.top_n_per_cluster)
        print(
            f"Auto-selected {len(lr_pairs)} L-R pair(s) from the top "
            f"{args.top_n_per_cluster} cluster-identity pair(s) per cluster "
            f"in adata.uns[{args.wilcox_key!r}]: {lr_pairs}"
        )
    missing_lr = [lr for lr in lr_pairs if lr not in adata.var_names]
    if missing_lr:
        print(f"Warning: L-R pair(s) not found in var_names, dropping: {missing_lr}")
        lr_pairs = [lr for lr in lr_pairs if lr not in missing_lr]
    if not lr_pairs:
        sys.exit("No valid L-R pairs left to test/plot.")

    print(
        f"\nRunning reference- vs treatment-condition comparison ({args.condition_col}: "
        f"1 vs 0) for {len(lr_pairs)} L-R pair(s) x {len(clusters)} cluster(s)..."
    )
    effect, padj = compute_comparison_stats(
        adata, clusters, args.groupby, args.condition_col, lr_pairs,
        min_cells_per_group=args.min_cells_per_group,
    )

    col_labels = clusters
    if not args.no_celltype_labels:
        label_map = dominant_celltype_labels(adata, args.groupby, clusters, args.celltype_col)
        if any(label_map[c] != c for c in clusters):
            col_labels = [label_map[c] for c in clusters]
            print(f"Labeling columns by dominant '{args.celltype_col}' per cluster: {col_labels}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy.ma as ma

    M = ma.masked_invalid(effect.values.astype(float))
    finite_vals = effect.values[np.isfinite(effect.values.astype(float))]
    vmax = np.abs(finite_vals).max() if finite_vals.size else 1.0
    vmax = vmax if vmax > 0 else 1.0

    fig_w = max(6, 0.6 * len(clusters) + 2)
    fig_h = max(4, 0.35 * len(lr_pairs) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad("lightgray")
    im = ax.imshow(M, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticks(range(len(lr_pairs)))
    ax.set_yticklabels(lr_pairs)

    for i in range(len(lr_pairs)):
        for j in range(len(clusters)):
            p = padj.values[i, j]
            if np.isfinite(p) and p < args.alpha:
                ax.text(j, i, "*", ha="center", va="center",
                         color="black", fontsize=12, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f"Cohen's d ({args.condition_col}: 1 vs 0)")
    ax.set_title(
        f"Reference- vs treatment-condition L-R score shift per cluster\n"
        f"(* = BH-adjusted p < {args.alpha}; gray = untested, "
        f"< {args.min_cells_per_group} cell pairs in a condition group)"
    )
    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight", dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
