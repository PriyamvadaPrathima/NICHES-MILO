#!/usr/bin/env python
"""
plot_lr_comparison_sidebyside.py -- two L-R heatmaps side by side, one per
condition, so you can directly compare reference- vs treatment-condition
signaling magnitude across clusters.

plot_lr_comparison_heatmap.py already answers "did this L-R pair's score
change between conditions, and is that change significant?" with a single
heatmap of the *difference* (Cohen's d). This script instead shows each
condition's own mean score as its own panel -- left panel = condition 0
(reference), right panel = condition 1 (treatment) -- so you can see the
actual signaling levels next to each other, not just the delta between them.

Both panels share:

  - The same L-R pairs (rows) and clusters (columns), in the same order.
  - The same color scale, so a given color means the same thing in both
    panels. By default (--scale row) each row is scaled by the largest
    value seen for that L-R pair across *both* panels combined -- this
    keeps different L-R pairs' very different score magnitudes readable on
    one colormap while still making the two conditions directly comparable
    (unlike plot_lr_heatmap.py's standard_scale='var', which scales each
    panel independently and would hide a pair that's simply higher
    everywhere in the treatment condition). Pass --scale none to plot raw
    NICHES scores instead, on one shared 0-to-global-max scale.

  - Significance markers on the treatment-condition (right) panel: an
    asterisk marks (L-R pair, cluster) cells where the reference-vs-treatment
    difference is significant -- the same Benjamini-Hochberg-adjusted
    Mann-Whitney U test plot_lr_comparison_heatmap.py uses (see
    _lib.compute_comparison_stats()).

Which clusters/L-R pairs to show, by default, works the same way as the
other two heatmap scripts: top --top-n-clusters clusters (ranked by
--rank-clusters-by), top --top-n-per-cluster L-R pairs per cluster (from
the cluster-identity Wilcoxon results in adata.uns[--wilcox-key]).

Run with `uv run` from the repository root:

    uv run python plot_lr_comparison_sidebyside.py \\
        --input data/my_dataset_MILO.h5ad \\
        --output data/my_dataset_lr_comparison_sidebyside.png \\
        --top-n-clusters 10

See PIPELINE_STEPS.md / README.md for the rest of the pipeline.
"""
import argparse
import sys

import numpy as np
import scanpy as sc

from _lib import (
    select_clusters, dominant_celltype_labels, print_cluster_significance,
    auto_select_lr_pairs, compute_comparison_stats, compute_group_means,
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
                         "written by run_milo_da.py).")
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
                    help="Drop clusters with fewer than this many cell "
                         "pairs before display -- a quality floor against "
                         "clusters Milo called significant on a handful of "
                         "cell pairs. Applied before --top-n-clusters, same "
                         "as --celltype-filter. Off by default.")
    p.add_argument("--lr-pairs", nargs="+", default=None,
                    help="Explicit list of L-R pair names (var_names) to display.")
    p.add_argument("--top-n-per-cluster", type=int, default=5,
                    help="If --lr-pairs isn't given, auto-select this many "
                         "top cluster-identity L-R pairs per displayed "
                         "cluster from the Wilcoxon results in .uns "
                         "(default: 5) -- same selection the other two "
                         "heatmap scripts use.")
    p.add_argument("--wilcox-key", default="louvain-wilc",
                    help="adata.uns key holding cluster-identity rank_genes_groups "
                         "results, used only for --lr-pairs auto-selection "
                         "(default: 'louvain-wilc').")
    p.add_argument("--scale", choices=["row", "none"], default="row",
                    help="'row' (default): scale each L-R pair's row by its "
                         "own max across BOTH panels combined, so the two "
                         "conditions stay comparable while different pairs' magnitudes "
                         "don't wash each other out. 'none': plot raw "
                         "NICHES scores on one shared 0-to-global-max scale.")
    p.add_argument("--alpha", type=float, default=0.05,
                    help="Significance threshold (BH-adjusted p-value) for "
                         "marking cells with an asterisk on the treatment-condition "
                         "panel (default: 0.05).")
    p.add_argument("--min-cells-per-group", type=int, default=3,
                    help="Minimum cell pairs required in each condition "
                         "group, per cluster, to compute a mean / run the "
                         "significance test (default: 3). Below this, that "
                         "cell is left blank/gray in both panels.")
    p.add_argument("--cmap", default="Reds",
                    help="Sequential matplotlib colormap (default: 'Reds') -- "
                         "these are magnitudes, not signed differences, so "
                         "unlike plot_lr_comparison_heatmap.py this isn't diverging.")
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
        sys.exit("No valid L-R pairs left to plot.")

    print(
        f"\nComputing reference/treatment mean NICHES scores for {len(lr_pairs)} L-R "
        f"pair(s) x {len(clusters)} cluster(s)..."
    )
    mean0, mean1 = compute_group_means(adata, clusters, args.groupby, args.condition_col, lr_pairs)
    _, padj = compute_comparison_stats(
        adata, clusters, args.groupby, args.condition_col, lr_pairs,
        min_cells_per_group=args.min_cells_per_group,
    )
    # compute_group_means() has no minimum-cell-count guard of its own (a
    # mean of 1 cell pair is still a valid mean) -- but a cell too small to
    # have been tested for significance is also too small to trust visually
    # next to properly-powered ones, so mask it out here the same way.
    for cl in clusters:
        mask_cl = (adata.obs[args.groupby].astype(str) == cl).values
        cond_sub = adata.obs.loc[mask_cl, args.condition_col]
        n0 = int((cond_sub == 0).sum())
        n1 = int((cond_sub == 1).sum())
        if n0 < args.min_cells_per_group:
            mean0[cl] = np.nan
        if n1 < args.min_cells_per_group:
            mean1[cl] = np.nan

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

    raw0 = mean0.values.astype(float)
    raw1 = mean1.values.astype(float)

    combined = np.hstack([raw0, raw1])

    if args.scale == "row":
        row_max = np.nanmax(combined, axis=1)
        row_max = np.where((row_max > 0) & np.isfinite(row_max), row_max, np.nan)
        M0 = raw0 / row_max[:, None]
        M1 = raw1 / row_max[:, None]
        vmin, vmax = 0.0, 1.0
        cbar_label = "NICHES score (row-scaled, shared reference/treatment max)"
    else:
        M0, M1 = raw0, raw1
        combined_max = np.nanmax(combined) if np.isfinite(combined).any() else 1.0
        vmin, vmax = 0.0, (combined_max if combined_max > 0 else 1.0)
        cbar_label = "Mean NICHES score (raw)"

    M0 = ma.masked_invalid(M0)
    M1 = ma.masked_invalid(M1)

    n_clusters, n_pairs = len(clusters), len(lr_pairs)
    fig_w = max(9, 1.1 * n_clusters + 4)
    fig_h = max(4, 0.35 * n_pairs + 2)
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h), sharey=True)
    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad("lightgray")

    im0 = axes[0].imshow(M0, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    im1 = axes[1].imshow(M1, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")

    axes[0].set_title("Condition 0 (reference)")
    axes[1].set_title("Condition 1 (treatment)")
    for ax in axes:
        ax.set_xticks(range(n_clusters))
        ax.set_xticklabels(col_labels, rotation=45, ha="right")
    axes[0].set_yticks(range(n_pairs))
    axes[0].set_yticklabels(lr_pairs)

    for i in range(n_pairs):
        for j in range(n_clusters):
            p = padj.values[i, j]
            if np.isfinite(p) and p < args.alpha:
                axes[1].text(j, i, "*", ha="center", va="center",
                             color="black", fontsize=12, fontweight="bold")

    fig.colorbar(im1, ax=axes, shrink=0.8, label=cbar_label)
    fig.suptitle(
        f"Reference- vs treatment-condition L-R signaling magnitude per cluster\n"
        f"(* on treatment panel = BH-adjusted p < {args.alpha} vs. reference; "
        f"gray = < {args.min_cells_per_group} cell pairs in that group)"
    )
    fig.savefig(args.output, bbox_inches="tight", dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
