#!/usr/bin/env python
"""
plot_lr_log2fc_heatmap.py -- log2(treatment/reference) L-R fold-change
heatmap, generalized from Bridges et al.'s figures/human_cd40ag_validation.py
(its Fig. 7G human validation heatmap; see milo_helpers.py's header for
attribution).

That original script is a standalone, curated-list analysis on an
independent human dataset (GSE244739): for 15 hand-picked L-R pairs and 3
hand-picked sender->receiver interactions, it computes
score = mean(ligand) * mean(receptor) separately in the pre- and
post-sotigalimab timepoints, then plots log2(post_score / pre_score) as a
diverging heatmap -- red = higher post-treatment, blue = higher
pre-treatment.

This version reuses that exact visual language (RdBu_r, TwoSlopeNorm
centered at 0, white gridlines, log2FC values printed in each cell) but
computes the underlying scores the way the rest of this pipeline does:
mean NICHES L-R score per (L-R pair, cluster) in each condition
(_lib.compute_group_means()), rather than a hand-picked sender/receiver
score computed from scratch. Columns are your NICHES+Milo sc_louvain
clusters (or a --celltype-filter/--top-n-clusters subset of them), not a
curated interaction list -- so this generalizes to any dataset/cell-type
combination instead of the original's fixed 3 interactions.

One addition beyond the original: an optional (--no-significance to turn
off) asterisk overlay using the same Benjamini-Hochberg-adjusted
Mann-Whitney U test plot_lr_comparison_heatmap.py uses -- the original figure
is purely descriptive (a ratio of two means, no test attached).

Run with `uv run` from the repository root:

    uv run python plot_lr_log2fc_heatmap.py \\
        --input data/my_dataset_MILO.h5ad \\
        --output data/my_dataset_lr_log2fc_heatmap.png \\
        --celltype-filter Macrophage

See PIPELINE_STEPS.md / README.md for the rest of the pipeline.
"""
import argparse
import sys

import numpy as np
import scanpy as sc

from _lib import (
    select_clusters, dominant_celltype_labels, print_cluster_significance,
    auto_select_lr_pairs, compute_group_means, compute_log2fc, compute_comparison_stats,
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
                         "--celltype-filter Macrophage for every macrophage- "
                         "involving cluster. Applied before --top-n-clusters.")
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
                         "(default: 5) -- same selection the other heatmap "
                         "scripts use.")
    p.add_argument("--wilcox-key", default="louvain-wilc",
                    help="adata.uns key holding cluster-identity rank_genes_groups "
                         "results, used only for --lr-pairs auto-selection "
                         "(default: 'louvain-wilc').")
    p.add_argument("--vmax", type=float, default=3.0,
                    help="Color scale runs -vmax to +vmax log2 units "
                         "(default: 3.0, same as the original Fig. 7G) -- "
                         "values beyond this just show fully saturated.")
    p.add_argument("--cap", type=float, default=4.0,
                    help="log2FC value substituted for a near-zero-baseline "
                         "pair (--zero-thresh) that can't get a literal ratio "
                         "(default: 4.0, same as the original). Deliberately "
                         "larger than --vmax, so these 'signal appeared/"
                         "vanished from nothing' cells always render fully "
                         "saturated rather than blending in with an ordinary "
                         "large-but-computable fold-change.")
    p.add_argument("--zero-thresh", type=float, default=1e-3,
                    help="Mean NICHES scores below this are treated as "
                         "'no signal' when computing log2FC, to avoid a "
                         "near-zero denominator blowing up the ratio "
                         "(default: 0.001, same as the original). See "
                         "_lib.compute_log2fc()'s docstring for exactly how "
                         "near-zero baselines are handled.")
    p.add_argument("--min-cells-for-mean", type=int, default=1,
                    help="Minimum cell pairs required in a condition group, "
                         "per cluster, to compute a mean at all (default: 1 "
                         "-- blank only when a group is completely empty).")
    p.add_argument("--no-significance", action="store_true",
                    help="Skip the Benjamini-Hochberg-adjusted Mann-Whitney "
                         "significance overlay (this pipeline's own addition "
                         "-- the original Fig. 7G is purely descriptive, no "
                         "test attached). Significant cells get a black dot.")
    p.add_argument("--alpha", type=float, default=0.05,
                    help="Significance threshold (BH-adjusted p-value) for "
                         "the dot marker, if not --no-significance (default: 0.05).")
    p.add_argument("--min-cells-per-group", type=int, default=3,
                    help="Minimum cell pairs required in each condition "
                         "group, per cluster, to run the significance test "
                         "(default: 3) -- separate from --min-cells-for-mean, "
                         "which only gates whether a log2FC value is shown "
                         "at all.")
    p.add_argument("--celltype-col", default="VectorType",
                    help="obs column for dominant-cell-type column labels "
                         "(default: 'VectorType'); pass --no-celltype-labels to disable.")
    p.add_argument("--no-celltype-labels", action="store_true",
                    help="Label columns with raw cluster numbers instead of "
                         "'<number>: <dominant cell type pair>'.")
    return p.parse_args(argv)


def _italicize_lr_label(label):
    """'Ligand—Receptor' -> italic gene names either side of the en-dash,
    matching the original figure's y-axis label style. Falls back to the
    plain label if it doesn't split cleanly on an em/en-dash.
    """
    for dash in ("—", "–", "-"):
        if dash in label:
            parts = label.split(dash)
            if len(parts) == 2:
                return f"$\\it{{{parts[0]}}}${dash}$\\it{{{parts[1]}}}$"
    return label


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
        f"\nComputing log2(treatment/reference) NICHES score fold-change for "
        f"{len(lr_pairs)} L-R pair(s) x {len(clusters)} cluster(s)..."
    )
    mean0, mean1 = compute_group_means(adata, clusters, args.groupby, args.condition_col, lr_pairs)
    for cl in clusters:
        mask_cl = (adata.obs[args.groupby].astype(str) == cl).values
        cond_sub = adata.obs.loc[mask_cl, args.condition_col]
        n0 = int((cond_sub == 0).sum())
        n1 = int((cond_sub == 1).sum())
        if n0 < args.min_cells_for_mean:
            mean0[cl] = np.nan
        if n1 < args.min_cells_for_mean:
            mean1[cl] = np.nan

    fc = compute_log2fc(mean0, mean1, cap=args.cap, zero_thresh=args.zero_thresh)

    padj = None
    if not args.no_significance:
        _, padj = compute_comparison_stats(
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
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import numpy.ma as ma

    n_rows, n_cols = len(lr_pairs), len(clusters)
    data = ma.masked_invalid(fc.values.astype(float))
    norm = mcolors.TwoSlopeNorm(vmin=-args.vmax, vcenter=0, vmax=args.vmax)

    fig_w = max(4.5, 0.75 * n_cols + 2.5)
    fig_h = max(3.5, 0.35 * n_rows + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("lightgray")
    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

    # White gridlines between cells, matching the original figure.
    for i in range(n_rows + 1):
        ax.axhline(i - 0.5, color="white", linewidth=0.5)
    for j in range(n_cols + 1):
        ax.axvline(j - 0.5, color="white", linewidth=0.5)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([_italicize_lr_label(lr) for lr in lr_pairs], fontsize=8)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, fontsize=8, rotation=30, ha="left")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    ax.set_title(
        "L–R interaction scores: treatment vs. reference condition\n"
        f"({args.condition_col}: 1 vs 0)",
        fontsize=10, pad=30,
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("log$_2$(treatment/reference)", fontsize=9)

    values = fc.values
    for i in range(n_rows):
        for j in range(n_cols):
            val = values[i, j]
            if not np.isfinite(val):
                continue
            if abs(val) > 0.01:
                text_color = "white" if abs(val) > 2.0 else "black"
                ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                        fontsize=6.5, color=text_color)
            if padj is not None:
                p = padj.values[i, j]
                if np.isfinite(p) and p < args.alpha:
                    ax.text(j, i + 0.33, "•", ha="center", va="center",
                            fontsize=9, color="black", fontweight="bold")

    subtitle = "value = log2(treatment/reference)"
    if padj is not None:
        subtitle += f"; • = BH-adjusted p < {args.alpha}"
    fig.text(0.5, 0.005, subtitle, ha="center", fontsize=7.5)

    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight", dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
