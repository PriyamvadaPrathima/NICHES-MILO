#!/usr/bin/env python
"""
plot_milo_volcano.py -- volcano plot of Milo's own differential-abundance
test: every neighborhood's logFC vs. -log10(SpatialFDR), colored by
sc_louvain cluster.

None of this pipeline's existing figures actually plot the Milo SpatialFDR
values themselves. `run_milo_da.py --plot-out` (via the bundled
`milo_helpers.plot_nhood_clusters()`) only uses SpatialFDR as an on/off
mask -- neighborhoods above --alpha are grayed out on a graph embedding,
but the SpatialFDR value itself never appears on an axis. The L-R heatmap
scripts print a per-cluster summary table (mean/max SpatialFDR) but don't
plot the underlying per-neighborhood distribution either. This script
plots that distribution directly, one point per neighborhood -- the same
"volcano plot" convention used for any differential test (Milo's own
tutorials plot DA results this way).

Reading it:

  - x-axis: logFC, Milo's per-neighborhood log fold-change of cell-pair
    abundance between your two conditions. Positive = more abundant in
    condition 1 (treatment), as encoded by prep_niches_input.py's
    --treatment-value/--reference-value. This is Milo's own statistic --
    unrelated to the L-R heatmaps' Cohen's d, which is a completely
    different test on completely different data (NICHES scores, not
    neighborhood abundance).

  - y-axis: -log10(SpatialFDR). Higher = more significant. The dashed
    horizontal line marks --alpha; points above it are the "significant"
    neighborhoods that get grouped into clusters at all (SpatialFDR is
    Milo's neighborhood-overlap-corrected FDR, computed once in Step 3 and
    read as-is here, not recomputed).

  - Point size: neighborhood size (how many cell pairs are in it), for
    HIGHLIGHTED neighborhoods only (--point-size-scale, capped at
    --max-point-size). Background neighborhoods are drawn as small,
    uniform, low-opacity dots (--background-size / --background-alpha)
    rather than sized individually -- with thousands of neighborhoods
    packed into a narrow band (Milo's SpatialFDR tracks logFC tightly),
    sizing every one of them just produces a solid gray mass with no
    visible structure; low-alpha fixed dots blend into a readable density
    gradient instead.

  - Color: which sc_louvain cluster the neighborhood belongs to (only for
    the clusters you're displaying -- --top-n-clusters etc., same
    selection as the other heatmap scripts). Everything else (clusters not
    displayed, plus every neighborhood that never cleared Step 3's
    --alpha) is plotted as light gray background for context.

Note --alpha here only draws the reference line -- it does NOT recompute
which neighborhoods belong to which cluster (that grouping is fixed by
whatever --alpha you passed to run_milo_da.py in Step 3). Pass the same
value here to see it reflected in the plot; a different value just moves
the dashed line without changing which points are colored.

Run with `uv run` from the repository root:

    uv run python plot_milo_volcano.py \\
        --input data/my_dataset_MILO.h5ad \\
        --output data/my_dataset_milo_volcano.png \\
        --top-n-clusters 10 --alpha 0.2

See PIPELINE_STEPS.md / README.md for the rest of the pipeline.
"""
import argparse
import sys

import numpy as np
import pandas as pd
import scanpy as sc

from _lib import select_clusters, dominant_celltype_labels, print_cluster_significance


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, help="NICHES+Milo h5ad from run_milo_da.py.")
    p.add_argument("--output", required=True, help="Output image path (e.g. volcano.png).")
    p.add_argument("--groupby", default="louvain_str",
                    help="obs column identifying cell-pair clusters, used "
                         "only to pick which clusters to highlight (default: "
                         "'louvain_str', written by run_milo_da.py).")
    p.add_argument("--clusters", nargs="+", default=None,
                    help="Which cluster labels to highlight in color "
                         "(default: all except '-1' = unassigned, sorted "
                         "numerically where possible).")
    p.add_argument("--include-unassigned", action="store_true",
                    help="Include the '-1' (unassigned) group in the "
                         "highlighted/colored set instead of only the gray "
                         "background.")
    p.add_argument("--top-n-clusters", type=int, default=None,
                    help="Narrow the highlighted cluster list to this many, "
                         "ranked by --rank-clusters-by. Ignored if --clusters is set.")
    p.add_argument("--rank-clusters-by", choices=["size", "spatialfdr"], default="size",
                    help="How to rank clusters for --top-n-clusters (default: 'size').")
    p.add_argument("--celltype-filter", default=None,
                    help="Restrict the highlighted set to clusters whose "
                         "DOMINANT --celltype-col value contains this text "
                         "(case-insensitive), e.g. --celltype-filter "
                         "Macrophage. Applied before --top-n-clusters, so "
                         "combining the two gives the top N *matching* "
                         "clusters.")
    p.add_argument("--min-cluster-size", type=int, default=None,
                    help="Drop clusters with fewer than this many cell "
                         "pairs from the highlighted set -- a quality floor "
                         "against clusters Milo called significant on a "
                         "handful of cell pairs. Applied before "
                         "--top-n-clusters, same as --celltype-filter. Off "
                         "by default.")
    p.add_argument("--alpha", type=float, default=0.1,
                    help="SpatialFDR threshold to draw as a reference line "
                         "(default: 0.1). Should normally match the --alpha "
                         "you passed to run_milo_da.py's Step 3 -- see the "
                         "module docstring for why this doesn't recompute "
                         "cluster membership.")
    p.add_argument("--point-size-scale", type=float, default=3.0,
                    help="Multiplier applied to each highlighted neighborhood's "
                         "Nhood_size for its marker size (default: 3.0). Only "
                         "affects the colored/highlighted clusters -- "
                         "everything else is drawn as small fixed-size dots "
                         "(see --background-size), since with thousands of "
                         "neighborhoods scaling every one by size just "
                         "produces an unreadable mass of overlapping circles.")
    p.add_argument("--max-point-size", type=float, default=250.0,
                    help="Cap on highlighted marker area (matplotlib 's' "
                         "units), so one unusually large neighborhood can't "
                         "dominate the plot (default: 250.0).")
    p.add_argument("--background-size", type=float, default=5.0,
                    help="Fixed marker size for the non-highlighted "
                         "background neighborhoods (default: 5.0).")
    p.add_argument("--background-alpha", type=float, default=0.2,
                    help="Opacity of background points (default: 0.2) -- low "
                         "opacity lets overlapping points read as a density "
                         "cloud instead of a solid gray mass.")
    p.add_argument("--cmap", default=None,
                    help="Qualitative matplotlib colormap for cluster colors. "
                         "Default: 'tab10' for up to 10 highlighted clusters, "
                         "'tab20' beyond that -- tab10's colors are all "
                         "distinct hues, while tab20 pairs a dark/light shade "
                         "per hue, which reads as visually confusing once you "
                         "have more than ~5-6 similar-hued pairs on screen.")
    p.add_argument("--celltype-col", default="VectorType",
                    help="obs column for dominant-cell-type legend labels "
                         "(default: 'VectorType'); pass --no-celltype-labels to disable.")
    p.add_argument("--no-celltype-labels", action="store_true",
                    help="Label the legend with raw cluster numbers instead "
                         "of '<number>: <dominant cell type pair>'.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    adata = sc.read_h5ad(args.input)

    if "nhood_adata" not in adata.uns:
        sys.exit(
            "adata.uns['nhood_adata'] not found -- this must be the direct "
            "h5ad output of run_milo_da.py (not a re-saved/subsetted copy)."
        )
    nhood_obs_full = adata.uns["nhood_adata"].obs
    missing_cols = [c for c in ("logFC", "SpatialFDR", "Nhood_size", "louvain") if c not in nhood_obs_full.columns]
    if missing_cols:
        sys.exit(
            f"adata.uns['nhood_adata'].obs is missing column(s) {missing_cols} "
            "-- this must be the direct h5ad output of run_milo_da.py."
        )
    if args.groupby not in adata.obs.columns:
        sys.exit(f"--groupby '{args.groupby}' not found in adata.obs. "
                  f"Available: {list(adata.obs.columns)}")

    clusters = select_clusters(
        adata, args.groupby, args.clusters, args.top_n_clusters,
        args.rank_clusters_by, args.include_unassigned,
        args.celltype_filter, args.celltype_col,
        min_cluster_size=args.min_cluster_size,
    )
    print(f"Highlighting {len(clusters)} cluster(s): {clusters}")
    print_cluster_significance(adata, clusters, args.groupby)

    nhood_obs = nhood_obs_full.copy()
    # A handful of neighborhoods can end up with NaN logFC/SpatialFDR (Milo
    # drops zero-count neighborhoods/samples during DA_nhoods_fixed's own
    # filtering, but the 'louvain' partition covers every neighborhood in
    # the original graph) -- drop those before plotting rather than letting
    # matplotlib silently skip NaN points.
    n_total = len(nhood_obs)
    nhood_obs = nhood_obs[np.isfinite(nhood_obs["logFC"]) & np.isfinite(nhood_obs["SpatialFDR"])]
    if len(nhood_obs) < n_total:
        print(f"Note: dropping {n_total - len(nhood_obs)} neighborhood(s) with no logFC/SpatialFDR "
              "(filtered out during Step 3's DA test) before plotting.")

    # THE FIX: np.where() evaluates both branches eagerly, so
    # `nhood_obs["louvain"].astype(int)` was being computed over the WHOLE
    # column -- including the `inf` entries used to mark non-significant
    # neighborhoods -- before np.where() ever got to select between them.
    # `inf` has no integer representation, so that raised
    # `IntCastingNaNError` regardless of which branch np.where() would have
    # picked. Filter to the finite rows first (same pattern _lib.py's
    # rank_clusters()/print_cluster_significance() already use), then only
    # cast those; everything else -- including any real NaN from the
    # earlier finite-logFC/SpatialFDR filter, which can't reach here since
    # louvain is never NaN, only inf or a real cluster id -- gets "-1".
    finite_louvain = nhood_obs["louvain"] < float("inf")
    cluster_str = pd.Series("-1", index=nhood_obs.index)
    cluster_str.loc[finite_louvain] = (
        nhood_obs.loc[finite_louvain, "louvain"].astype(int).astype(str)
    )
    nhood_obs["cluster_str"] = cluster_str
    nhood_obs["neglog10_fdr"] = -np.log10(nhood_obs["SpatialFDR"].clip(lower=1e-300))

    label_map = (
        dominant_celltype_labels(adata, args.groupby, clusters, args.celltype_col)
        if not args.no_celltype_labels else {c: c for c in clusters}
    )

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 6.5))

    # THE FIX (visual quality): with thousands of neighborhoods packed into
    # a narrow band -- Milo's SpatialFDR tracks logFC very tightly here, so
    # the background isn't a diffuse cloud, it's a dense curve -- sizing
    # every background point by Nhood_size and drawing it opaque just piles
    # up into a solid gray mass with no visible structure. A small FIXED
    # size plus low alpha lets overlapping points blend into a readable
    # density gradient instead (denser regions naturally look darker).
    background = nhood_obs[~nhood_obs["cluster_str"].isin(clusters)]
    ax.scatter(
        background["logFC"], background["neglog10_fdr"],
        s=args.background_size, c="lightgray", linewidths=0,
        alpha=args.background_alpha, label="other / not shown",
    )

    # tab20 pairs a dark/light shade per hue (e.g. blue, light-blue, orange,
    # light-orange, ...) -- with more than ~5 clusters on screen those pairs
    # become hard to tell apart. tab10 is 10 maximally distinct hues, a
    # better default whenever there are few enough clusters to fit it.
    cmap_name = args.cmap or ("tab10" if len(clusters) <= 10 else "tab20")
    cmap = plt.get_cmap(cmap_name)
    for i, cl in enumerate(clusters):
        sub = nhood_obs[nhood_obs["cluster_str"] == cl]
        if len(sub) == 0:
            continue
        sizes = (sub["Nhood_size"] * args.point_size_scale).clip(lower=6, upper=args.max_point_size)
        ax.scatter(
            sub["logFC"], sub["neglog10_fdr"],
            s=sizes, c=[cmap(i % cmap.N)], edgecolors="white", linewidths=0.5,
            alpha=0.85, label=label_map.get(cl, cl), zorder=3,
        )

    ax.axhline(-np.log10(args.alpha), linestyle="--", color="black", linewidth=1)
    ax.text(
        ax.get_xlim()[1], -np.log10(args.alpha), f"  alpha={args.alpha}",
        va="center", ha="left", fontsize=8,
    )
    ax.axvline(0, linestyle=":", color="gray", linewidth=0.8)

    ax.set_xlabel("logFC (condition 1 / treatment vs. condition 0 / reference)")
    ax.set_ylabel("-log10(SpatialFDR)")
    ax.set_title(
        f"Milo differential abundance per neighborhood ({len(nhood_obs)} "
        f"neighborhoods; point size = Nhood_size)"
    )
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, title="cluster", markerscale=1.5)
    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight", dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
