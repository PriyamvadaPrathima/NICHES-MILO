#!/usr/bin/env python
"""
plot_embeddings.py -- quick-look PCA/UMAP visualization for any h5ad at any
stage of the generalized pipeline (pre-NICHES cell-level data, or
post-NICHES cell-pair-level NICHES output).

Coordinates can come from three places, tried in this order unless a flag
overrides the choice:

  1. Already computed: .obsm['X_pca'] / .obsm['X_umap']. This is what you
     get for free on *_NICHES.h5ad -- run_niches_generic.R computes PCA and
     UMAP on the merged cell-pair object before writing it out.

  2. Two plain .obs columns holding precomputed 2D coordinates from an
     external source -- e.g. a published dataset's own 'UMAP1'/'UMAP2'
     columns, which survive Step 1's subsetting untouched since
     prep_niches_input.py only subsets rows and adds its own bookkeeping
     columns. Pass --umap-obs-cols UMAP1 UMAP2.

  3. Nothing precomputed: pass --compute to run scanpy's
     pp.pca + pp.neighbors + tl.umap fresh (slower, but works on any h5ad
     regardless of what stage it's from).

Run with `uv run` from the repository root:

    uv run python plot_embeddings.py \\
        --input data/my_dataset_LIM.h5ad \\
        --output data/my_dataset_LIM_embeddings.png \\
        --umap-obs-cols UMAP1 UMAP2

    uv run python plot_embeddings.py \\
        --input data/my_dataset_NICHES.h5ad \\
        --output data/my_dataset_NICHES_embeddings.png

High-cardinality columns (e.g. NICHES' own 'VectorType', which can have
dozens of sender-receiver pairs) overflow a single shared legend -- the
default grid above just drops the legend and keeps the colors
(--max-legend-categories). Pass --split-by instead to get one labeled PNG
per category: full embedding in light gray, that one category highlighted,
titled with its own name and cell-pair count -- no legend needed since each
file only ever shows one label.

    uv run python plot_embeddings.py \\
        --input data/my_dataset_NICHES.h5ad \\
        --output data/my_dataset_NICHES_embeddings.png \\
        --split-by VectorType

See PIPELINE_STEPS.md / README.md for the rest of the pipeline.
"""
import argparse
import os
import re
import sys

import scanpy as sc


_GUESS_CELLTYPE_COLS = [
    "celltype", "cell_type", "_grouping", "cluster", "VectorType", "celltype.Joint",
]
_GUESS_CONDITION_COLS = [
    "Condition", "_condition_binary", "condition", "treatment", "cond_binary",
]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, help="h5ad to plot.")
    p.add_argument("--output", required=True, help="Output image path (e.g. plot.png).")
    p.add_argument("--color-by", nargs="+", default=None,
                    help="obs column(s) to color points by, one column per "
                         "panel-column. Default: guesses one cell-type-like "
                         "and one condition-like column from common names "
                         "present in the file. Ignored when --split-by is "
                         "given -- use --split-color-by there instead.")
    p.add_argument("--umap-obs-cols", nargs=2, default=None, metavar=("XCOL", "YCOL"),
                    help="Use two existing plain .obs columns as 2D "
                         "coordinates instead of .obsm['X_umap'] (e.g. a "
                         "published dataset's own 'UMAP1' 'UMAP2' columns).")
    p.add_argument("--compute", action="store_true",
                    help="Compute PCA/neighbors/UMAP fresh with scanpy "
                         "defaults, ignoring anything already in .obsm or "
                         "--umap-obs-cols.")
    p.add_argument("--n-pcs", type=int, default=30, help="PCs for --compute (default: 30).")
    p.add_argument("--point-size", type=float, default=None, help="Passed to scanpy's size=.")
    p.add_argument("--title", default=None, help="Figure suptitle (default: --input's path).")
    p.add_argument("--max-legend-categories", type=int, default=15,
                    help="If a --color-by column has more unique values than "
                         "this, its legend is dropped instead of drawn "
                         "(default: 15). High-cardinality legends -- e.g. "
                         "'VectorType' on NICHES output can have dozens of "
                         "sender-receiver pairs, 'Condition' one entry per "
                         "split unit -- otherwise overflow their panel and "
                         "overlap the next one. The colors are still shown; "
                         "only the label key is dropped. Pass a larger "
                         "value, or 0 to never drop a legend.")
    p.add_argument("--split-by", default=None,
                    help="obs column to facet into one PNG per category "
                         "instead of a single shared-legend panel -- e.g. "
                         "--split-by VectorType. Each file shows the full "
                         "embedding in light gray with just that one "
                         "category highlighted, titled with the category "
                         "name and its cell-pair count, so every file is "
                         "self-labeled without needing a legend at all. "
                         "When given, --color-by/--max-legend-categories are "
                         "ignored and --output becomes a filename prefix: "
                         "one file per category is written alongside it as "
                         "'<stem>__<column>-<category><ext>'. Pass "
                         "--split-color-by to color the highlighted points "
                         "themselves by a second column instead of one flat "
                         "--highlight-color.")
    p.add_argument("--split-values", nargs="+", default=None,
                    help="If --split-by is given, only write files for "
                         "these specific category values instead of every "
                         "unique value present.")
    p.add_argument("--split-basis", choices=["auto", "umap", "pca"], default="auto",
                    help="Which single 2D embedding to use for --split-by "
                         "output (default: 'auto' -- prefers precomputed "
                         "UMAP, falls back to PCA; --umap-obs-cols/--compute "
                         "still override this the same way they do in the "
                         "normal multi-panel mode).")
    p.add_argument("--split-color-by", default=None,
                    help="obs column to color the HIGHLIGHTED points by in "
                         "--split-by output, instead of a single flat "
                         "--highlight-color -- e.g. --split-by VectorType "
                         "--split-color-by Condition to see, for each "
                         "sender-receiver pair, which patients/timepoints "
                         "its cell pairs actually come from. Colors are "
                         "assigned from this column's FULL set of values "
                         "across the whole dataset, so the same value (e.g. "
                         "the same patient) gets the same color in every "
                         "output file, not just within one. A small legend "
                         "is added to each file. Background (non-highlighted) "
                         "points stay gray either way.")
    p.add_argument("--min-count", type=int, default=1,
                    help="Skip --split-by categories with fewer than this "
                         "many cell pairs (default: 1 -- skip only "
                         "genuinely empty categories).")
    p.add_argument("--highlight-color", default="crimson",
                    help="Highlight color for the --split-by category in "
                         "each file (default: 'crimson'). Ignored if "
                         "--split-color-by is given.")
    return p.parse_args(argv)


def _slugify(value):
    """Filesystem-safe filename fragment from a category value. NICHES'
    own VectorType values contain em-dashes and spaces (e.g.
    'CD8_mem_T_cells—Macrophages'), which aren't all safe/readable in
    filenames -- this maps the sender-receiver em-dash to a plain '-to-'
    and collapses everything else non-alphanumeric to underscores.
    """
    s = str(value).replace("—", "-to-").replace("–", "-to-")
    s = re.sub(r"[^\w\-.]+", "_", s).strip("_")
    return s or "unnamed"


def resolve_split_basis(adata, args):
    """Resolve a single 2D embedding (already or about-to-be stored in
    .obsm) for --split-by output, following the same precedence as the
    normal multi-panel mode's coordinate resolution, but returning only one
    winner instead of every basis available (one shared embedding is used
    across every --split-by output file, so they stay visually comparable).
    """
    if args.compute:
        print("Computing PCA + neighbors + UMAP fresh (--compute)...")
        sc.pp.pca(adata, n_comps=args.n_pcs)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        return "umap", "UMAP (computed)"
    if args.umap_obs_cols:
        xcol, ycol = args.umap_obs_cols
        missing_uc = [c for c in (xcol, ycol) if c not in adata.obs.columns]
        if missing_uc:
            sys.exit(
                f"--umap-obs-cols column(s) not found: {missing_uc}. "
                f"Available: {list(adata.obs.columns)}"
            )
        adata.obsm["X_umap_obscols"] = adata.obs[[xcol, ycol]].to_numpy(dtype=float)
        return "umap_obscols", f"UMAP ({xcol}/{ycol})"
    if args.split_basis == "pca":
        if "X_pca" not in adata.obsm:
            sys.exit("--split-basis pca requested but .obsm['X_pca'] not found.")
        return "pca", "PCA (precomputed)"
    if args.split_basis in ("auto", "umap"):
        if "X_umap" in adata.obsm:
            return "umap", "UMAP (precomputed)"
        if args.split_basis == "umap":
            sys.exit("--split-basis umap requested but .obsm['X_umap'] not found.")
    if "X_pca" in adata.obsm:
        return "pca", "PCA (precomputed)"
    sys.exit(
        "No PCA/UMAP found in .obsm and no --umap-obs-cols given. Pass "
        "--compute to generate them fresh, or --umap-obs-cols if this file "
        "has 2D coordinates stored as plain obs columns. Available obs "
        f"columns: {list(adata.obs.columns)}"
    )


def do_split_plots(adata, args):
    if args.split_by not in adata.obs.columns:
        sys.exit(
            f"--split-by '{args.split_by}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )
    if args.split_color_by and args.split_color_by not in adata.obs.columns:
        sys.exit(
            f"--split-color-by '{args.split_color_by}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    basis, basis_title = resolve_split_basis(adata, args)
    coords = adata.obsm[f"X_{basis}"][:, :2]

    values = adata.obs[args.split_by].astype(str)
    counts = values.value_counts()

    if args.split_values:
        cats = [str(v) for v in args.split_values]
        missing_vals = [v for v in cats if v not in counts.index]
        if missing_vals:
            print(f"Warning: --split-values not present, skipping: {missing_vals}")
        cats = [v for v in cats if v in counts.index]
    else:
        cats = sorted(counts.index.tolist())

    cats = [c for c in cats if counts.get(c, 0) >= args.min_count]
    if not cats:
        sys.exit(
            f"No '{args.split_by}' categories with >= --min-count "
            f"{args.min_count} cell pairs to plot."
        )

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # --split-color-by: build ONE stable value->color map from the column's
    # full set of values across the whole dataset (not just what's present
    # within a single --split-by category), so e.g. a given patient is the
    # same color in every output file, not just internally consistent
    # within one.
    sub_color_map = None
    sub_cats_all = None
    if args.split_color_by:
        sub_series_all = adata.obs[args.split_color_by].astype(str)
        sub_cats_all = sorted(sub_series_all.unique().tolist())
        cmap_name = "tab10" if len(sub_cats_all) <= 10 else "tab20"
        cmap = plt.get_cmap(cmap_name)
        sub_color_map = {c: cmap(i % cmap.N) for i, c in enumerate(sub_cats_all)}
        print(
            f"--split-color-by {args.split_color_by!r}: {len(sub_cats_all)} "
            f"value(s), colors held constant across every output file: "
            f"{sub_cats_all}"
        )

    out_dir = os.path.dirname(args.output) or "."
    stem, ext = os.path.splitext(os.path.basename(args.output))
    ext = ext or ".png"
    os.makedirs(out_dir, exist_ok=True)

    # Same xlim/ylim on every file (rather than each panel auto-scaling to
    # just its own highlighted points) so the files stay visually
    # comparable -- a category's position within the whole embedding is
    # part of what you're trying to read off these.
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    xpad = 0.05 * ((x_max - x_min) or 1)
    ypad = 0.05 * ((y_max - y_min) or 1)

    written = []
    for cat in cats:
        mask = (values == cat).values
        n = int(mask.sum())
        fig, ax = plt.subplots(figsize=(7.4, 5.5) if sub_color_map else (6, 5.5))
        # Same density-cloud background trick used in plot_milo_volcano.py:
        # small fixed size, low alpha, so overlapping background points read
        # as a shape instead of a solid mass.
        ax.scatter(
            coords[~mask, 0], coords[~mask, 1],
            s=4, c="lightgray", linewidths=0, alpha=0.25, zorder=1,
        )
        if sub_color_map:
            sub_cat_vals = adata.obs.loc[mask, args.split_color_by].astype(str).values
            # Loop per sub-category (rather than one scatter with a color
            # array) so each gets its own legend handle/label.
            for subcat in sub_cats_all:
                subcat_mask = sub_cat_vals == subcat
                if not subcat_mask.any():
                    continue
                pts = coords[mask][subcat_mask]
                ax.scatter(
                    pts[:, 0], pts[:, 1],
                    s=12, c=[sub_color_map[subcat]], edgecolors="white",
                    linewidths=0.3, alpha=0.85, zorder=2, label=subcat,
                )
            ax.legend(
                bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7,
                title=args.split_color_by, markerscale=1.3, frameon=False,
            )
        else:
            ax.scatter(
                coords[mask, 0], coords[mask, 1],
                s=12, c=args.highlight_color, edgecolors="white", linewidths=0.3,
                alpha=0.85, zorder=2,
            )
        ax.set_xlim(x_min - xpad, x_max + xpad)
        ax.set_ylim(y_min - ypad, y_max + ypad)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{cat}\n({n} cell pair{'s' if n != 1 else ''})", fontsize=10)
        fig.suptitle(f"{args.split_by} -- {basis_title}", fontsize=9, y=0.98)
        fig.tight_layout()

        out_path = os.path.join(
            out_dir, f"{stem}__{_slugify(args.split_by)}-{_slugify(cat)}{ext}"
        )
        fig.savefig(out_path, bbox_inches="tight", dpi=150)
        plt.close(fig)
        written.append(out_path)
        print(f"Wrote {out_path} ({n} cell pairs)")

    print(f"\nWrote {len(written)} file(s) split by '{args.split_by}' into {out_dir}/")


def guess_color_cols(adata):
    cols = []
    for candidates in (_GUESS_CELLTYPE_COLS, _GUESS_CONDITION_COLS):
        for c in candidates:
            if c in adata.obs.columns:
                cols.append(c)
                break
    if not cols:
        sys.exit(
            "Could not guess any columns to color by. Pass --color-by "
            f"explicitly -- available obs columns: {list(adata.obs.columns)}"
        )
    return cols


def main(argv=None):
    args = parse_args(argv)
    adata = sc.read_h5ad(args.input)

    if args.split_by:
        do_split_plots(adata, args)
        return

    color_by = args.color_by or guess_color_cols(adata)
    missing = [c for c in color_by if c not in adata.obs.columns]
    if missing:
        sys.exit(
            f"--color-by column(s) not found: {missing}. "
            f"Available: {list(adata.obs.columns)}"
        )

    # Treat low-cardinality numeric/object columns (celltype codes, 0/1
    # condition flags, etc.) as categorical so scanpy uses discrete colors
    # instead of a continuous colormap.
    for c in color_by:
        if adata.obs[c].dtype.name != "category" and adata.obs[c].nunique() <= 20:
            adata.obs[c] = adata.obs[c].astype("category")

    bases = []  # list of (basis_name_without_X_prefix, panel_title)
    if args.compute:
        print("Computing PCA + neighbors + UMAP fresh (--compute)...")
        sc.pp.pca(adata, n_comps=args.n_pcs)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        bases = [("pca", "PCA (computed)"), ("umap", "UMAP (computed)")]
    else:
        if args.umap_obs_cols:
            xcol, ycol = args.umap_obs_cols
            missing_uc = [c for c in (xcol, ycol) if c not in adata.obs.columns]
            if missing_uc:
                sys.exit(
                    f"--umap-obs-cols column(s) not found: {missing_uc}. "
                    f"Available: {list(adata.obs.columns)}"
                )
            adata.obsm["X_umap_obscols"] = adata.obs[[xcol, ycol]].to_numpy(dtype=float)
            bases.append(("umap_obscols", f"UMAP ({xcol}/{ycol})"))
        elif "X_umap" in adata.obsm:
            bases.append(("umap", "UMAP (precomputed)"))
        if "X_pca" in adata.obsm:
            bases.append(("pca", "PCA (precomputed)"))
        if not bases:
            sys.exit(
                "No PCA/UMAP found in .obsm and no --umap-obs-cols given. "
                "Pass --compute to generate them fresh, or --umap-obs-cols "
                "if this file has 2D coordinates stored as plain obs "
                f"columns. Available obs columns: {list(adata.obs.columns)}"
            )

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_rows = len(bases)
    n_cols = len(color_by)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(4.5 * n_cols, 4.2 * n_rows), squeeze=False
    )
    for i, (basis, basis_title) in enumerate(bases):
        for j, color in enumerate(color_by):
            n_cats = adata.obs[color].nunique() if adata.obs[color].dtype.name == "category" else None
            drop_legend = (
                args.max_legend_categories > 0
                and n_cats is not None
                and n_cats > args.max_legend_categories
            )
            legend_loc = "none" if drop_legend else "right margin"
            if drop_legend:
                print(
                    f"Note: '{color}' has {n_cats} categories (> "
                    f"--max-legend-categories {args.max_legend_categories}) "
                    "-- dropping its legend to avoid overlap. Colors are "
                    "still shown; increase --max-legend-categories to force "
                    "it back, or pick a lower-cardinality --color-by column."
                )
            sc.pl.embedding(
                adata, basis=basis, color=color, ax=axes[i][j], show=False,
                size=args.point_size, title=f"{basis_title} -- {color}",
                legend_loc=legend_loc,
            )
    fig.suptitle(args.title or args.input)
    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight", dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
