# Pipeline Steps -- NICHES+MILO

Step-by-step reference for the scripts in this repo. See `README.md` for
setup and a runnable quickstart.

## At a glance

```
your_data.h5ad (QC'd + annotated)
        |
        |  Step 1: prep_niches_input.py            (Python)
        v
*_LIM.h5ad
        |
        |  Step 2: run_niches_generic.R             (R)
        v
*_NICHES.h5ad     (one row = one sender-receiver cell pair)
        |
        |  Step 3: run_milo_da.py                   (Python)
        v
*_MILO.h5ad  (+ optional Wilcoxon .xlsx, cluster plot)
        |
        v
   your figures / enrichment analysis
```

| # | Step | Script | Language |
|---|------|--------|----------|
| 1 | Subset + prep | `prep_niches_input.py` | Python |
| 2 | NICHES network generation | `run_niches_generic.R` | R |
| 3 | Milo differential abundance | `run_milo_da.py` | Python |

Steps 0 and 4 aren't scripts in this repo, but bookend the pipeline -- see
below.

---

## Step 0: Before you start (not part of this repo)

Your own QC / normalization / cell-type-annotation pipeline produces the
input file. This pipeline picks up after that point.

**Needs to exist:** a QC'd, normalized, cell-type-annotated `.h5ad`, with
these columns in `.obs`:
- a sample/batch ID column
- a condition/treatment column
- a cell-type column

---

## Step 1: Subset + prep

**Script:** `prep_niches_input.py` &nbsp;&nbsp;**Language:** Python

**What it does:** Restricts your annotated h5ad to exactly the two
conditions you're comparing and the cell types you want in the signaling
network, then adds three bookkeeping columns the rest of the pipeline
depends on.

**Input**
- Your annotated `.h5ad` (Step 0)
- Required flags: `--sample-col`, `--condition-col`, `--reference-value`,
  `--treatment-value`, `--celltype-col`, `--cell-types`
- Optional flags: `--grouping-col`, `--grouping-map`, `--replicate-col`

**Output**
- `*_LIM.h5ad`, with three new `.obs` columns added:
  - `_condition_binary` -- 0/1 encoding of the comparison
  - `_grouping` -- coarse cell-type bucket passed to NICHES
  - `_niches_split` -- condition-tagged sample unit, e.g. `"1__BD3"`

---

## Step 2: NICHES network generation

**Script:** `run_niches_generic.R` &nbsp;&nbsp;**Language:** R

**What it does:** Converts to Seurat, splits by `_niches_split`, imputes
each split with ALRA, runs `RunNICHES` (cell-to-cell only) on each split,
tags every resulting cell-pair with its split-unit string, merges all
splits back together, filters out low-connectivity pairs, and computes
PCA/UMAP for visualization.

**Input**
- `*_LIM.h5ad` (Step 1)
- Key flags (all optional, sensible defaults): `--split-col` (default
  `_niches_split`), `--celltype-col` (default `_grouping`), `--species`,
  `--lr-database`, `--min-features`

**Output**
- `*_NICHES.h5ad` -- one row per sender-receiver cell pair, with a
  `Condition` column carrying the split-unit string from Step 1

---

## Step 3: Milo differential abundance

**Script:** `run_milo_da.py` &nbsp;&nbsp;**Language:** Python

**What it does:** Recovers the 0/1 condition from `Condition`, builds a kNN
graph over cell pairs, tests neighborhoods for differential abundance
between conditions (design `~ cond_binary`), Louvain-clusters the
significant neighborhoods, and assigns every cell pair to a cluster
(`sc_louvain`). Optionally runs Wilcoxon ligand-receptor enrichment per
cluster.

**Input**
- `*_NICHES.h5ad` (Step 2)
- Key flags (all optional, sensible defaults): `--split-col` (default
  `Condition`), `--condition-delim` (default `__`), `--annotation-col`
  (default `celltype.Joint`), `--alpha` (SpatialFDR threshold for the
  cluster plot -- also determines how many neighborhoods count as
  "significant" enough to cluster at all; too strict a default can leave
  every cell pair unassigned, see below), `--min-connect`, `--max-difflfc`

**Output**
- `*_MILO.h5ad` -- final output, with `sc_louvain` cluster assignments
- Optional: `--wilcox-out` writes a per-cluster Wilcoxon `.xlsx`;
  `--plot-out` saves the cluster plot as PDF/PNG

**If every cell pair ends up in `sc_louvain = -1`** (unassigned): Milo
found zero neighborhoods significant at your `--alpha`. This is a real,
fairly common failure mode with the default `--alpha 0.1` on noisier or
smaller datasets -- rerun Step 3 with a looser `--alpha` (e.g. `0.2`).
Since Step 3 checkpoints right after the expensive R/edgeR call
(`<output>_da.pkl`, see **Resuming after a crash** in `README.md`), a
rerun with a different `--alpha` reuses that checkpoint and only redoes
clustering -- fast, no need to touch Steps 1-2 at all. `-Step Milo` /
`milo` in the pipeline wrapper does exactly this.

---

## Step 4: After the pipeline (not part of this repo's Steps 1-3, but included here)

`*_MILO.h5ad` is what you build figures and enrichment analyses from. This
repo includes ready-to-run plotting scripts built directly against that
file: `plot_lr_heatmap.py`, `plot_lr_comparison_heatmap.py`,
`plot_lr_comparison_sidebyside.py`, `plot_lr_log2fc_heatmap.py`,
`plot_milo_volcano.py`, and `plot_embeddings.py` (for any stage's h5ad, not
just the final one). See `README.md`'s "Making figures".

---

## Why the split unit encodes the condition

This is the one non-obvious design choice worth understanding before you
run the pipeline -- it's how Steps 1-3 pass the condition assignment across
the R/Python boundary without a separate lookup file.

1. **Step 1** builds `_niches_split = "{0_or_1}__{sample}"`, e.g.
   `"0__BD2"` for a reference-condition sample, `"1__BD3"` for a
   treatment-condition sample. (If you pass `--replicate-col`, the sample
   portion becomes `"{sample}-{replicate}"`.)
2. **Step 2** splits on `_niches_split`, imputes + runs NICHES separately on
   each split (so imputation doesn't bleed information across samples),
   then tags every cell-pair in each split's output with
   `Condition = names(imp.list)[i]` -- literally the `_niches_split` string
   that produced it. After merging, `Condition` is a per-cell-pair column
   whose value is still `"{0_or_1}__{sample}"`.
3. **Step 3** splits `Condition` on the delimiter (`__` by default) and
   casts the first token back to an integer, recovering exactly the
   `_condition_binary` value assigned in Step 1. No dictionary or join is
   needed, and the R and Python sides can't disagree about which sample
   belongs to which condition, because the mapping never left the data.

`Condition`'s *distinct values* do a second job in Step 3, beyond the 0/1
prefix: Milo counts how many cells from each `_niches_split` unit fall in
each neighborhood (`rep_code = Condition.cat.codes`), which is what gives
the differential-abundance test its replication structure. That's why
Step 2 needs at least two split units total, and why more samples or
replicates per condition -- rather than more cells in one sample -- is what
actually improves the statistical power of Step 3. If you get a
`NA dispersions not allowed` error from edgeR in Step 3, this is almost
always the cause: not enough split units per condition. Adding a
`--replicate-col` in Step 1 (if your data has one -- e.g. hashtag/CITE-seq
multiplexing IDs) is the fix, since it gives each condition several split
units instead of one per sample.

---

## Column reference

Columns `prep_niches_input.py` adds in Step 1, and where each one gets used
downstream:

- **`_condition_binary`** -- 0 = `--reference-value`, 1 =
  `--treatment-value`. Optionally passed through Step 2 via `--meta-cols`.
- **`_grouping`** -- coarse cell-type bucket. Used in Step 2 as NICHES's
  `cell_types=` argument (sender/receiver identity).
- **`_niches_split`** -- per-sample[-replicate] unit for ALRA imputation.
  Encodes `_condition_binary` as a prefix. Used in Step 2 as `split.by=`.

Columns NICHES creates in Step 2, consumed in Step 3:

- **`Condition`** -- per-cell-pair copy of the `_niches_split` value that
  produced it. Source of both the DA design covariate and the replicate
  counting. Read via `--split-col`.
- **`celltype.Joint`** (or **`VectorType`**) -- created automatically by
  NICHES from the mapped celltype field. Labels each neighborhood by its
  dominant sender-receiver pair. Read via `--annotation-col`.

  Watch out if your `--celltype-col` starts with an underscore (true of the
  default, `_grouping`): somewhere in the Seurat/NICHES round-trip, R
  renames it by prepending `X` (`data.frame` auto-mangles names that aren't
  syntactically valid R identifiers on their own), so the real column ends
  up `X_grouping.Joint`, not `_grouping.Joint` or `celltype.Joint`. Always
  check the actual `obs` columns in the NICHES output
  (`adata.obs.columns.tolist()` in Python) before trusting a guessed
  `--annotation-col` value -- or just use `VectorType`, which doesn't have
  this issue since NICHES names it directly.

Columns `run_milo_da.py` adds in Step 3, used downstream:

- **`sc_louvain`** / **`louvain_str`** -- final communication-cluster
  assignment per cell pair (`-1` = not part of any significant
  neighborhood at the `--alpha` you ran with).
- **`cond_binary`** -- the plain 0/1 condition encoding, recovered from
  `Condition`. This is what the plotting scripts' `--condition-col`
  reads for the reference/treatment split -- distinct from `Condition`
  itself, which still carries the sample ID too.
