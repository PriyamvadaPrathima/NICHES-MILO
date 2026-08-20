# NICHES+MILO

A flag-driven, dataset-agnostic cell-cell-communication pipeline: subset
your data to two conditions and the cell types you care about, build a
NICHES sender-receiver signaling network, then test which parts of that
network differ between conditions with Milo differential abundance. Comes
with a set of generic L-R heatmap/volcano figure scripts that read the
final output directly -- no dataset-specific code anywhere.

This started as a generalized reimplementation of the analysis pipeline
behind Bridges et al. (mapping intratumoral myeloid-T cell communication
with NICHES+Milo) -- see **Attribution** below -- rewritten to run on
**any** h5ad via command-line flags instead of hardcoded paths/columns/cell
types, and packaged here as a standalone repo.

For a step-by-step reference of exactly what each script reads and writes,
see **[PIPELINE_STEPS.md](PIPELINE_STEPS.md)**. For exact/recommended
dependency versions, see **[PACKAGE_VERSIONS.md](PACKAGE_VERSIONS.md)**.

## Scope: strict two-group comparisons

This pipeline tests one condition against one reference (e.g. treated vs.
control, mutant vs. wild-type, post- vs. pre-treatment) -- always exactly
two groups. It does not support >2-group designs or continuous covariates;
the Milo differential-abundance step uses the fixed design formula
`~ cond_binary`. If you need a more general design, `milo_helpers.py`'s
`run_milo()` is the pattern to extend.

## Files

| File | Language | Role |
|------|----------|------|
| `prep_niches_input.py` | Python | Step 1: subset to two conditions + cell types of interest, add bookkeeping columns |
| `run_niches_generic.R` | R | Step 2: per-sample ALRA imputation + NICHES cell-to-cell network generation |
| `run_milo_da.py` | Python | Step 3: Milo differential abundance, neighborhood clustering, optional Wilcoxon L-R export |
| `run_pipeline.ps1` / `run_pipeline.sh` | PowerShell / bash | Wrapper chaining all three steps -- edit the "EDIT THESE" block and run; `-Step`/first-arg lets you run just one step (see below) |
| `_lib.py` | Python | Shared helper functions (cluster selection, stats, patched milopy internals) used by Step 3 and all five plotting scripts |
| `milo_helpers.py` | Python | Bundled from Bridges et al. (see **Attribution**) -- Milo neighborhood clustering/plotting helpers |
| `plot_lr_heatmap.py` | Python | L-R enrichment heatmap per cluster -- which pairs define each cluster's identity |
| `plot_lr_comparison_heatmap.py` | Python | Reference- vs treatment-condition differential expression of L-R scores, per cluster (Cohen's d + significance) |
| `plot_lr_comparison_sidebyside.py` | Python | Two heatmaps side by side (reference vs treatment mean NICHES score), for direct visual comparison |
| `plot_lr_log2fc_heatmap.py` | Python | log2(treatment/reference) L-R fold-change heatmap per cluster, with significance overlay |
| `plot_milo_volcano.py` | Python | Volcano plot (logFC vs -log10 SpatialFDR) of Milo's own differential-abundance test, colored by cluster |
| `plot_embeddings.py` | Python | Quick-look PCA/UMAP plot for any h5ad, any pipeline stage |
| `install_niches.R` | R | One-off helper: remove + reinstall NICHES at the pinned tag (see **Known gotchas**) |
| `example_grouping_map.json` | -- | Template for the optional `--grouping-map` flag |
| `PIPELINE_STEPS.md` | -- | Step-by-step reference: exact inputs/outputs/columns for each step |
| `PACKAGE_VERSIONS.md` | -- | Every Python/R package this pipeline needs, with exact or recommended versions |
| `LICENSE-THIRD-PARTY.md` | -- | License/attribution for the bundled `milo_helpers.py` |

Every script also has `--help` (e.g. `uv run python prep_niches_input.py --help`)
listing every flag with its default.

## One-time setup

1. **Install `uv`, R, and Git**, and make sure all three are on PATH.
   - `uv`: `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`
     (Windows) or `curl -LsSf https://astral.sh/uv/install.sh | sh` (Linux/macOS).
   - R: install normally (CRAN installer, or RStudio); on Windows, add its
     `bin\x64` folder to PATH -- find the exact path from R itself with
     `normalizePath(R.home("bin"))` in the R/RStudio console.
   - Git: needed because several R packages below install from a GitHub
     commit/tag, not CRAN.
   - Restart your shell after each install for PATH changes to take effect.

2. **Sync the Python environment** from this repo's root:
   ```bash
   cd NICHES+MILO
   uv sync
   ```
   This reads `pyproject.toml`/`uv.lock` and builds a `.venv` with every
   Python dependency (scanpy, anndata, milopy, etc.) pinned to a known-good
   version. A venv built inside WSL is a *Linux* venv and can't be reused
   on native Windows (or vice versa) -- run `uv sync` from whichever
   environment you'll actually run the pipeline in.

3. **Pin `rpy2` to a known-good version.** A fresh `uv sync` resolves
   whatever's newest, which has a known S4-object round-tripping bug that
   breaks Step 3's call into R/edgeR (see **Known gotchas**):
   ```bash
   uv pip uninstall rpy2 rpy2-rinterface rpy2-robjects
   uv pip install "rpy2==3.5.10"
   ```
   This pin lives outside `pyproject.toml`/`uv.lock`, so a plain `uv run`
   would silently undo it by re-syncing first -- `run_pipeline.ps1`/
   `run_pipeline.sh` already set `UV_NO_SYNC=1` to prevent that. Don't run
   Step 1/3 with plain `uv run` outside the wrapper without that same
   environment variable set, or the pin will get reset.

4. **Install the R packages:**
   ```r
   install.packages("optparse")
   devtools::install_github("satijalab/seurat-wrappers")
   devtools::install_github("msraredon/NICHES", ref = "v1.2.5")
   devtools::install_github("saezlab/OmnipathR")
   devtools::install_github("cellgeni/sceasy")
   BiocManager::install(c("edgeR", "limma"))
   ```
   `sceasy` (used for the h5ad <-> Seurat conversion) drives Python's
   `anndata` package through `reticulate`. If the conversion step fails
   with a Python import error, point reticulate at this repo's uv
   environment: `reticulate::use_python("<this repo>/.venv/bin/python")`.

Once done, confirm Python can reach R:
```bash
uv run python -c "import rpy2.robjects as ro; print(ro.r('R.version.string'))"
```

## Adding your data

Put your QC'd, normalized, cell-type-annotated `.h5ad` file in a `data/`
subfolder here (create it if it doesn't exist):

```
NICHES+MILO/
  data/
    my_dataset.h5ad
  run_pipeline.ps1
  ...
```

Required columns in `.obs`, at minimum:
- a **sample/batch ID** column (one value per biological replicate)
- a **condition/treatment** column with at least the two values you want to compare
- a **cell type** column with the cell populations you want in the signaling network

This pipeline does not do QC, normalization, or cell-type annotation --
it picks up after that point.

## Configuring and running the pipeline

Edit the "EDIT THESE" block near the top of `run_pipeline.ps1` (Windows) or
`run_pipeline.sh` (Linux/WSL/macOS) -- see comments inline for what each
variable means, and "What to change to run a new dataset" below for the
full explanation of each one. Then run it from this repo's root:

```powershell
.\run_pipeline.ps1
```
```bash
./run_pipeline.sh
```

If PowerShell blocks the script with "running scripts is disabled on this
system," either run it once with
`powershell -ExecutionPolicy Bypass -File .\run_pipeline.ps1`, or fix it
permanently for your account with
`Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned`.

**Running just one step.** If you already have a `*_NICHES.h5ad` and just
want to rerun Milo with different settings (a different `--alpha`, say),
you don't need to redo Steps 1-2:
```powershell
.\run_pipeline.ps1 -Step Milo
```
```bash
./run_pipeline.sh milo
```
Valid values: `Full` (default), `Preprocess` (Step 1 only), `Niches` (Step
2 only), `Milo` (Step 3 only). Omit the argument entirely and the script
asks interactively instead. Each mode checks its own prerequisite file
exists first and fails with a clear message if not (e.g. `-Step Milo`
without an existing `*_NICHES.h5ad`).

Outputs land in `$OutDir`/`$OUT_DIR` (default `results/`), prefixed by
`$RunName`/`$RUN_NAME`: `*_LIM.h5ad` (Step 1), `*_NICHES.h5ad` (Step 2),
`*_MILO.h5ad` (Step 3, the final file the plotting scripts read),
`*_wilcoxon.xlsx`, and a timestamped `*_pipeline_log_*.txt` capturing that
run's full console output (warnings and info that don't make it into the
h5ad itself).

## Making figures

All five plotting scripts read the `*_MILO.h5ad` file Step 3 produces.
Run any of them with `--help` for the full flag list; the common ones:

```bash
# Which L-R pairs define each cluster's signaling identity
uv run python plot_lr_heatmap.py --input results/my_run_MILO.h5ad --output results/lr_heatmap.png

# Did a pair's score actually shift between conditions, per cluster (effect size + significance)
uv run python plot_lr_comparison_heatmap.py --input results/my_run_MILO.h5ad --output results/lr_comparison.png

# Raw reference/treatment mean scores side by side (not just the difference)
uv run python plot_lr_comparison_sidebyside.py --input results/my_run_MILO.h5ad --output results/lr_sidebyside.png

# log2(treatment/reference) fold-change heatmap
uv run python plot_lr_log2fc_heatmap.py --input results/my_run_MILO.h5ad --output results/lr_log2fc.png

# Volcano plot of Milo's own logFC vs -log10(SpatialFDR), colored by cluster
uv run python plot_milo_volcano.py --input results/my_run_MILO.h5ad --output results/milo_volcano.png
```

Useful flags across all five: `--clusters` (pick specific clusters by
number), `--top-n-clusters` (narrow to the N largest/most significant),
`--celltype-filter Macrophage` (only clusters dominated by a given cell
type), `--min-cluster-size` (drop small/unreliable clusters), and
`--lr-pairs` (an explicit L-R pair list instead of auto-selecting the top
cluster-defining pairs). Each script's own docstring (top of the file, or
`--help`) explains what it plots and why in more depth.

## What to change to run a new dataset

Everything dataset-specific lives in the `EDIT THESE` block near the top
of `run_pipeline.ps1`/`run_pipeline.sh` -- nothing else in this repo needs
to change. In order:

1. **Input path** -- path to your h5ad (relative to this repo's root, e.g.
   `data/my_dataset.h5ad`). Must already be QC'd, normalized, and
   cell-type-annotated.

2. **Sample column** -- the `.obs` column holding each cell's biological
   sample/batch/replicate ID. This is the unit NICHES imputes separately
   and Milo counts cells per, when testing for differential abundance --
   NOT the same thing as the condition/treatment column.

3. **Condition column, reference value, treatment value** -- the `.obs`
   column with your treatment/condition label, and which two of its
   values to compare (reference = encoded 0, treatment = encoded 1). This
   pipeline is strict two-group only. If your data has more than two
   conditions, pick exactly two per run -- rerun with different values for
   other pairwise comparisons.

4. **Cell-type column, cell types** -- the `.obs` column with your
   fine-grained per-cell-type labels, and which of its values to include
   as sender/receiver populations in the signaling network. Check
   `.obs[celltype_col].unique()` first and deliberately exclude any junk
   categories.

   Optionally, the grouping-map option (uncomment + point at a JSON file,
   see `example_grouping_map.json`) lets you collapse several fine cell
   types into one coarser NICHES sender/receiver category.

5. **Species** -- `'mouse'` or `'human'`, passed straight to NICHES.

6. **L-R database** -- ligand-receptor database passed to NICHES.
   `'omnipath'` (NICHES' own default) is used here. Fall back to
   `'fantom5'` if it acts up (connectivity issues, missing pairs).

7. **Output dir, run name** -- where outputs go and what they're
   prefixed with. No effect on the analysis itself.

8. **Annotation column** -- almost always needs to stay `'VectorType'`,
   not `run_milo_da.py`'s own default of `'celltype.Joint'`. NICHES only
   creates `celltype.Joint` when the mapped cell-type field passed to it
   is literally named `celltype`; this pipeline's internal `_grouping`
   column isn't, so `VectorType` (the sender-to-receiver pair NICHES
   always creates) is what you want. Leave this as `VectorType` unless
   you've specifically arranged for `celltype.Joint` to be created
   instead. See `PIPELINE_STEPS.md`'s "Column reference" for the full
   explanation.

9. **Replicate column** (optional) -- combines with the sample column to
   build a finer split unit (`<sample>-<replicate>`) when a single sample
   value actually pools multiple biological replicates (e.g.
   hashtag/CITE-seq-multiplexed samples). Leave unset if the sample column
   alone already identifies one biological replicate per value -- the
   common case. If Step 3 crashes with edgeR's `NA dispersions not
   allowed`, that means there weren't enough split units per condition for
   it to estimate a dispersion at all -- see `PIPELINE_STEPS.md`'s "Why the
   split unit encodes the condition" for the full explanation and whether
   this is your fix.

## Resuming after a crash

`run_niches_generic.R` (Step 2) is the slow step and checkpoints each
split's result to disk as it finishes, under a
`<output>_checkpoints_<species>_<lr-database>/` folder -- rerunning the
same command picks up where it left off rather than recomputing
everything. It also saves a full pre-write checkpoint
(`<output>_merged.rds`) right before the Seurat-to-h5ad conversion step.
`run_milo_da.py` (Step 3) checkpoints right after its call into R/edgeR
(`<output>_da.pkl`), since that's the step most likely to need retrying.
If a checkpoint file itself gets corrupted (e.g. the process was killed
mid-write), delete it and rerun -- the earlier, still-valid checkpoints
(per-split results, the DA pickle) are reused automatically, so you're
only redoing the step that actually failed.

## Known gotchas

- **Exactly two condition values.** If the condition column has more than
  two values, cells outside the reference/treatment values you pick are
  simply dropped by Step 1 -- not an error, just not part of the run.

- **`rpy2` must be pinned to `3.5.10`, not whatever `pyproject.toml`/
  `uv.lock` resolve by default.** Newer `rpy2` (confirmed problematic:
  3.6.x; confirmed problematic by others in the wild: 3.5.17) has had
  regressions in how R S4 objects round-trip through the Python/R
  boundary -- see
  [scverse/pertpy#681](https://github.com/scverse/pertpy/issues/681). See
  **One-time setup** above for the install command, and remember
  `UV_NO_SYNC=1` (or `uv run --no-sync`) for any ad-hoc command outside
  the pipeline wrapper, or a plain `uv run` will silently undo the pin.

- **`NICHES` should be pinned to `v1.2.5`, not `master`.** Earlier
  versions have a bug in the OmniPath connection path, fixed in this
  tagged release. `install_niches.R` reinstalls at this pin if you need to
  redo it.

- **Milo needs enough split units per condition to fit a dispersion.** If
  Step 3 crashes with edgeR's `NA dispersions not allowed`, you likely
  have too few biological replicates (split units) in one condition --
  see `PIPELINE_STEPS.md`'s "Why the split unit encodes the condition".

## Attribution

The overall NICHES+Milo statistical approach, and `milo_helpers.py`
(bundled here essentially unmodified -- see
[LICENSE-THIRD-PARTY.md](LICENSE-THIRD-PARTY.md)), come from:

> Bridges et al., **Bridgesetal-CCC**
> https://github.com/miller-jensen-lab/Bridgesetal-CCC

If you use this pipeline in published work, please cite that repository
and its associated paper.

This repository generalizes that original, dataset-specific analysis code
into flag-driven scripts that run on any h5ad, and adds the Milo-alpha
`-Step` runner, the comparison/log2FC/volcano figure scripts, and several
compatibility fixes for current `milopy`/`rpy2`/`networkx`/`pandas`
versions (see `_lib.py`'s module docstring for the full list).
