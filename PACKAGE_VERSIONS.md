# Package Versions

Reference list of everything this pipeline (Steps 1-3 + all five plotting
scripts) needs, with exact versions where they're pinned/locked, and
recommended versions where R doesn't give us a lockfile to point to. See
`README.md`'s "One-time setup" for the install commands themselves --
this file is just the version reference to check installs against.

## System tools

| Tool | Version | Notes |
|---|---|---|
| `uv` | latest | Installer script always grabs current; no known version sensitivity |
| R | >= 4.3 | Required by Seurat 5.x |
| Rscript / Rscript.exe | (same install as R) | Must be on PATH |
| Git | any recent | Needed because several R packages below install from a GitHub commit/tag, not CRAN |

## Python (managed by `uv`, no manual version-chasing needed)

Everything below is resolved and pinned by this repo's own `uv.lock` --
running `uv sync` from the repo root installs exactly these versions
automatically. Listed here so you can spot a mismatch (`uv pip list`) if
something behaves differently than expected.

| Package | Role |
|---|---|
| numpy | pinned `<2` -- some dependencies below aren't numpy-2-tested |
| pandas | core dataframes |
| scipy | stats (Mann-Whitney U, etc.) |
| matplotlib | all figure scripts |
| scanpy | single-cell I/O and analysis |
| anndata | h5ad I/O |
| scikit-learn | `euclidean_distances`, used by `milo_helpers.py` |
| networkx | neighborhood-graph clustering (`_lib.cluster_nhoods_fixed`) |
| python-louvain | Louvain community detection (imported as `community`) |
| milopy | installed from `github.com/emdann/milopy` @ commit `30646f5`, not PyPI |
| distinctipy | cluster color palettes (`milo_helpers.py`) |
| xlsxwriter | `--wilcox-out` Excel export (Step 3) |

**One deliberate exception, installed *outside* the lockfile:**

| Package | Version | Why it's pinned manually |
|---|---|---|
| `rpy2` | **3.5.10** (not whatever `uv.lock` would otherwise resolve, e.g. 3.6.x) | Newer `rpy2` (3.6.x, and reportedly 3.5.17) has a regression in how R S4 objects round-trip through the Python/R boundary, breaking Step 3's call into edgeR. See `README.md`'s "Known gotchas" for the full story and the [upstream issue](https://github.com/scverse/pertpy/issues/681). |

Install/re-pin it with:
```bash
uv pip uninstall rpy2 rpy2-rinterface rpy2-robjects
uv pip install "rpy2==3.5.10"
```
Because this sits outside `uv.lock`, a plain `uv run` re-syncs the venv
against the lockfile first and silently reverts it -- always run Step 1/3
with `UV_NO_SYNC=1` set (already done at the top of `run_pipeline.ps1` /
`run_pipeline.sh`), or pass `uv run --no-sync` yourself for any ad-hoc
command.

## R packages

R packages here aren't locked to exact versions the way `uv.lock` locks
Python (no `renv.lock` in this repo) -- install the current release of
each unless noted otherwise, and use `packageVersion("<pkg>")` in R to
check what you actually have if something behaves unexpectedly.

| Package | Source | Version | Notes |
|---|---|---|---|
| `optparse` | CRAN | latest | CLI flag parsing for `run_niches_generic.R` |
| `Seurat` | CRAN or GitHub | >= 5.0 | Needed for the `layer=`/`JoinLayers()` API this pipeline relies on; older Seurat 4.x will not work |
| `SeuratWrappers` | GitHub (`satijalab/seurat-wrappers`) | latest | Provides `RunALRA()` |
| `NICHES` | GitHub (`msraredon/NICHES`) | **`v1.2.5`** (pin via `ref = "v1.2.5"`) | Earlier versions have a bug in the OmniPath connection path -- fixed in this tagged release. `install_niches.R` in this repo reinstalls at this pin. Note: `packageVersion("NICHES")` may still report `1.2.4` even after installing this tag correctly (the maintainer didn't bump `DESCRIPTION`'s `Version:` field for this release) -- confirm via the installed commit hash in the install log instead if the version number looks stale. |
| `OmnipathR` | GitHub (`saezlab/OmnipathR`) | latest | LR database backend (`--lr-database omnipath`) |
| `sceasy` | GitHub (`cellgeni/sceasy`) | latest | h5ad <-> Seurat conversion |
| `edgeR` | Bioconductor | latest | Milo's differential-abundance NB-GLM fit (Step 3, via `rpy2`) |
| `limma` | Bioconductor | latest | edgeR dependency for the fit/test chain |

Install commands (see `README.md` for the full one-time setup sequence):
```r
install.packages("optparse")
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github("msraredon/NICHES", ref = "v1.2.5")
devtools::install_github("saezlab/OmnipathR")
devtools::install_github("cellgeni/sceasy")
BiocManager::install(c("edgeR", "limma"))
```

## Checking what's actually installed

Python:
```bash
uv pip list
```

R:
```r
sapply(c("optparse","Seurat","SeuratWrappers","NICHES","OmnipathR","sceasy","edgeR","limma"),
       function(p) as.character(packageVersion(p)))
```
