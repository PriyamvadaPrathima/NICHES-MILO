#!/usr/bin/env bash
# run_pipeline.sh -- chains the three NICHES+Milo pipeline steps (Linux/WSL/macOS).
#
# Edit the variables below to match your dataset, then run from inside this
# repo's root folder:
#   ./run_pipeline.sh
#
# To run only part of the pipeline (e.g. you already have a *_NICHES.h5ad
# and just want to rerun Milo with different settings), pass a step name:
#   ./run_pipeline.sh milo
# Valid values: full (default), preprocess, niches, milo. Omit it entirely
# and the script asks interactively instead.
#
# Each step can also be run and inspected individually -- see
# PIPELINE_STEPS.md for what each one reads/writes, and README.md for
# one-time setup (uv + the R packages run_niches_generic.R needs).
#
# The Python steps use `uv run --no-sync` rather than plain `uv run`: this
# pipeline needs rpy2 pinned to 3.5.10 (newer rpy2 has a numpy<->R
# conversion bug that breaks Milo's edgeR call -- see README's "Known
# gotchas"), which is installed as a manual override outside
# pyproject.toml/uv.lock. Plain `uv run` re-syncs the venv against the
# lockfile before every run and would silently undo that override;
# `--no-sync` skips that step.
#
# Step 2 (R) also needs Python, indirectly: run_niches_generic.R uses
# sceasy::convertFormat() to go from h5ad to Seurat, which calls Python's
# `anndata` package under the hood via R's reticulate package. Without
# guidance, reticulate can auto-provision its own separate, empty Python
# environment (via uv) that doesn't have anndata installed, instead of
# reusing this repo's own venv that Step 1/3 already use (and which already
# has anndata/scanpy/pandas from `uv sync`) -- causing a confusing
# `ModuleNotFoundError: No module named 'anndata'` deep inside an R error,
# with no obvious R-side fix. Below, RETICULATE_PYTHON is pointed explicitly
# at that same venv's python so reticulate reuses it instead of creating a
# new one.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ---- Which step(s) to run ---------------------------------------------------
STEP="${1:-}"
if [[ -z "$STEP" ]]; then
  echo ""
  echo "Which step(s) do you want to run?"
  echo "  1) Full pipeline (preprocess + NICHES + Milo)  [default]"
  echo "  2) NICHES only (Step 2)  -- needs an existing *_LIM.h5ad"
  echo "  3) Milo only (Step 3)    -- needs an existing *_NICHES.h5ad"
  echo "  4) Preprocessing only (Step 1)"
  read -r -p "Enter 1-4: " choice
  case "$choice" in
    2) STEP="niches" ;;
    3) STEP="milo" ;;
    4) STEP="preprocess" ;;
    *) STEP="full" ;;
  esac
fi
case "$STEP" in
  full|preprocess|niches|milo) ;;
  *) echo "ERROR: unrecognized step '$STEP' -- expected full, preprocess, niches, or milo." >&2; exit 1 ;;
esac
RUN_STEP1=false; RUN_STEP2=false; RUN_STEP3=false
[[ "$STEP" == "full" || "$STEP" == "preprocess" ]] && RUN_STEP1=true
[[ "$STEP" == "full" || "$STEP" == "niches" ]] && RUN_STEP2=true
[[ "$STEP" == "full" || "$STEP" == "milo" ]] && RUN_STEP3=true
echo "Running mode: $STEP"

# ============================== EDIT THESE ==================================
INPUT_H5AD="data/your_dataset.h5ad"   # QC'd, normalized, cell-type-annotated
SAMPLE_COL="sample"                   # obs column: biological sample/batch
CONDITION_COL="condition"             # obs column: the two groups you're comparing
REFERENCE_VALUE="control"             # condition value encoded as 0
TREATMENT_VALUE="treated"             # condition value encoded as 1
CELLTYPE_COL="cell_type"              # obs column: fine cell type labels
CELL_TYPES=("CellTypeA" "CellTypeB" "CellTypeC")   # which values of CELLTYPE_COL to include
# GROUPING_MAP="example_grouping_map.json"  # uncomment to collapse celltypes
# REPLICATE_COL="hashing"               # uncomment if SAMPLE_COL doesn't already
                                         # give one row per biological replicate
SPECIES="human"                       # 'mouse' or 'human', passed to NICHES
LR_DATABASE="omnipath"                # ligand-receptor DB passed to NICHES
OUT_DIR="results"
RUN_NAME="my_dataset"
ANNOTATION_COL="VectorType"           # Step 3: celltype.Joint (the default) isn't
                                       # always created -- see PIPELINE_STEPS.md

# Set true to skip Step 1 entirely when its output file already exists --
# a plain file-existence check, NOT re-validated against SAMPLE_COL/
# CELL_TYPES/etc. above. If you change anything above Step 1 uses, delete
# $LIM_H5AD (or leave this false) rather than risk a stale mismatch.
SKIP_STEP1_IF_EXISTS=false
# ==============================================================================

if ! command -v Rscript >/dev/null 2>&1; then
  if $RUN_STEP2; then
    echo "ERROR: 'Rscript' not found on PATH -- Step 2 (NICHES network generation) needs R installed and on PATH." >&2
    echo "  If R lives in a separate conda environment, activate it first, e.g.: conda activate r-env" >&2
    exit 1
  fi
fi

RETICULATE_VENV_PYTHON="$SCRIPT_DIR/.venv/bin/python"
if [[ -x "$RETICULATE_VENV_PYTHON" ]]; then
  export RETICULATE_PYTHON="$RETICULATE_VENV_PYTHON"
elif $RUN_STEP2; then
  echo "WARNING: no Python venv found at $RETICULATE_VENV_PYTHON -- R's reticulate (used by Step 2's sceasy conversion) may fall back to provisioning its own Python and fail with \"ModuleNotFoundError: No module named 'anndata'\". Run 'uv sync' first if Step 2 fails with that error." >&2
fi

LIM_H5AD="${OUT_DIR}/${RUN_NAME}_LIM.h5ad"
NICHES_H5AD="${OUT_DIR}/${RUN_NAME}_NICHES.h5ad"
MILO_H5AD="${OUT_DIR}/${RUN_NAME}_MILO.h5ad"
WILCOX_XLSX="${OUT_DIR}/${RUN_NAME}_wilcoxon.xlsx"

GROUPING_ARGS=()
if [[ -n "${GROUPING_MAP:-}" ]]; then
  GROUPING_ARGS=(--grouping-map "$GROUPING_MAP")
fi
REPLICATE_ARGS=()
if [[ -n "${REPLICATE_COL:-}" ]]; then
  REPLICATE_ARGS=(--replicate-col "$REPLICATE_COL")
fi

echo "=== Step 1/3: subset + prep (Python) ==="
if ! $RUN_STEP1; then
  echo "Skipping (mode: $STEP)."
elif [[ "$SKIP_STEP1_IF_EXISTS" == true && -f "$LIM_H5AD" ]]; then
  echo "Skipping -- $LIM_H5AD already exists and SKIP_STEP1_IF_EXISTS=true."
  echo "  (Not re-checked against SAMPLE_COL/CELL_TYPES/etc. above -- delete the file if those changed.)"
else
  uv run --no-sync python prep_niches_input.py \
      --input "$INPUT_H5AD" --output "$LIM_H5AD" \
      --sample-col "$SAMPLE_COL" \
      --condition-col "$CONDITION_COL" \
      --reference-value "$REFERENCE_VALUE" \
      --treatment-value "$TREATMENT_VALUE" \
      --celltype-col "$CELLTYPE_COL" \
      --cell-types "${CELL_TYPES[@]}" \
      "${GROUPING_ARGS[@]}" "${REPLICATE_ARGS[@]}"
fi

if $RUN_STEP2 && [[ ! -f "$LIM_H5AD" ]]; then
  echo "ERROR: Step 2 needs '$LIM_H5AD', which doesn't exist -- run with 'full' or 'preprocess' first (or point an existing prepped file at that path)." >&2
  exit 1
fi
echo "=== Step 2/3: NICHES network generation (R) ==="
if $RUN_STEP2; then
  Rscript run_niches_generic.R \
      --input "$LIM_H5AD" --output "$NICHES_H5AD" --species "$SPECIES" --lr-database "$LR_DATABASE"
else
  echo "Skipping (mode: $STEP)."
fi

if $RUN_STEP3 && [[ ! -f "$NICHES_H5AD" ]]; then
  echo "ERROR: Step 3 needs '$NICHES_H5AD', which doesn't exist -- run with 'full' or 'niches' first (or point an existing NICHES h5ad at that path)." >&2
  exit 1
fi
echo "=== Step 3/3: Milo differential abundance (Python) ==="
if $RUN_STEP3; then
  uv run --no-sync python run_milo_da.py \
      --input "$NICHES_H5AD" --output "$MILO_H5AD" \
      --annotation-col "$ANNOTATION_COL" --wilcox-out "$WILCOX_XLSX"
else
  echo "Skipping (mode: $STEP)."
fi

echo ""
echo "Done."
$RUN_STEP1 && echo "  Prepped input:          $LIM_H5AD"
$RUN_STEP2 && echo "  NICHES network:         $NICHES_H5AD"
if $RUN_STEP3; then
  echo "  Milo+cluster output:    $MILO_H5AD"
  echo "  Wilcoxon L-R tables:    $WILCOX_XLSX"
fi
