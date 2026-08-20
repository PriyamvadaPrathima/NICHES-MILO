# run_pipeline.ps1 -- chains the three NICHES+Milo pipeline steps entirely
# on native Windows (no WSL).
#
# This is a TEMPLATE -- the $InputH5ad/$SampleCol/etc. values below are
# placeholders. Edit the "EDIT THESE" block to match your dataset (see
# README.md's "Running on a new dataset" section for exactly what each one
# means), then run from inside this repo's root folder:
#   .\run_pipeline.ps1
#
# To run only part of the pipeline (e.g. you already have a *_NICHES.h5ad
# and just want to rerun Milo with different settings), pass -Step:
#   .\run_pipeline.ps1 -Step Milo
# Valid values: Full (default), Preprocess, Niches, Milo. Omit -Step
# entirely and the script asks interactively instead.
#
# One-time setup this script assumes is already done on whatever machine
# it's running on (see README.md):
#   - uv, R (Rscript.exe), and Git all installed and on PATH
#   - This repo's venv synced via `uv sync` from this folder (a
#     Windows-native .venv, NOT a WSL one -- Linux and Windows venvs are not
#     interchangeable, so a venv built on one machine can't just be copied
#     to another either; each machine needs its own `uv sync`)
#   - rpy2 pinned to 3.5.10 in that venv (`uv pip install "rpy2==3.5.10"`)
#   - R packages installed: optparse, Seurat (5.0+), SeuratWrappers, NICHES,
#     sceasy, OmnipathR, edgeR, limma

param(
    # Which step(s) to run. Omit this to get an interactive menu instead.
    #   Full       -- steps 1+2+3, the whole pipeline
    #   Preprocess -- step 1 only (subset + prep)
    #   Niches     -- step 2 only (NICHES network generation) -- needs an
    #                 existing *_LIM.h5ad from a prior Preprocess/Full run
    #   Milo       -- step 3 only (Milo differential abundance) -- needs an
    #                 existing *_NICHES.h5ad from a prior Niches/Full run
    # Example: .\run_pipeline.ps1 -Step Milo
    [ValidateSet("Full", "Preprocess", "Niches", "Milo")]
    [string]$Step
)

$ErrorActionPreference = "Stop"

if (-not $Step) {
    Write-Host ""
    Write-Host "Which step(s) do you want to run?"
    Write-Host "  1) Full pipeline (preprocess + NICHES + Milo)  [default]"
    Write-Host "  2) NICHES only (Step 2)  -- needs an existing *_LIM.h5ad"
    Write-Host "  3) Milo only (Step 3)    -- needs an existing *_NICHES.h5ad"
    Write-Host "  4) Preprocessing only (Step 1)"
    $choice = Read-Host "Enter 1-4"
    $Step = switch ($choice) {
        "2" { "Niches" }
        "3" { "Milo" }
        "4" { "Preprocess" }
        default { "Full" }
    }
}
$RunStep1 = $Step -in @("Full", "Preprocess")
$RunStep2 = $Step -in @("Full", "Niches")
$RunStep3 = $Step -in @("Full", "Milo")
Write-Host "Running mode: $Step"

# ============================== EDIT THESE ==================================
$InputH5ad      = "data\your_dataset.h5ad"    # QC'd, normalized, cell-type-annotated
$SampleCol      = "sample"                    # obs column: biological sample/batch
$ConditionCol   = "condition"                 # obs column: the two groups you're comparing
$ReferenceValue = "control"                   # condition value encoded as 0
$TreatmentValue = "treated"                   # condition value encoded as 1
$CelltypeCol    = "cell_type"                 # obs column: fine cell type labels
$CellTypes      = @(
    "CellTypeA", "CellTypeB", "CellTypeC"     # which values of $CelltypeCol to include
)
# $GroupingMap  = "example_grouping_map.json"  # uncomment to collapse celltypes into
                                               # coarser sender/receiver buckets
# $ReplicateCol = "hashing"                    # uncomment if $SampleCol doesn't already
                                               # give one row per biological replicate
$Species        = "human"                     # 'mouse' or 'human', passed to NICHES
$LrDatabase     = "omnipath"                  # ligand-receptor DB passed to NICHES --
                                               # its own default; fall back to "fantom5"
                                               # if the omnipath connection acts up
$OutDir         = "results"
$RunName        = "my_dataset"
$AnnotationCol  = "VectorType"                # Step 3: celltype.Joint (the default) isn't
                                               # always created -- see PIPELINE_STEPS.md's
                                               # "Column reference" for why VectorType is
                                               # the safe default

# Set $true to skip Step 1 entirely when its output file already exists --
# a plain file-existence check, NOT re-validated against SampleCol/CellTypes/
# etc. above. If you change anything above Step 1 uses, delete the *_LIM.h5ad
# file (or leave this $false) rather than risk a stale mismatch.
$SkipStep1IfExists = $false
# ==============================================================================

# rpy2 pin protection (see README's "Known gotchas") -- plain `uv sync`
# would otherwise resolve whatever's newest and silently undo the pin.
$env:UV_NO_SYNC = "1"

# ---- Logging ----------------------------------------------------------------
# Every run's full console output -- Step 1/2/3 stdout+stderr, R package
# startup messages, NICHES/Milo warnings, dropped-cell counts, checkpoint
# hits, everything -- is captured to a timestamped log file under $OutDir, in
# addition to being shown live in this window.
if (-not (Test-Path $OutDir)) {
    New-Item -ItemType Directory -Path $OutDir | Out-Null
}
$LogFile = Join-Path $OutDir "${RunName}_pipeline_log_$(Get-Date -Format 'yyyyMMdd_HHmmss').txt"
Start-Transcript -Path $LogFile | Out-Null
Write-Host "Logging full output to: $LogFile"

try {

# ---- Preflight checks ------------------------------------------------------
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Error "ERROR: 'uv' not found on PATH."
    exit 1
}
if ($RunStep2 -and -not (Get-Command Rscript.exe -ErrorAction SilentlyContinue)) {
    Write-Error "ERROR: 'Rscript.exe' not found on PATH -- Step 2 needs R installed and on PATH."
    exit 1
}
$PyProjectFile = Join-Path $PSScriptRoot "pyproject.toml"
if (-not (Test-Path $PyProjectFile)) {
    Write-Error "ERROR: no pyproject.toml found next to this script ('$PyProjectFile') -- run this from inside the repo, not a copy of just this file."
    exit 1
}

$LimH5ad    = Join-Path $OutDir "${RunName}_LIM.h5ad"
$NichesH5ad = Join-Path $OutDir "${RunName}_NICHES.h5ad"
$MiloH5ad   = Join-Path $OutDir "${RunName}_MILO.h5ad"
$WilcoxXlsx = Join-Path $OutDir "${RunName}_wilcoxon.xlsx"

$GroupingArgs = @()
if ($GroupingMap) {
    $GroupingArgs = @("--grouping-map", $GroupingMap)
}
$ReplicateArgs = @()
if ($ReplicateCol) {
    $ReplicateArgs = @("--replicate-col", $ReplicateCol)
}

Write-Host "=== Step 1/3: subset + prep (Python) ==="
if (-not $RunStep1) {
    Write-Host "Skipping (mode: $Step)."
} elseif ($SkipStep1IfExists -and (Test-Path $LimH5ad)) {
    Write-Host "Skipping -- $LimH5ad already exists and SkipStep1IfExists=`$true."
    Write-Host "  (Not re-checked against SampleCol/CellTypes/etc. above -- delete the file if those changed.)"
} else {
    $step1Args = @(
        "--input", $InputH5ad, "--output", $LimH5ad,
        "--sample-col", $SampleCol,
        "--condition-col", $ConditionCol,
        "--reference-value", $ReferenceValue,
        "--treatment-value", $TreatmentValue,
        "--celltype-col", $CelltypeCol,
        "--cell-types"
    ) + $CellTypes + $GroupingArgs + $ReplicateArgs

    uv run --no-sync python prep_niches_input.py @step1Args
    if ($LASTEXITCODE -ne 0) { Write-Error "Step 1 failed (exit $LASTEXITCODE)"; exit 1 }
}

if ($RunStep2 -and -not (Test-Path $LimH5ad)) {
    Write-Error "Step 2 needs '$LimH5ad', which doesn't exist -- run with -Step Full or -Step Preprocess first (or point an existing prepped file at that path)."
    exit 1
}
Write-Host "=== Step 2/3: NICHES network generation (R) ==="
if ($RunStep2) {
    Rscript.exe run_niches_generic.R --input $LimH5ad --output $NichesH5ad --species $Species --lr-database $LrDatabase
    if ($LASTEXITCODE -ne 0) { Write-Error "Step 2 failed (exit $LASTEXITCODE)"; exit 1 }
} else {
    Write-Host "Skipping (mode: $Step)."
}

if ($RunStep3 -and -not (Test-Path $NichesH5ad)) {
    Write-Error "Step 3 needs '$NichesH5ad', which doesn't exist -- run with -Step Full or -Step Niches first (or point an existing NICHES h5ad at that path)."
    exit 1
}
Write-Host "=== Step 3/3: Milo differential abundance (Python) ==="
if ($RunStep3) {
    uv run --no-sync python run_milo_da.py --input $NichesH5ad --output $MiloH5ad --annotation-col $AnnotationCol --wilcox-out $WilcoxXlsx
    if ($LASTEXITCODE -ne 0) { Write-Error "Step 3 failed (exit $LASTEXITCODE)"; exit 1 }
} else {
    Write-Host "Skipping (mode: $Step)."
}

Write-Host ""
Write-Host "Done."
if ($RunStep1) { Write-Host "  Prepped input:          $LimH5ad" }
if ($RunStep2) { Write-Host "  NICHES network:         $NichesH5ad" }
if ($RunStep3) {
    Write-Host "  Milo+cluster output:    $MiloH5ad"
    Write-Host "  Wilcoxon L-R tables:    $WilcoxXlsx"
}
Write-Host "  Full run log:           $LogFile"

} finally {
    Stop-Transcript | Out-Null
}
