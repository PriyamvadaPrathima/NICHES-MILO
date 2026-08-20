#!/usr/bin/env Rscript
# run_niches_generic.R -- Step 2 of the generalized NICHES+Milo pipeline.
#
# Reads the h5ad written by prep_niches_input.py, converts it to Seurat,
# splits by the `_niches_split` column (condition-tagged sample[-replicate]
# unit), imputes each split with ALRA, runs NICHES cell-to-cell signaling on
# each, merges the results, and writes the merged network back out as h5ad
# for run_milo_da.py.
#
# R package setup (see README.md for the full list, including edgeR/limma
# needed by Step 3) -- install once:
#
#   Rscript -e 'install.packages("optparse")'
#   Rscript -e 'devtools::install_github("satijalab/seurat-wrappers")'
#   Rscript -e 'devtools::install_github("msraredon/NICHES", ref = "v1.2.5")'   # see
#     # README.md's "Known gotchas" for why this exact tag, not master/latest
#   Rscript -e 'devtools::install_github("saezlab/OmnipathR")'
#   Rscript -e 'devtools::install_github("cellgeni/sceasy")'
#
# sceasy additionally needs a Python env with `anndata` importable via
# reticulate -- if conversion fails with a Python import error, point
# reticulate at this project's uv env, e.g.:
#   Rscript -e 'reticulate::use_python("<this repo>/.venv/bin/python")'
#
# Usage:
#   Rscript run_niches_generic.R \
#     --input data/my_dataset_LIM.h5ad \
#     --output data/my_dataset_NICHES.h5ad \
#     --species mouse
#
# See PIPELINE_STEPS.md for the full step-by-step reference and README.md
# for setup instructions.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratWrappers)
  library(NICHES)
  library(sceasy)
})

# sceasy::convertFormat(..., from="seurat", to="anndata") dispatches to
# sceasy's own seurat2anndata() (see
# https://github.com/cellgeni/sceasy/blob/master/R/functions.R), which calls
# Seurat::GetAssayData(object = obj, assay = assay, slot = main_layer). The
# `slot` argument was renamed to `layer` in SeuratObject 5.0+ and is now
# fully defunct (errors, not just deprecated) -- sceasy (last released
# before that rename) hasn't been updated for it, so the final write below
# would otherwise crash every time with "The 'slot' argument ... is now
# defunct." This is a copy of seurat2anndata() with that one call fixed to
# use `layer=` instead, plus two further changes:
#   1. `assay` defaults to `Seurat::DefaultAssay(obj)` instead of sceasy's
#      hardcoded "RNA" -- our merged object's real (and only) assay is
#      "CellToCell" (NICHES' own naming), not "RNA", so the hardcoded
#      default would have failed differently (assay not found) right after
#      fixing the `slot`/`layer` issue.
#   2. `Seurat::GetAssay(obj, assay = assay)@meta.features` is wrapped in a
#      defensive tryCatch. Seurat5's default Assay5 class stores per-feature
#      metadata differently than the legacy Assay class sceasy was written
#      against, and this hasn't been tested against every Assay5 case -- a
#      failure here degrades to empty feature metadata instead of crashing
#      the whole write (`var` is otherwise cosmetic: h5ad requires some
#      var/gene-level table, but nothing downstream in this pipeline reads
#      per-L-R-pair metadata from it).
.write_seurat_as_h5ad <- function(obj, outFile, assay = Seurat::DefaultAssay(obj), main_layer = "data") {
  main_layer <- match.arg(main_layer, c("data", "counts", "scale.data"))

  if (compareVersion(as.character(obj@version), "3.0.0") < 0) {
    obj <- Seurat::UpdateSeuratObject(object = obj)
  }

  # THE FIX: `layer=`, not the removed-in-SeuratObject-5.0+ `slot=`.
  X <- Seurat::GetAssayData(object = obj, assay = assay, layer = main_layer)

  obs <- sceasy:::.regularise_df(obj@meta.data, drop_single_values = TRUE)

  var <- tryCatch(
    sceasy:::.regularise_df(Seurat::GetAssay(obj, assay = assay)@meta.features, drop_single_values = TRUE),
    error = function(e) {
      cat("  Note: couldn't read per-feature metadata from the Seurat5 assay object (",
          conditionMessage(e), ") -- writing empty var metadata instead.\n", sep = "")
      data.frame(row.names = rownames(obj[[assay]]))
    }
  )

  obsm <- NULL
  reductions <- names(obj@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(obj, reduction = name)),
      simplify = FALSE
    )
    names(obsm) <- paste0("X_", tolower(names(obj@reductions)))
  }

  anndata <- reticulate::import("anndata", convert = FALSE)

  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm
  )

  if (!is.null(outFile)) {
    adata$write(outFile, compression = "gzip")
  }

  adata
}

option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Input h5ad from prep_niches_input.py [required]"),
  make_option("--output", type = "character", default = NULL,
              help = "Output h5ad path for the merged NICHES cell-to-cell network [required]"),
  make_option("--split-col", type = "character", default = "_niches_split",
              help = "obs column to split on before per-group ALRA imputation [default: %default]"),
  make_option("--celltype-col", type = "character", default = "_grouping",
              help = "obs column with the coarse cell-type grouping passed to NICHES as cell_types= [default: %default]"),
  make_option("--species", type = "character", default = "mouse",
              help = "'mouse' or 'human', passed to RunNICHES [default: %default]"),
  make_option("--lr-database", type = "character", default = "omnipath",
              help = "Ligand-receptor database passed to RunNICHES [default: %default]"),
  make_option("--min-features", type = "integer", default = 5,
              help = "Minimum nFeature_CellToCell to retain a cell-pair connection [default: %default]"),
  make_option("--meta-cols", type = "character", default = NULL,
              help = "Comma-separated obs columns to map onto each cell-pair (each becomes <col>.Sending / <col>.Receiving, and <col>.Joint where they agree). Default: --celltype-col plus '_condition_binary'."),
  make_option("--checkpoint-dir", type = "character", default = NULL,
              help = "Directory to cache each split unit's ALRA+NICHES result to disk as soon as it finishes. If the script crashes later (a later split, merging, PCA/UMAP, or the final write), rerunning the exact same command loads any already-cached split from here instead of recomputing it. Default: '<output-without-extension>_checkpoints_<species>_<lr-database>/' -- species/lr-database are baked into the default path so switching either doesn't silently reuse a stale cache; --celltype-col/--meta-cols/--split-col are NOT keyed, so delete the directory yourself if you change those. Pass '' to disable."),
  make_option("--seurat-checkpoint", type = "character", default = NULL,
              help = "Path to cache the fully merged, filtered, post-PCA/UMAP Seurat object right before the Seurat->h5ad conversion (sceasy) -- the step most likely to crash from a Python/reticulate/SeuratObject-version mismatch. If this file already exists, the script skips splitting/ALRA/NICHES/merging/filtering/PCA/UMAP entirely, loads the object from here, and retries only the write. Default: '<output-without-extension>_merged.rds'. Pass '' to disable.")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "Step 2 of the generalized NICHES+Milo pipeline: per-sample ALRA imputation + NICHES cell-to-cell network generation."
))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("--input and --output are required. Run with --help for usage.", call. = FALSE)
}

split_col <- opt$`split-col`
celltype_col <- opt$`celltype-col`
out_base <- sub("\\.[^.]*$", "", opt$output)

checkpoint_dir <- opt$`checkpoint-dir`
if (is.null(checkpoint_dir)) {
  checkpoint_dir <- sprintf("%s_checkpoints_%s_%s", out_base, opt$species, opt$`lr-database`)
}
if (identical(checkpoint_dir, "")) checkpoint_dir <- NULL

seurat_checkpoint <- opt$`seurat-checkpoint`
if (is.null(seurat_checkpoint)) {
  seurat_checkpoint <- paste0(out_base, "_merged.rds")
}
if (identical(seurat_checkpoint, "")) seurat_checkpoint <- NULL

if (!is.null(seurat_checkpoint) && file.exists(seurat_checkpoint)) {

  cat("Found post-merge checkpoint, loading instead of recomputing:", seurat_checkpoint, "\n")
  merged <- readRDS(seurat_checkpoint)

} else {

  cat("Converting", opt$input, "to Seurat...\n")
  seurat_obj <- sceasy::convertFormat(opt$input, from = "anndata", to = "seurat", outFile = NULL)

  if (!(split_col %in% colnames(seurat_obj@meta.data))) {
    stop(sprintf("--split-col '%s' not found in the object's metadata. Available columns: %s",
                 split_col, paste(colnames(seurat_obj@meta.data), collapse = ", ")), call. = FALSE)
  }
  if (!(celltype_col %in% colnames(seurat_obj@meta.data))) {
    stop(sprintf("--celltype-col '%s' not found in the object's metadata. Available columns: %s",
                 celltype_col, paste(colnames(seurat_obj@meta.data), collapse = ", ")), call. = FALSE)
  }

  meta_cols <- if (!is.null(opt$`meta-cols`)) {
    trimws(strsplit(opt$`meta-cols`, ",")[[1]])
  } else {
    unique(c(celltype_col, "_condition_binary"))
  }
  missing_meta <- setdiff(meta_cols, colnames(seurat_obj@meta.data))
  if (length(missing_meta) > 0) {
    stop(sprintf("--meta-cols reference column(s) not present in the object: %s",
                 paste(missing_meta, collapse = ", ")), call. = FALSE)
  }

  cat("Splitting by", split_col, "...\n")
  data.list <- SplitObject(seurat_obj, split.by = split_col)
  cat(sprintf("  %d split unit(s): %s\n", length(data.list), paste(names(data.list), collapse = ", ")))
  if (length(data.list) < 2) {
    stop("Only one split unit found -- NICHES needs at least two samples ",
         "(one per condition, minimum) to generate a comparable network.",
         call. = FALSE)
  }

  if (!is.null(checkpoint_dir) && !dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE)
  }
  if (!is.null(checkpoint_dir)) {
    cat("Per-split checkpoints:", checkpoint_dir, "\n")
  }

  imp.list <- list()
  for (i in seq_along(data.list)) {
    split_name <- names(data.list)[i]
    ckpt_file <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, paste0(split_name, ".rds")) else NULL

    if (!is.null(ckpt_file) && file.exists(ckpt_file)) {
      cat(sprintf("  [%d/%d] split '%s': loading cached result from %s\n",
                  i, length(data.list), split_name, ckpt_file))
      imp.list[[i]] <- readRDS(ckpt_file)
    } else {
      cat(sprintf("  [%d/%d] ALRA + NICHES on split '%s' (%d cells)\n",
                  i, length(data.list), split_name, ncol(data.list[[i]])))
      imputed <- SeuratWrappers::RunALRA(data.list[[i]])
      imp.list[[i]] <- RunNICHES(imputed,
                                  LR.database = opt$`lr-database`,
                                  species = opt$species,
                                  assay = "alra",
                                  cell_types = celltype_col,
                                  meta.data.to.map = meta_cols,
                                  SystemToCell = TRUE,
                                  CellToCell = TRUE)
      if (!is.null(ckpt_file)) {
        saveRDS(imp.list[[i]], ckpt_file)
        cat("    -> cached to", ckpt_file, "\n")
      }
    }
  }
  names(imp.list) <- names(data.list)

  cat("Merging CellToCell networks...\n")
  temp.list <- list()
  for (i in seq_along(imp.list)) {
    temp.list[[i]] <- imp.list[[i]]$CellToCell
    # 'Condition' carries the exact split-unit value (e.g. "1__sampleA"), which
    # run_milo_da.py parses to recover the 0/1 condition encoding -- see
    # PIPELINE_STEPS.md, "Why the split unit encodes the condition".
    temp.list[[i]]$Condition <- names(imp.list)[i]
  }

  # NOTE: the original run-niches.R merges list elements with
  # `for (k in length(temp.list)-1) merge(...)`, which is a for-loop over a
  # single scalar and silently merges only the first and last elements. That
  # was invisible in the original analyses because they always had exactly two
  # split units per run. This generalized version supports an arbitrary number
  # of samples/replicates per condition, so it merges all of them properly:
  merged <- temp.list[[1]]
  if (length(temp.list) > 1) {
    merged <- merge(merged, y = temp.list[2:length(temp.list)])
  }

  cat(sprintf("Filtering to nFeature_CellToCell > %d...\n", opt$`min-features`))
  merged <- subset(merged, nFeature_CellToCell > opt$`min-features`)
  cat(sprintf("  %d cell-pairs retained\n", ncol(merged)))

  cat("Computing PCA/UMAP for visualization...\n")
  merged <- ScaleData(merged)
  merged <- FindVariableFeatures(merged, selection.method = "disp")
  merged <- RunPCA(merged, npcs = 50)
  merged <- RunUMAP(merged, dims = 1:25)

  if (!is.null(seurat_checkpoint)) {
    saveRDS(merged, seurat_checkpoint)
    cat("Saved post-merge checkpoint:", seurat_checkpoint, "\n")
  }
}

# `merged` was built by merging separate per-split Seurat objects.
# SeuratObject 5's merge() keeps each source's data as its OWN layer
# (data.1, data.2, ...) rather than combining them -- deliberately, so
# layer-aware functions like ScaleData/FindVariableFeatures/RunPCA/RunUMAP
# above can still operate correctly without joining first (and did: this is
# why the checkpoint above already has working PCA/UMAP even though
# JoinLayers() hasn't run yet). GetAssayData() inside
# .write_seurat_as_h5ad(), by contrast, asks for one specific layer by
# name ("data") and fails with "GetAssayData doesn't work for multiple
# layers in v5 assay" once there's more than one. JoinLayers() consolidates
# them into a single "data" layer right before the one call that actually
# needs that. Safe to call even if --seurat-checkpoint loaded an object
# from an older run that predates this fix, and safe/harmless if the
# layers are already joined.
merged <- SeuratObject::JoinLayers(merged)

cat("Writing", opt$output, "...\n")
.write_seurat_as_h5ad(merged, opt$output)
cat("Done.\n")
