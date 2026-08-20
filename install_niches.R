#!/usr/bin/env Rscript
# reinstall_niches.R -- one-off helper to pick up NICHES' upstream omnipath
# fix. Remove and reinstall from GitHub master, since NICHES isn't on CRAN.
#
# Run with:
#   Rscript.exe reinstall_niches.R
#
# (Written as a file rather than `Rscript.exe -e '...'` because PowerShell's
# argument passing to native executables can strip the double quotes around
# string literals in an inline -e command -- R then sees a bare, undefined
# `NICHES` symbol instead of the string "NICHES".)

if ("NICHES" %in% rownames(installed.packages())) {
  cat("Removing existing NICHES install...\n")
  remove.packages("NICHES")
} else {
  cat("NICHES not currently installed -- skipping removal.\n")
}

cat("Installing NICHES from GitHub (msraredon/NICHES@v1.2.5)...\n")
devtools::install_github("msraredon/NICHES", ref = "v1.2.5", force = TRUE)

cat("Done. Installed version:\n")
print(packageVersion("NICHES"))
