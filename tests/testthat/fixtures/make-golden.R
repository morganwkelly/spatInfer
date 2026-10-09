# Record expected outputs from the current implementation.
# Run from the package root:  Rscript tests/testthat/fixtures/make-golden.R
# Optionally pass case names to regenerate only those.

pkgload::load_all(quiet = TRUE)
library(testthat)
source("tests/testthat/helper-golden.R")

cases <- commandArgs(trailingOnly = TRUE)
if (length(cases) == 0) cases <- names(golden_cases)
if (!requireNamespace("BRISC", quietly = TRUE)) cases <- setdiff(cases, "placebo_brisc")

out_dir <- file.path("tests", "testthat", "fixtures", "golden")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (name in cases) {
  message("Recording ", name)
  saveRDS(golden_cases[[name]](), file.path(out_dir, paste0(name, ".rds")), version = 3)
}
