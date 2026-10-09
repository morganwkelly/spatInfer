# Shared fixtures for the characterization tests.
#
# Expected outputs ("golden" files) were recorded from the implementation at
# commit ed80437 by running, from the package root:
#   Rscript tests/testthat/fixtures/make-golden.R
# Regenerate them only when a change in results is intended and explained.

# Deterministic 99-row subset: no random sampling, so it does not depend on
# how dplyr::slice_sample() uses the RNG.
fixture_data <- function(weighted = FALSE) {
  d <- spatInfer::opportunity[seq(1, 693, by = 7), ]
  if (weighted) d$weights <- 0.5 + (seq_len(nrow(d)) %% 4) / 2
  d
}

fixture_fm <- mobility ~ single_mothers + short_commute + gini + dropout_rate +
  social_cap + dropout_na

# Small settings to keep each simulation call to a couple of seconds.
fixture_args <- list(splines = 4, pc_num = 4, nSim = 20, max_clus = 5)

# fixest notes on collinear PCs within clusters, and non-PSD VCOV warnings,
# are expected with small clusters; silence them so test output stays readable.
quietly <- function(expr) suppressWarnings(suppressMessages(expr))

run_sim <- function(fun, ..., data = fixture_data()) {
  args <- utils::modifyList(c(list(fm = fixture_fm, df = data), fixture_args), list(...))
  quietly(do.call(fun, args))
}

summarise_basis <- function(ob) {
  list(
    title = ob$patches$annotation$title,
    subtitle = ob$patches$annotation$subtitle,
    bic = ob[[1]]$data,
    r2 = ob[[2]]$data
  )
}

summarise_feols <- function(m) {
  list(coeftable = m$coeftable, nobs = stats::nobs(m), r2 = fixest::r2(m, "r2"))
}

# Each case returns a plain R object that is compared with its golden file.
golden_cases <- list(
  optimal_basis = function() {
    summarise_basis(quietly(optimal_basis(fixture_fm, fixture_data(), max_splines = 5)))
  },
  # Known issue: with max_splines = 3 the loop 2:mx runs backwards, adding a
  # 4x4 tensor and duplicating the 3x3 results (BIC_3.x, BIC_4, BIC_3.y).
  optimal_basis_max3 = function() {
    summarise_basis(quietly(optimal_basis(fixture_fm, fixture_data(), max_splines = 3)))
  },
  basis_regression_bch = function() {
    summarise_feols(quietly(basis_regression(fixture_fm, fixture_data(),
      splines = 4, pc_num = 4, clusters = 4, cov = "BCH")))
  },
  basis_regression_hc = function() {
    summarise_feols(quietly(basis_regression(fixture_fm, fixture_data(),
      splines = 4, pc_num = 4, clusters = 4, cov = "HC")))
  },
  placebo = function() run_sim(placebo, Parallel = FALSE),
  placebo_im = function() run_sim(placebo_im, Parallel = FALSE),
  synth = function() run_sim(synth, Parallel = FALSE),
  synth_im = function() run_sim(synth_im, Parallel = FALSE),
  placebo_weighted = function() {
    run_sim(placebo, Parallel = FALSE, weights = TRUE, data = fixture_data(TRUE))
  },
  synth_im_weighted = function() {
    run_sim(synth_im, Parallel = FALSE, weights = TRUE, data = fixture_data(TRUE))
  },
  placebo_clara = function() run_sim(placebo, Parallel = FALSE, k_medoids = FALSE),
  placebo_brisc = function() run_sim(placebo, Parallel = FALSE, exact_cholesky = FALSE)
)

golden_path <- function(name) {
  testthat::test_path("fixtures", "golden", paste0(name, ".rds"))
}

expect_golden <- function(name) {
  path <- golden_path(name)
  if (!file.exists(path)) testthat::skip(paste("No golden file for", name))
  testthat::expect_identical(golden_cases[[name]](), readRDS(path))
}
