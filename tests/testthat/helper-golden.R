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
  if (weighted) d$wt <- 0.5 + (seq_len(nrow(d)) %% 4) / 2
  d
}

fixture_fm <- mobility ~ single_mothers + short_commute + gini + dropout_rate +
  social_cap + dropout_na

# Small settings to keep each simulation call to a couple of seconds.
fixture_args <- list(splines = 4, pc_num = 4, nSim = 20, max_clus = 5)

# fixest notes on collinear PCs within clusters, and non-PSD VCOV warnings,
# are expected with small clusters; silence them so test output stays readable.
quietly <- function(expr) suppressWarnings(suppressMessages(expr))

# Silence fixest notes only, so that the package's own warnings can be tested.
quietly_notes <- function(expr) suppressMessages(expr)

run_sim <- function(fun, ..., data = fixture_data()) {
  args <- utils::modifyList(c(list(fm = fixture_fm, df = data), fixture_args), list(...))
  quietly(do.call(fun, args))
}

# SCPC bounds for the least squares golden cases.
fixture_avc <- seq(0.02, 0.1, by = 0.02)

# IV fixture: deterministic 121-row subset of the 2000-2007 period of autor_china.
iv_fixture_data <- function() {
  d <- spatInfer::autor_china[spatInfer::autor_china$t2 == 1, ]
  d <- d[seq(1, nrow(d), by = 6), ]
  d
}

iv_fixture_fm <- d_sh_empl_mfg ~ l_shind_manuf_cbp + l_sh_popedu_c | d_tradeusch_pw ~ d_tradeotch_pw_lag

# SCPC bounds for the Anderson-Rubin tests of the IV golden cases.
iv_avc <- c(0.05, 0.1)

run_iv <- function(fun, ..., data = iv_fixture_data()) {
  args <- utils::modifyList(c(list(fm = iv_fixture_fm, df = data), fixture_args), list(...))
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

# The coefficient table is kept as a plain matrix: fixest labels it with attributes
# whose names change between versions ('type' before 0.14, 'vcov_type' after).
plain_coeftable <- function(ct) {
  ct[, , drop = FALSE]
}

summarise_feols <- function(m) {
  list(coeftable = plain_coeftable(m$coeftable), nobs = stats::nobs(m), r2 = fixest::r2(m, "r2"))
}

# Each case returns a plain R object that is compared with its golden file.
golden_cases <- list(
  optimal_basis = function() {
    summarise_basis(quietly(optimal_basis(fixture_fm, fixture_data(), max_splines = 5)))
  },
  optimal_basis_max3 = function() {
    summarise_basis(quietly(optimal_basis(fixture_fm, fixture_data(), max_splines = 3)))
  },
  basis_regression_bch = function() {
    summarise_feols(quietly(bch_regression(fixture_fm, fixture_data(),
      splines = 4, pc_num = 4, clusters = 4, cov = "BCH")))
  },
  basis_regression_hc = function() {
    summarise_feols(quietly(bch_regression(fixture_fm, fixture_data(),
      splines = 4, pc_num = 4, clusters = 4, cov = "HC")))
  },
  basis_regression_weighted = function() {
    summarise_feols(quietly(bch_regression(fixture_fm, fixture_data(TRUE),
      splines = 4, pc_num = 4, clusters = 4, weights = "wt")))
  },
  # Least squares placebo and synth cases were recorded when the BCH, IM and SCPC versions were
  # combined into placebo() and synth(); their rows matched the earlier separate versions exactly.
  placebo = function() run_sim(placebo, Parallel = FALSE, avc = fixture_avc),
  synth = function() run_sim(synth, Parallel = FALSE, avc = fixture_avc),
  # The weighted cases were recorded when weighted conditional SCPC rescaled the regressors by the
  # square root of the weights, now the option scpc_weighting = "sqrt", which they check.
  placebo_weighted = function() {
    run_sim(placebo, Parallel = FALSE, avc = fixture_avc, weights = "wt", data = fixture_data(TRUE),
            scpc_weighting = "sqrt")
  },
  synth_weighted = function() {
    run_sim(synth, Parallel = FALSE, avc = fixture_avc, weights = "wt", data = fixture_data(TRUE),
            scpc_weighting = "sqrt")
  },
  # Weighted cases with the default scpc_weighting = "scpcR", recorded when it became the default.
  placebo_weighted_scpcr = function() {
    run_sim(placebo, Parallel = FALSE, avc = fixture_avc, weights = "wt", data = fixture_data(TRUE))
  },
  synth_weighted_scpcr = function() {
    run_sim(synth, Parallel = FALSE, avc = fixture_avc, weights = "wt", data = fixture_data(TRUE))
  },
  placebo_clara = function() run_sim(placebo, Parallel = FALSE, k_medoids = FALSE),
  placebo_brisc = function() run_sim(placebo, Parallel = FALSE, exact_cholesky = FALSE),
  # IV cases were recorded when the Anderson-Rubin tests were redone by the reduced form
  # approach, with IM and SCPC inference added.
  iv_placebo = function() run_iv(iv_placebo, Parallel = FALSE, avc = iv_avc),
  iv_synth = function() run_iv(iv_synth, Parallel = FALSE, avc = iv_avc),
  iv_placebo_weighted = function() {
    run_iv(iv_placebo, Parallel = FALSE, weights = "timepwt48")
  },
  iv_placebo_weighted_scpc = function() {
    run_iv(iv_placebo, Parallel = FALSE, weights = "timepwt48", avc = iv_avc)
  }
)

golden_path <- function(name) {
  testthat::test_path("fixtures", "golden", paste0(name, ".rds"))
}

# Golden comparisons allow for rounding error of about 1e-10 (relative). Updating
# a dependency can change the last digits of results without any change in the
# method: fixest 0.14.2, for example, differs from 0.13.2 by about 1e-14 (relative).
# Any change in a method moves results far more than this.
golden_tolerance <- 1e-10

expect_matches_golden <- function(object, name) {
  testthat::expect_equal(object, readRDS(golden_path(name)), tolerance = golden_tolerance)
}

expect_golden <- function(name) {
  path <- golden_path(name)
  if (!file.exists(path)) testthat::skip(paste("No golden file for", name))
  expect_matches_golden(golden_cases[[name]](), name)
}

# With 99 observations, the 5x5 tensor with 23 or more PCs is nearly saturated, and
# which near-collinear columns are dropped depends on the BLAS and package versions:
# its BIC and R2 have moved by about 0.5% between installations (316.37, 316.77 and
# 317.22 for the BIC with 23 PCs) while every other value agreed to 1e-12. Those rows
# are compared with a 1% tolerance and everything else, including the chosen basis,
# with golden_tolerance.
saturated_tolerance <- 0.01

expect_basis_golden <- function(name, saturated_splines = 5, saturated_pcs = 23) {
  path <- golden_path(name)
  if (!file.exists(path)) testthat::skip(paste("No golden file for", name))
  object <- golden_cases[[name]]()
  expected <- readRDS(path)
  for (curve in c("bic", "r2")) {
    saturated <- object[[curve]]$index >= saturated_pcs &
      grepl(paste0("_", saturated_splines, "$"), object[[curve]]$name)
    testthat::expect_identical(saturated,
      expected[[curve]]$index >= saturated_pcs &
        grepl(paste0("_", saturated_splines, "$"), expected[[curve]]$name))
    testthat::expect_equal(object[[curve]][saturated, ], expected[[curve]][saturated, ],
                           tolerance = saturated_tolerance)
    object[[curve]] <- object[[curve]][!saturated, ]
    expected[[curve]] <- expected[[curve]][!saturated, ]
  }
  testthat::expect_equal(object, expected, tolerance = golden_tolerance)
}
