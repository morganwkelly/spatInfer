# SCPC inference (Mueller and Watson) in placebo_scpc() and synth_scpc().

test_that("placebo_scpc() and synth_scpc() results are unchanged", {
  expect_golden("placebo_scpc")
  expect_golden("synth_scpc")
  expect_golden("placebo_scpc_weighted")
})

default_avc <- seq(0.02, 0.1, by = 0.02)

test_that("SCPC results have the documented components and columns", {
  plbo <- readRDS(golden_path("placebo_scpc"))
  expect_null(attr(plbo, "class"))
  expect_named(plbo, c("Results", "Spatial_Params", "Estimates"))
  expect_named(plbo$Results,
    c("SE", "avc", "Est p", "Plac p", "Plac 5%", "CI Width", "CI", "Pseudo SE"))
  expect_equal(plbo$Results$SE, c("HC", rep("SCPC", 5)))
  expect_equal(plbo$Results$avc, c(".", as.character(default_avc)))
  expect_named(plbo$Estimates, c("avc", "Coef", "Std_Err", "t", "p", "lower", "upper"))
  expect_equal(plbo$Estimates$avc, default_avc)

  syn <- readRDS(golden_path("synth_scpc"))
  expect_s3_class(syn, "synth_scpc")
  expect_named(syn$Results,
    c("SE", "avc", "Est p", "Synth p", "CI Width", "CI", "Pseudo SE"))
  expect_named(syn$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "N", "Splines", "PCs"))
})

# The SCPC code is ported from scpcR 0.1.3. scpcR_reference.rds holds the results that
# scpcR::scpc() gave for the regressions below, so that they can be checked without scpcR.
scpcR_reference <- function() readRDS(golden_path("scpcR_reference"))

test_that("SCPC estimates match scpcR on the basis regression", {
  plbo <- readRDS(golden_path("placebo_scpc"))
  expect_equal(unname(as.matrix(plbo$Estimates[, -1])), scpcR_reference()$basis, tolerance = 1e-10)
  expect_equal(plbo$Results$`Est p`[-1], round(plbo$Estimates$p, 3))
})

test_that("Pseudo SE is one quarter of the width of the 95% confidence interval", {
  plbo <- readRDS(golden_path("placebo_scpc"))
  est <- plbo$Estimates
  expect_equal(plbo$Results$`Pseudo SE`[-1], signif((est$upper - est$lower) / 4, 3))

  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eq <- build_formulas(prep$rhs, prep$pc, "explan_var")$eq_est
  hc <- quietly(fixest::feols(eq, data = prep$df, weights = ~wts, vcov = "hetero"))
  expect_equal(plbo$Results$`Pseudo SE`[1], signif((confint(hc)[2, 2] - confint(hc)[2, 1]) / 4, 3))
})

test_that("avc values are checked and sorted", {
  for (fun in list(placebo_scpc, synth_scpc)) {
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, avc = 0), "between 0.001 and 0.99")
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, avc = c(0.02, 1)), "between 0.001 and 0.99")
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, avc = "0.02"), "between 0.001 and 0.99")
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, avc = c(0.02, 0.02)), "repeated")
  }
  expect_identical(check_avc(c(0.1, 0.02, 0.05)), c(0.02, 0.05, 0.1))
})

test_that("SCPC functions check coordinates and weights like the others", {
  d <- fixture_data()
  no_x <- d[, setdiff(names(d), "X")]
  for (fun in list(placebo_scpc, synth_scpc)) {
    expect_error(fun(fixture_fm, no_x, 4, 4), "named X and Y")
    expect_error(fun(fixture_fm, d, 4, 4, weights = "pop"), "no variable called pop")
  }
})

test_that("simulations where SCPC fails are left out, with a warning", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  avc <- c(0.02, 0.1)
  scpc_out <- data.frame(avc_0.02 = c(0.01, NA, 0.5, 0.5), avc_0.1 = c(0.01, NA, 0.5, 0.9))
  hc_out <- data.frame(hc_p = rep(0.5, 4))
  warnings <- capture_warnings(summ <- quietly_notes(summary_scpc(prep$df, eqs$eq_est, avc, scpc_out, hc_out)))
  expect_match(warnings, "SCPC could not be computed in 1 of 4 simulations", all = FALSE)
  expect_equal(summ$Results$sim_05[-1], round(c(1, 1) / 3, 3))

  # A failing simulation gives NA p values instead of stopping the run. A constant
  # placebo is dropped by feols as collinear, so it must not be replaced by the
  # next coefficient.
  Sim <- matrix(0, nrow = nrow(prep$df), ncol = 1)
  expect_identical(unname(unlist(quietly(scpc_sim(1, Sim, eqs$eq_sim, prep$df, avc)))), c(NA_real_, NA_real_))
  fit <- quietly(fixest::feols(eqs$eq_sim, data = cbind.data.frame(sim = Sim[, 1], prep$df), weights = ~wts))
  expect_error(scpc_stats(fit, prep$df, avc), "dropped from the regression")
})

test_that("tables of SCPC results leave out pseudo SEs and describe avc", {
  plac_tab <- placebo_table(readRDS(golden_path("placebo_scpc")))
  expect_identical(names(plac_tab@data),
    c("Adj", "avc", "Est p", "Plac p", "Plac 5%", "CI Width", "CI"))
  expect_identical(plac_tab@data$Adj, c("HC", "SCPC", rep("", 4)))
  expect_no_match(unlist(plac_tab@notes), "Pseudo SE")
  expect_match(unlist(plac_tab@notes), "average correlation bound avc", all = FALSE)

  syn_tab <- synth_table(readRDS(golden_path("synth_scpc")))
  expect_identical(names(syn_tab@data),
    c("Adj", "avc", "Est p", "Synth p", "CI Width", "CI"))
  expect_match(unlist(syn_tab@notes), "average correlation bound avc", all = FALSE)

  # Tables of cluster results keep their original note.
  expect_no_match(unlist(synth_table(readRDS(golden_path("synth")))@notes), "avc")
})

test_that("SCPC functions do not change the user's random number state", {
  set.seed(42)
  before <- .Random.seed
  run_scpc(placebo_scpc, Parallel = FALSE, nSim = 5)
  expect_identical(.Random.seed, before)
})

test_that("SCPC with saved spatial weights matches scpcR in simulated regressions", {
  # This regression is weighted, where spatInfer departs from scpcR only through
  # scpc_cond_weights(); with scpcR's choice everything else must match.
  local_mocked_bindings(scpc_cond_weights = function(w) rep(1, length(w)))
  prep <- prepare_spatial_data(fixture_fm, fixture_data(TRUE), 4, 4, "wt")
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  set.seed(1)
  df1 <- cbind.data.frame(sim = rnorm(nrow(prep$df)), prep$df)
  fit <- quietly(fixest::feols(eqs$eq_sim, data = df1, weights = ~wts))
  treat <- which(names(coef(fit)) == "sim")
  avc <- c(0.02, 0.1)
  direct <- scpcR_reference()$sim_exact
  expect_equal(unname(scpc_stats(fit, df1, avc)), unname(direct), tolerance = 1e-10)
  expect_equal(unname(unlist(scpc_sim(1, matrix(df1$sim), eqs$eq_sim, prep$df, avc))), direct[, 4],
               tolerance = 1e-10)

  # The approximate method, used by scpcR above 4500 observations.
  spc <- scpc_spatial(as.matrix(df1[, c("X", "Y")]), 0.04, large_n = TRUE)
  direct <- scpcR_reference()$approx
  expect_equal(unname(scpc_coef_stats(fit, spc, treat)), direct[treat, , drop = FALSE],
               tolerance = 1e-10)
  expect_equal(unname(scpc_coef_stats(fit, spc, seq_along(coef(fit)))), direct, tolerance = 1e-10)
})

test_that("weighted conditional SCPC rescales the regressors by the square root of the weights", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(TRUE), 4, 4, "wt")
  eq <- build_formulas(prep$rhs, prep$pc, "explan_var")$eq_est
  avc <- c(0.02, 0.1)
  fit_w <- quietly(fixest::feols(eq, data = prep$df, weights = ~wts))
  weighted <- scpc_stats(fit_w, prep$df, avc)
  expect_false(isTRUE(all.equal(weighted, scpcR_reference()$sim_exact)))

  # The rescaling is applied: replacing it with scpcR's gives scpcR's results instead.
  local({
    local_mocked_bindings(scpc_cond_weights = function(w) rep(1, length(w)))
    expect_false(isTRUE(all.equal(scpc_stats(fit_w, prep$df, avc), weighted)))
  })

  # Equal weights give the unweighted results.
  df2 <- prep$df
  df2$wts <- 2
  fit_2 <- quietly(fixest::feols(eq, data = df2, weights = ~wts))
  fit_1 <- quietly(fixest::feols(eq, data = df2))
  expect_equal(scpc_stats(fit_2, df2, avc), scpc_stats(fit_1, df2, avc), tolerance = 1e-10)
})

test_that("Gauss-Legendre quadrature integrates polynomials exactly", {
  expect_equal(sum(scpc_gq$w), 1)
  expect_equal(sum(scpc_gq$w * scpc_gq$x^79), 1 / 80)
})
