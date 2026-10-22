# SCPC inference (Mueller and Watson) in placebo() and synth().

test_that("SCPC rows follow the HC, BCH and IM rows, with avc in Clusters", {
  plbo <- readRDS(golden_path("placebo"))
  n_clus <- fixture_args$max_clus - 2
  expect_equal(plbo$Results$SE, c("HC", rep("BCH", n_clus), rep("IM", n_clus), rep("SCPC", 5)))
  expect_equal(tail(plbo$Results$Clusters, 5), as.character(fixture_avc))

  syn <- readRDS(golden_path("synth"))
  expect_equal(syn$Results$SE, plbo$Results$SE)
  expect_equal(syn$Results$Clusters, plbo$Results$Clusters)

  # Without avc there are no SCPC rows.
  expect_false("SCPC" %in% readRDS(golden_path("placebo_clara"))$Results$SE)
})

# The SCPC code is ported from scpcR 0.1.3. scpcR_reference.rds holds the results that
# scpcR::scpc() gave for the regressions below, so that they can be checked without scpcR.
scpcR_reference <- function() readRDS(golden_path("scpcR_reference"))

test_that("SCPC estimates match scpcR on the basis regression", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eq <- build_formulas(prep$rhs, prep$pc, "explan_var")$eq_est
  fit <- quietly(fixest::feols(eq, data = prep$df, weights = ~wts))
  est <- scpc_stats(fit, prep$df, fixture_avc)
  expect_equal(unname(est), scpcR_reference()$basis, tolerance = 1e-10)

  plbo <- readRDS(golden_path("placebo"))
  scpc <- plbo$Results$SE == "SCPC"
  expect_equal(plbo$Results$`Est p`[scpc], round(est[, "p"], 3), ignore_attr = TRUE)
  fmt <- function(lo, hi) paste0("[", signif(lo, 3), ", ", signif(hi, 3), "]")
  expect_identical(plbo$Results$CI[scpc], fmt(est[, "lower"], est[, "upper"]))
})

test_that("avc values are checked and sorted", {
  for (fun in list(placebo, synth)) {
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, avc = 0), "between 0.001 and 0.99")
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, avc = c(0.02, 1)), "between 0.001 and 0.99")
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, avc = "0.02"), "between 0.001 and 0.99")
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, avc = c(0.02, 0.02)), "repeated")
  }
  expect_identical(check_avc(c(0.1, 0.02, 0.05)), c(0.02, 0.05, 0.1))
})

test_that("simulations where SCPC fails are left out, with a warning", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  avc <- c(0.02, 0.1)
  sim_out <- data.frame(hc = rep(0.5, 4), scpc_0.02 = c(0.01, NA, 0.5, 0.5), scpc_0.1 = c(0.01, NA, 0.5, 0.9))
  warnings <- capture_warnings(res <- quietly_notes(
    ols_summary(prep$df, eqs$eq_est, list(), sim_out, avc, scpc_setup(prep$df, avc))))
  expect_match(warnings, "SCPC could not be computed in 1 of 4 simulations", all = FALSE)
  expect_equal(res$sim_05[-1], round(c(1, 1) / 3, 3))

  # A failing simulation gives NA SCPC p values instead of stopping the run.
  local_mocked_bindings(scpc_stats = function(...) stop("SCPC failed"))
  set.seed(1)
  Sim <- matrix(rnorm(nrow(prep$df)), ncol = 1)
  out <- quietly(ols_sim(1, Sim, eqs$eq_sim, prep$df, list(), avc))
  expect_named(out, c("hc", "scpc_0.02", "scpc_0.1"))
  expect_identical(unname(unlist(out[-1])), c(NA_real_, NA_real_))
})

test_that("a dropped treatment is an error in SCPC, not the next coefficient", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  # A constant placebo is dropped by feols as collinear.
  fit <- quietly(fixest::feols(eqs$eq_sim, data = cbind.data.frame(sim = 0, prep$df), weights = ~wts))
  expect_error(scpc_stats(fit, prep$df, c(0.02, 0.1)), "dropped from the regression")
})

test_that("tables with SCPC rows describe avc", {
  plac_tab <- placebo_table(readRDS(golden_path("placebo")))
  expect_identical(names(plac_tab@data)[2], "Clusters/avc")
  expect_identical(plac_tab@data$Adj[plac_tab@data$Adj != ""], c("HC", "BCH", "IM", "SCPC"))
  expect_match(unlist(plac_tab@notes), "average correlation bound avc", all = FALSE)

  syn_tab <- synth_table(readRDS(golden_path("synth")))
  expect_identical(names(syn_tab@data),
    c("Adj", "Clusters/avc", "Est p", "Synth p", "CI Width", "CI"))
  expect_match(unlist(syn_tab@notes), "average correlation bound avc", all = FALSE)

  # Tables without SCPC rows do not mention avc.
  clara_tab <- placebo_table(readRDS(golden_path("placebo_clara")))
  expect_identical(names(clara_tab@data)[2], "Clusters")
  expect_no_match(unlist(clara_tab@notes), "avc")
})

test_that("SCPC functions do not change the user's random number state", {
  set.seed(42)
  before <- .Random.seed
  run_sim(placebo, Parallel = FALSE, nSim = 5, avc = 0.02)
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
  sim_p <- quietly(ols_sim(1, matrix(df1$sim), eqs$eq_sim, prep$df, list(), avc))
  expect_equal(unname(unlist(sim_p[c("scpc_0.02", "scpc_0.1")])), direct[, 4], tolerance = 1e-10)

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
