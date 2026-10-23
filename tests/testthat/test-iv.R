# Instrumental variable placebo and synthetic outcome tests: iv_placebo() and iv_synth().

iv_se_labels <- function(n_clus, n_avc = 0) {
  c("HC", rep("BCH", n_clus), rep("IM", n_clus), rep("SCPC", n_avc),
    "AR HC", rep("AR BCH", n_clus), rep("AR IM", n_clus), rep("AR SCPC", n_avc))
}

iv_clusters_col <- function(avc = NULL) {
  clus <- as.character(3:fixture_args$max_clus)
  c(".", clus, clus, as.character(avc), ".", clus, clus, as.character(avc))
}

test_that("iv_placebo() and iv_synth() results are unchanged", {
  expect_golden("iv_placebo")
  expect_golden("iv_synth")
  expect_golden("iv_placebo_weighted")
  expect_golden("iv_placebo_weighted_scpc")
})

test_that("IV results have the documented components and columns", {
  n_clus <- fixture_args$max_clus - 2
  plbo <- readRDS(golden_path("iv_placebo"))
  expect_null(attr(plbo, "class"))
  expect_named(plbo, c("Results", "Spatial_Params"))
  expect_named(plbo$Results,
    c("SE", "Clusters", "F", "Est p", "Plac p", "Plac 5%", "CI Width", "CI"))
  expect_equal(plbo$Results$SE, iv_se_labels(n_clus, length(iv_avc)))
  expect_equal(plbo$Results$Clusters, iv_clusters_col(iv_avc))
  weighted <- readRDS(golden_path("iv_placebo_weighted"))$Results
  expect_equal(weighted$SE, iv_se_labels(n_clus))   #no SCPC rows without avc
  expect_equal(weighted$Clusters, iv_clusters_col())
  expect_named(plbo$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "R2_Inst", "Effective_Range_Inst",
      "Structure_Inst", "N", "Splines", "PCs"))

  syn <- readRDS(golden_path("iv_synth"))
  expect_null(attr(syn, "class"))
  expect_named(syn$Results, c("SE", "Clusters", "F", "Est p", "Synth p", "CI Width", "CI"))
  expect_equal(syn$Results$SE, iv_se_labels(n_clus, length(iv_avc)))
  expect_named(syn$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "N", "Splines", "PCs"))

  # Estimates with the real data do not depend on the simulations.
  shared <- c("SE", "Clusters", "F", "Est p", "CI Width", "CI")
  expect_identical(plbo$Results[shared], syn$Results[shared])
})

# F must follow the inference and number of clusters of each row: HC and BCH from the first stage
# t statistic of the instrument. IM rows have no F. AR HC and BCH rows share the F of their inference;
# AR IM and SCPC rows have none.
test_that("first stage F statistics follow the inference and number of clusters", {
  res <- readRDS(golden_path("iv_placebo"))$Results
  prep <- quietly(prepare_iv_data(iv_fixture_fm, iv_fixture_data(), 4, 4, NULL, 5))
  df <- prep$df
  eqs <- build_iv_formulas(prep$rhs, prep$pc)
  fs <- quietly(fixest::feols(eqs$eq_fs, df, weights = ~wts, vcov = "hetero"))
  fit <- quietly(fixest::feols(eqs$eq_est, df, weights = ~wts, vcov = "hetero"))
  hold_clus <- quietly(sim_clusters(df, TRUE, 5))
  expect_equal(res$F[res$SE == "HC"], round(fixest::fitstat(fit, "ivwald1", simplify = TRUE)$stat, 3))
  for (clust in hold_clus) {
    G <- as.character(nlevels(clust))
    bch_fit <- quietly(fixest::feols(eqs$eq_est, df, weights = ~wts, cluster = clust))
    expect_equal(res$F[res$SE == "BCH" & res$Clusters == G],
                 round(fixest::fitstat(bch_fit, "ivwald1", simplify = TRUE)$stat, 3))
  }
  expect_true(all(is.na(res$F[res$SE %in% c("IM", "SCPC", "AR IM", "AR SCPC")])))
  for (se in c("HC", "BCH"))
    expect_identical(res$F[res$SE == paste("AR", se)], res$F[res$SE == se])
  expect_gt(length(unique(stats::na.omit(res$F))), 1)
})

test_that("IV formulas are parsed into outcome, controls, treatment and instrument", {
  ivf <- parse_iv_formula(y ~ a + b | x ~ z)
  expect_identical(ivf, list(outcome = "y", treatment = "x", instrument = "z", rhs = "+a+b"))
  expect_identical(parse_iv_formula(y ~ 1 | x ~ z)$rhs, "")
  expect_identical(deparse1(iv_ols_formula(y ~ a + b | x ~ z)), "y ~ x + a + b")
})

test_that("malformed IV formulas and data are rejected", {
  d <- iv_fixture_data()
  expect_error(iv_placebo(d_sh_empl_mfg ~ d_tradeusch_pw, d, 4, 4), "must have the form")
  expect_error(iv_synth(d_sh_empl_mfg ~ l_sh_popedu_c | d_tradeusch_pw ~ d_tradeotch_pw_lag + yr, d, 4, 4),
    "one instrument")
  expect_error(iv_placebo(d_sh_empl_mfg ~ d_tradeotch_pw_lag | d_tradeusch_pw ~ d_tradeotch_pw_lag, d, 4, 4),
    "cannot also be controls")
  expect_error(iv_placebo(d_sh_empl_mfg ~ 1 | d_tradeusch_pw ~ not_a_variable, d, 4, 4), "not_a_variable")
  na_d <- d
  na_d$d_tradeotch_pw_lag[1] <- NA
  for (fun in list(iv_placebo, iv_synth)) {
    expect_error(fun(iv_fixture_fm, na_d, 4, 4), "missing values")
    expect_error(fun(iv_fixture_fm, d[, setdiff(names(d), "X")], 4, 4), "named X and Y")
    expect_error(fun(iv_fixture_fm, d, 4, 4, max_clus = 2), "greater than 2")
    expect_error(fun(iv_fixture_fm, d, 4, 4, weights = "pop"), "no variable called pop")
  }
})

# The AR confidence sets are found analytically from the reduced forms. At their finite endpoints a
# direct test of the instrument in a regression of outcome - b*treatment, with the same inference,
# must give a p value of 0.05, and at b = 0 the p value of no effect. SCPC critical values are found by
# bisection, so its endpoints are checked to 1e-3. The full 2000-2007 data have a strong instrument
# and bounded HC and BCH sets; the weak instrument of the fixture gives rays.
ar_endpoint_p <- function(data, splines, pc_num, clust_index, avc) {
  prep <- quietly(prepare_iv_data(iv_fixture_fm, data, splines, pc_num, "timepwt48", 6))
  df <- prep$df
  eqs <- build_iv_formulas(prep$rhs, prep$pc)
  eq_diff <- stats::update(eqs$eq_rf, ar_diff ~ .)
  clust <- quietly(sim_clusters(df, TRUE, 6))[[clust_index]]
  spc <- scpc_setup(df, avc)[[1]]
  rf <- ar_data_rf(df, eqs)
  stats <- list(
    HC = ar_hc(rf),
    BCH = ar_bch(rf, clust),
    IM = quietly(ar_im(eqs$eq_rf, df, "inst_var", cbind(y = df$dep_var, x = df$explan_var), clust)),
    SCPC = ar_scpc(rf, spc))
  p_at <- function(b, inference) {
    df$ar_diff <- df$dep_var - b * df$explan_var
    fit <- quietly(fixest::feols(eq_diff, df, weights = ~wts, vcov = "hetero"))
    switch(inference,
      HC = fixest::pvalue(fit)[["inst_var"]],
      BCH = fixest::pvalue(fit, cluster = clust)[["inst_var"]],
      IM = im_t_test(list(estimate = quietly(iv_cluster_coefs(eq_diff, df, clust, "inst_var"))))$p.value,
      SCPC = scpc_coef_stats(fit, spc, match("inst_var", names(coef(fit))))[, "p"])
  }
  lapply(stats::setNames(names(stats), names(stats)), function(inference) {
    set <- ar_set(stats[[inference]])
    list(set = set, p = vapply(set[is.finite(set)], p_at, numeric(1), inference = inference),
         p0 = stats[[inference]]$p, direct_p0 = p_at(0, inference))
  })
}

expect_ar_endpoints <- function(res) {
  for (inference in names(res)) {
    r <- res[[inference]]
    tol <- if (inference == "SCPC") 1e-3 else 1e-8
    expect_true(all(abs(r$p - 0.05) < tol), label = paste(inference, "endpoint p values"))
    expect_equal(r$p0, r$direct_p0, tolerance = 1e-8, label = paste(inference, "p value of no effect"))
  }
}

test_that("AR confidence sets have p values of 0.05 at their endpoints for every inference", {
  strong <- spatInfer::autor_china[spatInfer::autor_china$t2 == 1, ]
  res <- ar_endpoint_p(strong, 6, 8, 2, 0.05)
  for (inference in c("HC", "BCH")) expect_length(res[[inference]]$set, 2)
  expect_ar_endpoints(res)
  weak <- ar_endpoint_p(iv_fixture_data(), 4, 4, 1, 0.05)
  expect_true(any(vapply(weak, function(r) length(r$set) == 4, logical(1))))
  expect_ar_endpoints(weak)
})

# The variances of the reduced form coefficients are those fixest gives for HC and BCH, and those
# of scpc_coef_stats() for SCPC, for the outcome, the treatment and any combination of the two.
test_that("reduced form variances agree with direct regressions", {
  prep <- quietly(prepare_iv_data(iv_fixture_fm, iv_fixture_data(), 4, 4, "timepwt48", 5))
  df <- prep$df
  eqs <- build_iv_formulas(prep$rhs, prep$pc)
  clust <- quietly(sim_clusters(df, TRUE, 5))[[2]]
  spc <- scpc_setup(df, 0.05)[[1]]
  rf <- ar_data_rf(df, eqs)
  df$comb <- df$dep_var - 0.7 * df$explan_var
  eq_comb <- stats::update(eqs$eq_rf, comb ~ .)
  fit <- quietly(fixest::feols(eq_comb, df, weights = ~wts))
  v <- c(1, -0.7)
  expect_equal(drop(v %*% ar_hc(rf)$V %*% v), fixest::se(fit, vcov = "hetero")[["inst_var"]]^2, tolerance = 1e-10)
  expect_equal(drop(v %*% ar_bch(rf, clust)$V %*% v), fixest::se(fit, cluster = clust)[["inst_var"]]^2,
               tolerance = 1e-10)
  sc <- scpc_coef_stats(fit, spc, match("inst_var", names(coef(fit))))
  expect_equal(drop(v %*% ar_scpc(rf, spc)$V %*% v), sc[, "Std_Err"]^2, tolerance = 1e-10)
  expect_equal(drop(v %*% rf$coef), unname(coef(fit)["inst_var"]), tolerance = 1e-10)
  fs <- quietly(fixest::feols(eqs$eq_fs, df, weights = ~wts))
  expect_equal(ar_F(ar_bch(rf, clust)), fixest::tstat(fs, cluster = clust)[["inst_var"]]^2, tolerance = 1e-10)
})

test_that("ar_conf_set() gives intervals, pairs of rays and the whole line", {
  # (a - b c)^2 <= crit * (Vyy - 2 b Vxy + b^2 Vxx) checked on a grid
  in_set <- function(b, a, c, Vyy, Vxy, Vxx, crit) (a - b * c)^2 <= crit * (Vyy - 2 * b * Vxy + b^2 * Vxx)
  cases <- list(
    interval = c(a = 2, c = 1, Vyy = 0.1, Vxy = 0.02, Vxx = 0.05, crit = 4),
    rays = c(a = 2, c = 0.3, Vyy = 0.1, Vxy = 0.02, Vxx = 0.05, crit = 4),
    line = c(a = 0.2, c = 0.1, Vyy = 0.1, Vxy = 0.02, Vxx = 0.05, crit = 4))
  grid <- seq(-200, 200, by = 0.01)
  for (name in names(cases)) {
    k <- as.list(cases[[name]])
    set <- do.call(ar_conf_set, k)
    inside <- do.call(in_set, c(list(b = grid), k))
    expected <- switch(name,
      interval = grid >= set[1] & grid <= set[2],
      rays = grid <= set[2] | grid >= set[3],
      line = rep(TRUE, length(grid)))
    expect_length(set, switch(name, interval = 2, rays = 4, line = 2))
    expect_identical(inside, expected)
  }
  expect_identical(format_conf_set(c(-Inf, -1, 2, Inf)), "(-Inf, -1] U [2, Inf)")
  expect_identical(format_conf_set(c(-Inf, Inf)), "(-Inf, Inf)")
  expect_identical(format_conf_set(c(-1.23456, 2)), "[-1.23, 2]")
})

test_that("unbounded AR sets give an infinite width", {
  row <- conf_set_row("AR HC", 0, 0.3, c(-Inf, -1, 2, Inf), c(0.1, 0.5))
  expect_identical(row$width_ci, Inf)
  expect_identical(row$CI, "(-Inf, -1] U [2, Inf)")
  expect_equal(row$sim_p, 0.5)
})

test_that("optimal_basis() and plot_basis() accept IV formulas", {
  ols_fm <- d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp + l_sh_popedu_c
  expect_equal(summarise_basis(quietly(optimal_basis(iv_fixture_fm, iv_fixture_data(), max_splines = 4))),
               summarise_basis(quietly(optimal_basis(ols_fm, iv_fixture_data(), max_splines = 4))))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(quietly(plot_basis(iv_fixture_fm, iv_fixture_data(), splines = 4)))
})

test_that("IV simulations leave the user's random number state unchanged", {
  for (fun in list(iv_placebo, iv_synth)) {
    set.seed(42)
    expected <- runif(3)
    set.seed(42)
    run_iv(fun, Parallel = FALSE, nSim = 5)
    expect_identical(runif(3), expected)
  }
})

test_that("placebo_table() and synth_table() show 2SLS and AR results in one table", {
  plbo <- readRDS(golden_path("iv_placebo"))
  n_clus <- fixture_args$max_clus - 2
  rows <- table_rows(plbo$Results)
  block <- c("HC", "BCH", rep("", n_clus - 1), "IM", rep("", n_clus - 1), "SCPC", rep("", length(iv_avc) - 1))
  expect_equal(rows$rows$SE, c(block, block))
  expect_equal(rows$n_wald, 1 + 2 * n_clus + length(iv_avc))
  is_im_scpc <- sub("^AR ", "", plbo$Results$SE) %in% c("IM", "SCPC")
  expect_identical(rows$rows$F[is_im_scpc], rep("", sum(is_im_scpc)))
  expect_false(any(rows$rows$F[!is_im_scpc] == ""))
  expect_identical(rows$rows$CI, plbo$Results$CI)

  tab <- placebo_table(plbo)
  expect_s4_class(tab, "tinytable")
  expect_identical(nrow(tab@data), nrow(plbo$Results))
  expect_identical(names(tab@data),
    c("Adj", "Clusters/avc", "F", "Est p", "Plac p", "Plac 5%", "CI Width", "CI"))
  notes <- unlist(tab@notes)
  expect_true(any(grepl("Noise instrument", notes)))
  expect_true(any(grepl("Anderson-Rubin", notes)))
  expect_true(any(grepl("not given for IM or SCPC", notes)))

  syn_tab <- synth_table(readRDS(golden_path("iv_synth")))
  expect_identical(nrow(syn_tab@data), nrow(tab@data))
  expect_false(any(grepl("Noise instrument", unlist(syn_tab@notes))))
  # The 2SLS and Anderson-Rubin headings are orange, as in placebo tables.
  expect_match(tinytable::save_tt(syn_tab, "html"), "#FFA500")

  # Without avc there are no SCPC rows or notes.
  weighted <- placebo_table(readRDS(golden_path("iv_placebo_weighted")))
  expect_identical(names(weighted@data)[2], "Clusters")
  expect_false(any(grepl("SCPC", unlist(weighted@notes))))
})

# iv_regression() and iv_im_regression() use the same clusters and tests as iv_placebo(), so their
# results for the real data must agree with the iv_placebo() rows for the same inference and clusters.
fit_iv <- function(cov, clusters = 4, data = iv_fixture_data()) {
  if (cov == "IM") return(quietly(iv_im_regression(iv_fixture_fm, data, 4, 4, clusters = clusters)))
  if (cov == "SCPC") return(quietly(iv_scpc_regression(iv_fixture_fm, data, 4, 4, avc = iv_avc[1])))
  quietly(iv_regression(iv_fixture_fm, data, 4, 4, clusters = if (cov == "HC") NULL else clusters, cov = cov))
}

test_that("iv_regression(), iv_im_regression() and iv_scpc_regression() agree with iv_placebo()", {
  plbo <- readRDS(golden_path("iv_placebo"))$Results
  row <- function(se, clusters) plbo[plbo$SE == se & plbo$Clusters == clusters, ]
  for (cov in c("BCH", "IM", "HC", "SCPC")) {
    reg <- fit_iv(cov)
    k <- switch(cov, HC = ".", SCPC = as.character(iv_avc[1]), "4")
    wald <- row(cov, k)
    treat <- reg$Coefs[reg$Coefs$term == "fit_d_tradeusch_pw", ]
    expect_equal(round(treat$p.value, 3), wald$`Est p`)
    expect_identical(paste0("[", signif(treat$conf.low, 3), ", ", signif(treat$conf.high, 3), "]"), wald$CI)
    expect_equal(round(reg$F, 3), wald$F)   #both NA for IM and SCPC
    ar <- row(paste("AR", cov), k)
    expect_equal(round(reg$AR$p.value, 3), ar$`Est p`)
    expect_identical(reg$AR$CI, ar$CI)
  }
})

test_that("iv_regression() and iv_im_regression() return the documented components", {
  bch <- fit_iv("BCH")
  expect_s3_class(bch, "iv_regression")
  expect_false(inherits(bch, "iv_im_regression"))
  expect_named(bch, c("Coefs", "AR", "AR_Set", "F", "Inference", "Clusters", "N", "Splines", "PCs", "Fit"))
  expect_identical(bch$Coefs$term, c("(Intercept)", "fit_d_tradeusch_pw", "l_shind_manuf_cbp", "l_sh_popedu_c"))
  expect_s3_class(bch$Fit, "fixest")
  expect_equal(bch$N, nrow(iv_fixture_data()))
  im <- fit_iv("IM")
  expect_s3_class(im, c("iv_im_regression", "iv_regression"), exact = TRUE)
  expect_named(im, c("Coefs", "AR", "AR_Set", "F", "Inference", "Clusters", "N", "Splines", "PCs",
                     "Cluster_Coefs", "Cluster_N"))
  expect_identical(dim(im$Cluster_Coefs), c(4L, 5L))
  expect_equal(sum(im$Cluster_N), nrow(iv_fixture_data()))
  expect_output(print(im), "IM inference with 4 clusters")
  expect_output(print(im), "Coefficients in each cluster")
  expect_output(print(bch), "Anderson-Rubin test")
  expect_output(print(im), "Anderson-Rubin test")
  for (type in c("coefs", "clusters")) {
    tab <- im_table(im, type = type)
    expect_s4_class(tab, "tinytable")
    expect_true(any(grepl("2SLS in each cluster", unlist(tab@notes))))
    expect_false(any(grepl("Anderson-Rubin", unlist(tab@notes))))
  }
})

test_that("iv_regression() and iv_im_regression() check their arguments", {
  d <- iv_fixture_data()
  expect_error(iv_regression(iv_fixture_fm, d, 4, 4), "whole number")
  expect_error(iv_regression(iv_fixture_fm, d, 4, 4, clusters = 4, cov = "IM"), "should be one of")
  expect_error(iv_im_regression(iv_fixture_fm, d, 4, 4, clusters = 1.5), "whole number")
})

test_that("optimal_basis() and plot_basis() accept IV formulas without controls", {
  fm <- d_sh_empl_mfg ~ 1 | d_tradeusch_pw ~ d_tradeotch_pw_lag
  expect_equal(summarise_basis(quietly(optimal_basis(fm, iv_fixture_data(), max_splines = 4))),
               summarise_basis(quietly(optimal_basis(d_sh_empl_mfg ~ d_tradeusch_pw, iv_fixture_data(),
                                                     max_splines = 4))))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(quietly(plot_basis(fm, iv_fixture_data(), splines = 4)))
})

test_that("iv_regression_table() gives coefficients, first stage F and AR results by column", {
  skip_if_not_installed("modelsummary")
  d <- iv_fixture_data()
  orig <- quietly(fixest::feols(iv_fixture_fm, d, vcov = "hetero"))
  bch <- fit_iv("BCH")
  im <- fit_iv("IM")
  tab <- iv_regression_table(orig = orig, bch = bch, im = im, output = "data.frame")
  expect_identical(names(tab)[4:6], c("Original", "BCH", "IM"))
  row <- function(term) unlist(tab[tab$term == term, 4:6])
  expect_equal(unname(row("First stage F")),
               c(format(signif(fixest::fitstat(orig, "ivwald1", simplify = TRUE)$stat, 3)),
                 format(signif(bch$F, 3)), ""))
  expect_true(is.na(im$F))
  expect_false(any(grepl("First stage F", capture.output(print(im)))))
  expect_identical(unname(row("AR p")),
                   c("", format(signif(bch$AR$p.value, 3)), format(signif(im$AR$p.value, 3))))
  expect_identical(unname(row("AR 95% set")), c("", bch$AR$CI, im$AR$CI))
  expect_true("fit_d_tradeusch_pw" %in% tab$term)
  expect_false(any(grepl("^PC[0-9]", tab$term)))
  hc <- iv_regression_table(bch = fit_iv("HC"), output = "data.frame")
  expect_identical(names(hc)[4], "HC")

  expect_error(iv_regression_table(), "at least one")
  expect_error(iv_regression_table(bch = im), "iv_regression()")
  expect_error(iv_regression_table(im = bch), "iv_im_regression()")
  expect_error(regression_table(bch = bch), "iv_regression_table")
  expect_error(regression_table(im = im), "iv_regression_table")
})

test_that("iv_regression_plot() draws AR sets, with arrows for unbounded parts", {
  bch <- quietly(iv_regression(iv_fixture_fm, iv_fixture_data(), 4, 4, clusters = 3, weights = "timepwt48"))
  hc <- quietly(iv_regression(iv_fixture_fm, iv_fixture_data(), 4, 4, cov = "HC", weights = "timepwt48"))
  im <- quietly(iv_im_regression(iv_fixture_fm, iv_fixture_data(), 4, 4, clusters = 4, weights = "timepwt48"))
  expect_length(bch$AR_Set, 4)               # two rays
  expect_identical(hc$AR_Set, c(-Inf, Inf))  # whole line
  plot_segs <- function(p) purrr::list_rbind(lapply(p$layers, function(l) if (inherits(l$geom, "GeomSegment")) l$data))
  p <- iv_regression_plot(bch = bch, im = im)
  expect_s3_class(p, "ggplot")
  segs <- plot_segs(p)
  expect_true(all(is.finite(c(segs$x, segs$xend))))
  expect_equal(nrow(segs), nrow(ar_segments(bch$AR_Set, "b", 0)) + nrow(ar_segments(im$AR_Set, "b", 0)))
  expect_true(all(c("first", "last") %in% segs$arrow))
  expect_identical(levels(segs$term)[1], "fit_d_tradeusch_pw")
  expect_identical(plot_segs(iv_regression_plot(bch = hc))$arrow, "both")
  expect_no_error(iv_regression_plot(orig = quietly(fixest::feols(iv_fixture_fm, iv_fixture_data()))))
})

test_that("ar_segments() splits AR sets into line segments", {
  expect_equal(ar_segments(c(-1, 2), "b", 1)[, c("x", "xend")], data.frame(x = -1, xend = 2))
  expect_equal(ar_segments(c(-Inf, -1, 2, Inf), "b", 1)[, c("x", "xend")],
               data.frame(x = c(-Inf, 2), xend = c(-1, Inf)))
})

# With the treatment as its own instrument 2SLS is least squares, so SCPC inference on 2SLS must give
# the results of scpc_regression(). fixest rejects an exact copy, so the instrument is perturbed by
# about 1e-4 of its scale.
test_that("SCPC inference on 2SLS reduces to scpc_regression() when the treatment is its own instrument", {
  d <- iv_fixture_data()
  d$own <- d$d_tradeusch_pw * (1 + 1e-4 * cos(seq_len(nrow(d))))
  for (weights in list(NULL, "timepwt48")) {
    iv <- quietly(iv_scpc_regression(d_sh_empl_mfg ~ l_shind_manuf_cbp + l_sh_popedu_c | d_tradeusch_pw ~ own,
                                     d, 4, 4, avc = 0.05, weights = weights))
    ols <- quietly(scpc_regression(d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp + l_sh_popedu_c,
                                   d, 4, 4, avc = 0.05, weights = weights))
    ols_coefs <- ols$Coefs[match(sub("^fit_", "", iv$Coefs$term), ols$Coefs$term), ]
    expect_equal(iv$Coefs[, -1], ols_coefs[, -1], tolerance = 1e-3, ignore_attr = TRUE)
  }
})

test_that("iv_scpc_regression() returns the documented components and checks avc", {
  sc <- fit_iv("SCPC")
  expect_s3_class(sc, c("iv_scpc_regression", "iv_regression"), exact = TRUE)
  expect_named(sc, c("Coefs", "AR", "AR_Set", "F", "Inference", "avc", "N", "Splines", "PCs", "Fit"))
  expect_identical(sc$Coefs$term, c("(Intercept)", "fit_d_tradeusch_pw", "l_shind_manuf_cbp", "l_sh_popedu_c"))
  expect_true(is.na(sc$F))
  expect_output(print(sc), "SCPC inference with avc = 0.05")
  expect_output(print(sc), "Anderson-Rubin test")
  expect_error(iv_scpc_regression(iv_fixture_fm, iv_fixture_data(), 4, 4, avc = c(0.05, 0.1)), "single number")
  expect_error(iv_scpc_regression(iv_fixture_fm, iv_fixture_data(), 4, 4, avc = 2), "between")
})

test_that("iv_regression_table() and iv_regression_plot() take iv_scpc_regression() results", {
  skip_if_not_installed("modelsummary")
  sc <- fit_iv("SCPC")
  bch <- fit_iv("BCH")
  tab <- iv_regression_table(bch = bch, scpc = sc, output = "data.frame")
  expect_identical(names(tab)[4:5], c("BCH", "SCPC"))
  row <- function(term) unname(unlist(tab[tab$term == term, 4:5]))
  expect_identical(row("First stage F"), c(format(signif(bch$F, 3)), ""))
  expect_identical(row("AR 95% set"), c(bch$AR$CI, sc$AR$CI))
  expect_error(iv_regression_table(bch = sc), "iv_regression()")
  expect_error(iv_regression_table(scpc = bch), "iv_scpc_regression()")
  expect_error(regression_table(scpc = sc), "iv_regression_table")
  p <- iv_regression_plot(bch = bch, scpc = sc)
  expect_s3_class(p, "ggplot")
  expect_true("SCPC" %in% p$layers[[2]]$data$model)
})

test_that("iv_combined_table() picks the chosen 2SLS and AR rows and adds synth p values", {
  plbo <- readRDS(golden_path("iv_placebo"))
  syn <- readRDS(golden_path("iv_synth"))
  pick <- function(res, se, level) res[res$SE == se & res$Clusters == level, ]
  tab <- iv_combined_table(plbo, syn, select = c(BCH = 4, IM = 5, SCPC = 0.05), caption = "x")
  expect_s4_class(tab, "tinytable")
  expect_identical(names(tab@data),
    c("Adj", "Clusters/avc", "F", "Est p", "Plac p", "Plac 5%", "Synth p", "CI Width", "CI"))
  expect_identical(tab@data$Adj, rep(c("HC", "BCH", "IM", "SCPC"), 2))
  expect_identical(tab@data$`Clusters/avc`, rep(c(".", "4", "5", "0.05"), 2))
  expect_identical(tab@data$F[c(3, 4, 7, 8)], rep("", 4))   #no F for IM and SCPC
  expect_equal(tab@data$`Est p`[c(2, 6)], c(pick(plbo$Results, "BCH", "4")$`Est p`,
                                            pick(plbo$Results, "AR BCH", "4")$`Est p`))
  expect_equal(tab@data$`Synth p`[c(4, 8)], c(pick(syn$Results, "SCPC", "0.05")$`Synth p`,
                                              pick(syn$Results, "AR SCPC", "0.05")$`Synth p`))
  expect_identical(tab@data$CI[5:8], plbo$Results$CI[match(paste("AR", c("HC", "BCH", "IM", "SCPC"), c(".", "4", "5", "0.05")),
                                                           paste(plbo$Results$SE, plbo$Results$Clusters))])
  expect_true(any(grepl("Instrument: Structure", unlist(tab@notes))))

  # AR rows can use their own choices of inference, clusters and avc.
  ar <- iv_combined_table(plbo, syn, select = c(BCH = 4), ar_select = c(IM = 3))
  expect_identical(ar@data$Adj, c("HC", "BCH", "HC", "IM"))
  expect_identical(names(ar@data)[2], "Clusters")
  ar <- iv_combined_table(plbo, syn, select = c(BCH = 4, SCPC = 0.05), ar_select = c(BCH = 5, SCPC = 0.1))
  expect_identical(ar@data$`Clusters/avc`, c(".", "4", "0.05", ".", "5", "0.1"))
  expect_equal(ar@data$`Est p`[5:6], c(pick(plbo$Results, "AR BCH", "5")$`Est p`,
                                       pick(plbo$Results, "AR SCPC", "0.1")$`Est p`))
  expect_error(iv_combined_table(plbo, syn, select = c(BCH = 4), ar_select = c(IM = 9)), "Available values")

  expect_error(iv_combined_table(plbo, syn, select = c(BCH = 9)), "Available values: 3, 4, 5")
  expect_error(iv_combined_table(plbo, syn, select = c(4, 5)), "named vector")
  expect_error(iv_combined_table(plbo, syn, select = c(OLS = 4)), "named vector")
  expect_error(iv_combined_table(syn, syn, select = c(BCH = 4)), "not the output of iv_placebo")
  expect_error(iv_combined_table(plbo, readRDS(golden_path("synth")), select = c(BCH = 4)),
               "not the output of iv_synth")
  expect_error(iv_combined_table(plbo, readRDS(golden_path("iv_placebo_weighted")), select = c(BCH = 4)),
               "not the output of iv_synth")
})
