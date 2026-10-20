# Instrumental variable placebo and synthetic outcome tests: iv_placebo() and iv_synth().

iv_se_labels <- function(n_clus) {
  c("HC", rep("BCH", n_clus), rep("IM", n_clus), "AR HC", rep("AR BCH", n_clus))
}

test_that("iv_placebo() and iv_synth() results are unchanged", {
  expect_golden("iv_placebo")
  expect_golden("iv_synth")
  expect_golden("iv_placebo_weighted")
})

test_that("IV results have the documented components and columns", {
  n_clus <- fixture_args$max_clus - 2
  plbo <- readRDS(golden_path("iv_placebo"))
  expect_null(attr(plbo, "class"))
  expect_named(plbo, c("Results", "Spatial_Params"))
  expect_named(plbo$Results,
    c("SE", "Clusters", "F", "Est p", "Plac p", "Plac 5%", "CI Width", "CI"))
  expect_equal(plbo$Results$SE, iv_se_labels(n_clus))
  expect_equal(plbo$Results$Clusters,
    c(".", rep(as.character(3:fixture_args$max_clus), 2), ".", as.character(3:fixture_args$max_clus)))
  expect_named(plbo$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "R2_Inst", "Effective_Range_Inst",
      "Structure_Inst", "N", "Splines", "PCs"))

  syn <- readRDS(golden_path("iv_synth"))
  expect_null(attr(syn, "class"))
  expect_named(syn$Results, c("SE", "Clusters", "F", "Est p", "Synth p", "CI Width", "CI"))
  expect_equal(syn$Results$SE, iv_se_labels(n_clus))
  expect_named(syn$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "N", "Splines", "PCs"))

  # Estimates with the real data do not depend on the simulations.
  shared <- c("SE", "Clusters", "F", "Est p", "CI Width", "CI")
  expect_identical(plbo$Results[shared], syn$Results[shared])
})

# F must follow the inference and number of clusters of each row: HC and BCH from the first stage
# t statistic of the instrument. IM rows have no F. AR rows (HC and BCH only) share the F of their inference.
test_that("first stage F statistics follow the inference and number of clusters", {
  res <- readRDS(golden_path("iv_placebo"))$Results
  prep <- quietly(prepare_iv_data(iv_fixture_fm, iv_fixture_data(), 4, 4, NULL, 5))
  df <- prep$df
  eqs <- build_iv_formulas(prep$rhs, prep$pc)
  fs <- quietly(fixest::feols(eqs$eq_fs, df, weights = ~wts, vcov = "hetero"))
  fit <- quietly(fixest::feols(eqs$eq_est, df, weights = ~wts, vcov = "hetero"))
  hold_clus <- quietly(iv_clusters(df, TRUE, 5))
  expect_equal(res$F[res$SE == "HC"], round(fixest::fitstat(fit, "ivwald1", simplify = TRUE)$stat, 3))
  for (clust in hold_clus) {
    G <- as.character(nlevels(clust))
    bch_fit <- quietly(fixest::feols(eqs$eq_est, df, weights = ~wts, cluster = clust))
    expect_equal(res$F[res$SE == "BCH" & res$Clusters == G],
                 round(fixest::fitstat(bch_fit, "ivwald1", simplify = TRUE)$stat, 3))
  }
  expect_true(all(is.na(res$F[res$SE == "IM"])))
  expect_false("AR IM" %in% res$SE)
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

# The AR confidence set is found analytically; at its finite endpoints a direct regression of
# outcome - b*treatment on the instrument must give a p value of exactly 0.05. The full 2000-2007
# data have a strong instrument and bounded sets; the weak instrument of the fixture gives rays.
ar_endpoint_p <- function(data, splines, pc_num, clust_index) {
  prep <- quietly(prepare_iv_data(iv_fixture_fm, data, splines, pc_num, "timepwt48", 6))
  df <- prep$df
  eqs <- build_iv_formulas(prep$rhs, prep$pc)
  fits <- quietly(ar_fits(df, eqs))
  clust <- quietly(iv_clusters(df, TRUE, 6))[[clust_index]]
  p_at <- function(b, inference) {
    df$ar_diff <- df$dep_var - b * df$explan_var
    fit <- quietly(fixest::feols(eqs$eq_diff, df, weights = ~wts, vcov = "hetero"))
    if (inference == "HC") fixest::pvalue(fit)[["inst_var"]] else fixest::pvalue(fit, cluster = clust)[["inst_var"]]
  }
  sets <- list(
    HC = ar_set_fixest(fits),
    BCH = ar_set_fixest(fits, clust))
  lapply(names(sets), function(inference) {
    set <- sets[[inference]]
    list(set = set, p = vapply(set[is.finite(set)], p_at, numeric(1), inference = inference))
  })
}

test_that("AR confidence sets have p values of 0.05 at their endpoints", {
  strong <- spatInfer::autor_china[spatInfer::autor_china$t2 == 1, ]
  for (res in ar_endpoint_p(strong, 6, 8, 2)) {
    expect_length(res$set, 2)
    expect_equal(res$p, c(0.05, 0.05), tolerance = 1e-8)
  }
  weak <- ar_endpoint_p(iv_fixture_data(), 4, 4, 1)
  expect_true(any(vapply(weak, function(res) length(res$set) == 4, logical(1))))
  for (res in weak) expect_equal(res$p, rep(0.05, length(res$p)), tolerance = 1e-8)
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

test_that("placebo_table() and synth_table() split IV results into BCH, IM and AR tables", {
  plbo <- readRDS(golden_path("iv_placebo"))
  n_clus <- fixture_args$max_clus - 2
  expect_equal(iv_table_rows(plbo$Results, "BCH")$SE, c("HC", rep("BCH", n_clus)))
  expect_equal(iv_table_rows(plbo$Results, "IM")$SE, c("HC", rep("IM", n_clus)))
  ar <- iv_table_rows(plbo$Results, "AR")
  expect_equal(ar$SE, c("HC", rep("BCH", n_clus)))
  expect_identical(ar$CI, plbo$Results$CI[startsWith(plbo$Results$SE, "AR ")])
  for (type in c("BCH", "IM", "AR")) {
    tab <- placebo_table(plbo, type = type)
    expect_s4_class(tab, "tinytable")
    expect_equal(nrow(tab@data), 1 + n_clus)
    expect_true("F" %in% names(tab@data))
    expect_true(any(grepl("Noise instrument", unlist(tab@notes))))
    syn_tab <- synth_table(readRDS(golden_path("iv_synth")), type = type)
    expect_identical(nrow(syn_tab@data), nrow(tab@data))
    expect_false(any(grepl("Noise instrument", unlist(syn_tab@notes))))
  }
  expect_true(any(grepl("Anderson-Rubin", unlist(placebo_table(plbo, type = "AR")@notes))))
  expect_error(placebo_table(plbo, type = "SCPC"))
  # Other results are unaffected by type.
  expect_identical(placebo_table(readRDS(golden_path("placebo")), type = "AR")@data,
                   placebo_table(readRDS(golden_path("placebo")))@data)
})

# iv_regression() and iv_im_regression() use the same clusters and tests as iv_placebo(), so their
# results for the real data must agree with the iv_placebo() rows for the same inference and clusters.
fit_iv <- function(cov, clusters = 4, data = iv_fixture_data()) {
  if (cov == "IM") return(quietly(iv_im_regression(iv_fixture_fm, data, 4, 4, clusters = clusters)))
  quietly(iv_regression(iv_fixture_fm, data, 4, 4, clusters = if (cov == "HC") NULL else clusters, cov = cov))
}

test_that("iv_regression() and iv_im_regression() agree with iv_placebo()", {
  plbo <- readRDS(golden_path("iv_placebo"))$Results
  row <- function(se, clusters) plbo[plbo$SE == se & plbo$Clusters == clusters, ]
  for (cov in c("BCH", "IM", "HC")) {
    reg <- fit_iv(cov)
    k <- if (cov == "HC") "." else "4"
    wald <- row(cov, k)
    treat <- reg$Coefs[reg$Coefs$term == "fit_d_tradeusch_pw", ]
    expect_equal(round(treat$p.value, 3), wald$`Est p`)
    expect_identical(paste0("[", signif(treat$conf.low, 3), ", ", signif(treat$conf.high, 3), "]"), wald$CI)
    expect_equal(round(reg$F, 3), wald$F)   #both NA for IM
    if (cov == "IM") next   #IM is not valid for AR
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
  expect_named(im, c("Coefs", "F", "Inference", "Clusters", "N", "Splines", "PCs",
                     "Cluster_Coefs", "Cluster_N"))
  expect_identical(dim(im$Cluster_Coefs), c(4L, 5L))
  expect_equal(sum(im$Cluster_N), nrow(iv_fixture_data()))
  expect_output(print(im), "IM inference with 4 clusters")
  expect_output(print(im), "Coefficients in each cluster")
  expect_output(print(bch), "Anderson-Rubin test")
  expect_false(any(grepl("Anderson-Rubin", capture.output(print(im)))))
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
  expect_identical(unname(row("AR p")), c("", format(signif(bch$AR$p.value, 3)), ""))
  expect_identical(unname(row("AR 95% set")), c("", bch$AR$CI, ""))
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
  expect_setequal(segs$arrow, c("first", "last"))   # no AR set for IM
  expect_identical(levels(segs$term)[1], "fit_d_tradeusch_pw")
  expect_identical(plot_segs(iv_regression_plot(bch = hc))$arrow, "both")
  expect_no_error(iv_regression_plot(orig = quietly(fixest::feols(iv_fixture_fm, iv_fixture_data()))))
})

test_that("ar_segments() splits AR sets into line segments", {
  expect_equal(ar_segments(c(-1, 2), "b", 1)[, c("x", "xend")], data.frame(x = -1, xend = 2))
  expect_equal(ar_segments(c(-Inf, -1, 2, Inf), "b", 1)[, c("x", "xend")],
               data.frame(x = c(-Inf, 2), xend = c(-1, Inf)))
})
