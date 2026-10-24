# Returned object structures that downstream code and the table functions rely on.

test_that("optimal_basis(max_splines = 3) uses only the 3x3 tensor", {
  expect_golden("optimal_basis_max3")
  ob3 <- readRDS(golden_path("optimal_basis_max3"))
  expect_identical(unique(ob3$bic$name), "BIC_3")
  expect_identical(unique(ob3$r2$name), "R2_3")
  expect_match(ob3$subtitle, "3x3 Linear Tensor")
})

test_that("optimal_basis() has a distinct colour for every tensor up to 12x12", {
  for (n in 1:10) {
    colours <- basis_colours(n)
    expect_length(colours, n)
    expect_false(anyNA(colours))
    expect_identical(anyDuplicated(colours), 0L)
  }
  # Plots with up to eight tensors keep their colours.
  expect_identical(basis_colours(8),
    c("#999999", "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00"))
})

test_that("placebo-type results have the documented components and columns", {
  plbo <- readRDS(golden_path("placebo"))
  expect_type(plbo, "list")
  expect_null(attr(plbo, "class"))
  expect_named(plbo, c("Results", "Spatial_Params"))
  expect_named(plbo$Results,
    c("SE", "Clusters", "Est p", "Plac p", "Plac 5%", "CI Width", "CI"))
  expect_named(plbo$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "N", "Splines", "PCs"))
  n_clus <- fixture_args$max_clus - 2
  expect_equal(plbo$Results$SE,
    c("HC", rep("BCH", n_clus), rep("IM", n_clus), rep("SCPC", length(fixture_avc))))
  expect_equal(plbo$Results$Clusters,
    c(".", rep(as.character(3:fixture_args$max_clus), 2), as.character(fixture_avc)))
})

test_that("synth-type results have the documented components and columns", {
  syn <- readRDS(golden_path("synth"))
  expect_null(attr(syn, "class"))
  expect_named(syn, c("Results", "Spatial_Params"))
  expect_named(syn$Results,
    c("SE", "Clusters", "Est p", "Synth p", "CI Width", "CI"))
  expect_named(syn$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "N", "Splines", "PCs"))
})

test_that("each simulation is estimated once for HC, BCH, IM and SCPC p values", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  hold_clus <- sim_clusters(prep$df, TRUE, 4)
  set.seed(1)
  Sim <- matrix(rnorm(nrow(prep$df)), ncol = 1)
  out <- quietly(ols_sim(1, Sim, eqs$eq_sim, prep$df, hold_clus, 0.05))
  expect_named(out, c("hc", "bch_3", "im_3", "bch_4", "im_4", "scpc_0.05"))
  df1 <- cbind.data.frame(sim = Sim[, 1], prep$df)
  hc <- quietly(fixest::feols(eqs$eq_sim, data = df1, weights = ~wts, vcov = "hetero"))
  expect_equal(out$hc, unname(fixest::pvalue(hc)["sim"]))
  im <- quietly(im_t_test(im_cluster_coefs(eqs$eq_sim, df1, hold_clus[[2]])))
  expect_equal(out$im_4, im$p.value)
  expect_equal(out$scpc_0.05, unname(scpc_stats(hc, df1, 0.05)[, "p"]))
})

test_that("bch_regression() returns a fixest object with saved data", {
  m <- quietly(bch_regression(fixture_fm, fixture_data(),
    splines = 4, pc_num = 4, clusters = 4))
  expect_s3_class(m, "fixest")
  expect_identical(rownames(m$coeftable)[2], "single_mothers")
  expect_true(all(paste0("PC", 1:4) %in% rownames(m$coeftable)))
})

test_that("bch_regression() applies the weights it is given", {
  d <- fixture_data(weighted = TRUE)
  unweighted <- quietly(bch_regression(fixture_fm, d, 4, 4, 4))
  weighted <- quietly(bch_regression(fixture_fm, d, 4, 4, 4, weights = "wt"))
  expect_false(isTRUE(all.equal(coef(weighted), coef(unweighted))))
  expect_equal(unname(weights(weighted)), d$wt)

  # Same coefficients as a weighted lm() fit on the regression's own data.
  check <- lm(formula(weighted), data = weighted$data, weights = wt)
  expect_equal(coef(weighted), coef(check))
})

test_that("placebo_table() and synth_table() return tinytables", {
  plac_tab <- placebo_table(readRDS(golden_path("placebo")), caption = "x")
  expect_s4_class(plac_tab, "tinytable")
  expect_identical(names(plac_tab@data),
    c("Adj", "Clusters/avc", "Est p", "Plac p", "Plac 5%", "CI Width", "CI"))
  expect_identical(nrow(plac_tab@data), nrow(readRDS(golden_path("placebo"))$Results))

  syn_tab <- synth_table(readRDS(golden_path("placebo_clara")) |> (\(x) {
    names(x$Results)[4] <- "Synth p"
    x$Results$`Plac 5%` <- NULL
    x
  })())
  expect_s4_class(syn_tab, "tinytable")
  expect_identical(names(syn_tab@data),
    c("Adj", "Clusters", "Est p", "Synth p", "CI Width", "CI"))
})

test_that("combined_table() picks the chosen rows and adds their synth p values", {
  g <- function(f) readRDS(golden_path(f))
  plbo <- g("placebo")
  syn <- g("synth")
  tab <- combined_table(plbo, syn, select = c(BCH = 4, IM = 5, SCPC = 0.04), caption = "x")
  expect_s4_class(tab, "tinytable")
  expect_identical(names(tab@data),
    c("Adj", "Clusters/avc", "Est p", "Plac p", "Plac 5%", "Synth p", "CI Width", "CI"))
  expect_identical(tab@data$Adj, c("HC", "BCH", "IM", "SCPC"))
  expect_identical(tab@data$`Clusters/avc`, c(".", "4", "5", "0.04"))
  rows <- c(1, which(plbo$Results$SE == "BCH" & plbo$Results$Clusters == "4"),
            which(plbo$Results$SE == "IM" & plbo$Results$Clusters == "5"),
            which(plbo$Results$SE == "SCPC" & plbo$Results$Clusters == "0.04"))
  expect_equal(tab@data$`Est p`, round(plbo$Results$`Est p`[rows], 3))
  expect_equal(tab@data$`Synth p`, round(syn$Results$`Synth p`[rows], 3))

  # Rows follow the order of select.
  one <- combined_table(plbo, syn, c(IM = 3, BCH = 4))
  expect_identical(names(one@data)[2], "Clusters")
  expect_identical(one@data$Adj, c("HC", "IM", "BCH"))
  expect_identical(names(combined_table(plbo, syn, c(SCPC = 0.02))@data)[2], "avc")

  expect_error(combined_table(plbo, syn, c(BCH = 9)), "Available values: 3, 4, 5")
  expect_error(combined_table(plbo, syn, 4), "named vector")
  expect_error(combined_table(plbo, syn, c(BCH = 4, BCH = 5)), "named vector")
  expect_error(combined_table(syn, syn, c(BCH = 4)), "is not the output of placebo")
  expect_error(combined_table(plbo, plbo, c(BCH = 4)), "is not the output of synth")
  other_basis <- syn
  other_basis$Spatial_Params$PCs <- 5
  expect_error(combined_table(plbo, other_basis, c(BCH = 4)), "same regression")
  expect_error(combined_table(g("iv_placebo"), g("iv_synth"), c(BCH = 4)), "use iv_combined_table")
})

# Results objects with the HC row and n_clus cluster rows, as returned by placebo() and synth().
fake_placebo <- function(n_clus) {
  n <- n_clus + 1
  list(
    Results = data.frame(SE = c("HC", rep("BCH", n_clus)), Clusters = c(".", 2 + seq_len(n_clus)),
      `Est p` = 0.01, `Plac p` = 0.2, `Plac 5%` = 0.05, `CI Width` = 1, CI = "[-1, 0]",
      check.names = FALSE),
    Spatial_Params = data.frame(Moran = 1, R2 = 0.5, Effective_Range = 0.1, Structure = 0.9,
      N = 99, Splines = 4, PCs = 4)
  )
}

test_that("tables show each SE label once and never add rows", {
  for (n_clus in 1:6) {
    plbo <- fake_placebo(n_clus)
    n <- n_clus + 1
    expected_labels <- c("HC", "BCH", rep("", n - 2))
    plac_tab <- placebo_table(plbo)
    expect_identical(plac_tab@data$Adj, expected_labels)
    expect_false(anyNA(plac_tab@data))

    syn <- plbo
    names(syn$Results)[4] <- "Synth p"
    syn$Results$`Plac 5%` <- NULL
    syn_tab <- synth_table(syn)
    expect_identical(syn_tab@data$Adj, expected_labels)
    expect_false(anyNA(syn_tab@data))
  }
})

test_that("the Moran test handles identical coordinates without jittering", {
  dup <- fixture_data()
  dup$X[2:6] <- dup$X[1]
  dup$Y[2:6] <- dup$Y[1]

  prep <- prepare_spatial_data(fixture_fm, dup, 4, 4, NULL, 5)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  expect_no_warning(m1 <- moran(eqs$eq_est, prep$df))
  expect_true(is.finite(m1))
  # No random noise is added, so the statistic is reproducible.
  expect_identical(moran(eqs$eq_est, prep$df), m1)
})

test_that("plot_basis() passes theta and phi to vis.gam()", {
  local_mocked_bindings(vis.gam = function(x, ...) list(...), .package = "mgcv")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  args <- quietly(plot_basis(fixture_fm, fixture_data(), splines = 4, theta = 120, phi = 15))
  expect_identical(args$theta, 120)
  expect_identical(args$phi, 15)
  defaults <- quietly(plot_basis(fixture_fm, fixture_data(), splines = 4))
  expect_identical(c(defaults$theta, defaults$phi), c(30, 50))
})

test_that("IM leaves out clusters without a treatment estimate, with a warning", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL, 5)
  df <- prep$df
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  hold <- sim_clusters(df, TRUE, 5)
  # The treatment does not vary within cluster 1 of the four cluster partition.
  df$explan_var[hold[[2]] == 1] <- 0
  sim_out <- as.data.frame(as.list(c(hc = 0.5, setNames(rep(0.5, 6), paste0(c("bch_", "im_"), rep(3:5, each = 2))))))
  sim_out <- sim_out[rep(1, 5), ]

  warnings <- capture_warnings(quietly_notes(ols_summary(df, eqs$eq_est, hold, sim_out)))
  expect_match(warnings, "IM with 4 clusters: .* in 1 cluster", all = FALSE)

  coefs <- quietly(im_cluster_coefs(eqs$eq_est, df, hold[[2]]))
  expect_equal(nrow(coefs), 3)
  expect_identical(im_t_test(coefs)$p.value, t.test(coefs$estimate)$p.value)
  expect_identical(im_t_test(tibble::tibble(estimate = c(1, 2, NA, 4)))$p.value,
    t.test(c(1, 2, 4))$p.value)
})

test_that("plot_basis() draws without error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(quietly(plot_basis(fixture_fm, fixture_data(), splines = 4)))
})

test_that("confidence intervals are in the units of the treatment coefficient", {
  fmt <- function(ci) paste0("[", signif(ci[1], 3), ", ", signif(ci[2], 3), "]")
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eq <- build_formulas(prep$rhs, prep$pc, "explan_var")$eq_est
  hc <- quietly(fixest::feols(eq, data = prep$df, weights = ~wts, vcov = "hetero"))
  hc_ci <- confint(hc)[2, ]

  plbo <- readRDS(golden_path("placebo"))
  expect_identical(plbo$Results$CI[1], fmt(unlist(hc_ci)))
  expect_equal(plbo$Results$`CI Width`[1], signif(hc_ci[[2]] - hc_ci[[1]], 3))

})

test_that("small coefficients keep significant digits in confidence intervals", {
  row <- ci_summary_row("HC", 0, 0.01, c(0.000123, 0.000456), c(0.5, 0.5))
  expect_identical(row$CI, "[0.000123, 0.000456]")
  expect_equal(finish_summary(list(row))$width_ci, 0.000333)
})
