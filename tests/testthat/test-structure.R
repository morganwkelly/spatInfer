# Returned object structures that downstream code and the table functions rely on.

test_that("optimal_basis(max_splines = 3) uses only the 3x3 tensor", {
  expect_golden("optimal_basis_max3")
  ob3 <- readRDS(golden_path("optimal_basis_max3"))
  expect_identical(unique(ob3$bic$name), "BIC_3")
  expect_identical(unique(ob3$r2$name), "R2_3")
  expect_match(ob3$subtitle, "3x3 Linear Tensor")
})

test_that("placebo-type results have the documented components and columns", {
  plbo <- readRDS(golden_path("placebo"))
  expect_type(plbo, "list")
  expect_null(attr(plbo, "class"))
  expect_named(plbo, c("Results", "Spatial_Params"))
  expect_named(plbo$Results,
    c("SE", "Clusters", "Est p", "Plac p", "Plac 5%", "CI Width", "CI"))
  expect_named(plbo$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "Splines", "PCs"))
  expect_equal(plbo$Results$SE, c("HC", rep("BCH", fixture_args$max_clus - 2)))

  pim <- readRDS(golden_path("placebo_im"))
  expect_named(pim, c("Results", "Spatial_Params", "Coefs"))
  expect_length(pim$Coefs, fixture_args$max_clus - 1)
  expect_equal(pim$Results$SE, c("HC", rep("IM", fixture_args$max_clus - 2)))
})

test_that("synth-type results have the documented components and columns", {
  syn <- readRDS(golden_path("synth"))
  expect_null(attr(syn, "class"))
  expect_named(syn, c("Results", "Spatial_Params"))
  expect_named(syn$Results,
    c("SE", "Clusters", "Est p", "Synth p", "CI Width", "CI"))
  expect_named(syn$Spatial_Params,
    c("Moran", "R2", "Effective_Range", "Structure", "N", "Splines", "PCs"))

  # synth_im() is the only function that sets a class.
  syn_im <- readRDS(golden_path("synth_im"))
  expect_s3_class(syn_im, "synth_im")
  expect_named(syn_im$Results, names(syn$Results))
})

test_that("basis_regression() returns a fixest object with saved data", {
  m <- quietly(basis_regression(fixture_fm, fixture_data(),
    splines = 4, pc_num = 4, clusters = 4))
  expect_s3_class(m, "fixest")
  expect_identical(rownames(m$coeftable)[2], "single_mothers")
  expect_true(all(paste0("PC", 1:4) %in% rownames(m$coeftable)))
})

test_that("basis_regression() applies weights = TRUE", {
  d <- fixture_data(weighted = TRUE)
  unweighted <- quietly(basis_regression(fixture_fm, d, 4, 4, 4))
  weighted <- quietly(basis_regression(fixture_fm, d, 4, 4, 4, weights = TRUE))
  expect_false(isTRUE(all.equal(coef(weighted), coef(unweighted))))
  expect_equal(unname(weights(weighted)), d$weights)

  # Same coefficients as a weighted lm() fit on the regression's own data.
  check <- lm(formula(weighted), data = weighted$data, weights = weights)
  expect_equal(coef(weighted), coef(check))
})

test_that("placebo_table() and synth_table() return tinytables", {
  plac_tab <- placebo_table(readRDS(golden_path("placebo")), caption = "x")
  expect_s4_class(plac_tab, "tinytable")
  expect_identical(names(plac_tab@data),
    c("Adj", "Clusters", "Est p", "Plac p", "Plac 5%", "CI Width", "CI"))

  syn_tab <- synth_table(readRDS(golden_path("synth_im")))
  expect_s4_class(syn_tab, "tinytable")
  expect_identical(names(syn_tab@data),
    c("Adj", "Clusters", "Est p", "Synth p", "CI Width", "CI"))
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

test_that("jitter_coords controls jittering of identical coordinates in the Moran test", {
  dup <- fixture_data()
  dup$X[2:6] <- dup$X[1]
  dup$Y[2:6] <- dup$Y[1]

  jittered <- run_sim(placebo, Parallel = FALSE, data = dup, jitter_coords = TRUE)
  unjittered <- run_sim(placebo, Parallel = FALSE, data = dup, jitter_coords = FALSE)
  expect_false(identical(jittered$Spatial_Params$Moran, unjittered$Spatial_Params$Moran))
  # Only the Moran test uses the jittered coordinates.
  expect_identical(jittered$Results, unjittered$Results)
  expect_identical(jittered$Spatial_Params[-1], unjittered$Spatial_Params[-1])

  prep <- prepare_spatial_data(fixture_fm, dup, 4, 4, FALSE, 5)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  warnings <- capture_warnings(moran(eqs$eq_est, prep$df, jitter_coords = FALSE))
  expect_match(warnings, "identical points", all = FALSE)
  expect_no_warning(moran(eqs$eq_est, prep$df, jitter_coords = TRUE))

  # Without identical coordinates the setting makes no difference.
  expect_identical(run_sim(placebo, Parallel = FALSE, jitter_coords = FALSE),
    readRDS(golden_path("placebo")))
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
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, FALSE, 5)
  df <- prep$df
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  hold <- generate_clusters(df, TRUE, 5)
  # The treatment does not vary within cluster 1 of the four cluster partition.
  df$explan_var[hold[, 3] == 1] <- 0
  im_out <- as.data.frame(matrix(0.5, nrow = 5, ncol = 4))
  hc_out <- data.frame(hc_p = rep(0.5, 5))

  warnings <- capture_warnings(quietly_notes(summary_im(df, eqs$eq_est, hold, 5, im_out, hc_out)))
  expect_match(warnings, "IM with 4 clusters: .* in 1 cluster", all = FALSE)

  coefs <- quietly(im_cluster_coefs(eqs$eq_est, df, hold[, 3]))
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
