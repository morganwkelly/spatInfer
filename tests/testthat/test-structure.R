# Returned object structures that downstream code and the table functions rely on.

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

test_that("plot_basis() draws without error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(quietly(plot_basis(fixture_fm, fixture_data(), splines = 4)))
})
