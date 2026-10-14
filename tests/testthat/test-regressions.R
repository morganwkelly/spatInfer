# im_regression(), scpc_regression() and regression_table().

fit_im <- function(...) quietly(im_regression(fixture_fm, fixture_data(), 4, 4, clusters = 4, ...))
fit_scpc <- function(...) quietly(scpc_regression(fixture_fm, fixture_data(), 4, 4, avc = 0.04, ...))
fixture_terms <- c("(Intercept)", "single_mothers", "short_commute", "gini",
                   "dropout_rate", "social_cap", "dropout_na")

test_that("im_regression() tests every coefficient on its cluster estimates", {
  im <- fit_im()
  expect_s3_class(im, "im_regression")
  expect_identical(im$Coefs$term, fixture_terms)
  expect_named(im$Cluster_Coefs, c("cluster", fixture_terms))
  expect_identical(im$Cluster_Coefs$cluster, as.character(1:4))
  expect_equal(sum(im$Cluster_N), nrow(fixture_data()))

  # The treatment matches the IM t-test of placebo_im().
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, FALSE)
  eq <- build_formulas(prep$rhs, prep$pc)$eq_est
  clust <- factor(cluster::pam(as.matrix(fixture_data()[, c("X", "Y")]), k = 4)$clustering)
  tt <- im_t_test(quietly(im_cluster_coefs(eq, prep$df, clust)))
  expect_equal(im$Coefs$p.value[2], tt$p.value)
  expect_equal(im$Coefs$conf.low[2], tt$conf.int[1])

  # Pseudo standard errors are a quarter of the CI width; pseudo t uses them.
  expect_equal(im$Coefs$std.error, (im$Coefs$conf.high - im$Coefs$conf.low) / 4)
  expect_equal(im$Coefs$statistic, im$Coefs$estimate / im$Coefs$std.error)

  # Each IM estimate is the mean of the cluster estimates.
  for (k in seq_along(fixture_terms)) {
    est <- im$Cluster_Coefs[[fixture_terms[k]]]
    expect_equal(im$Coefs$estimate[k], mean(est[is.finite(est)]))
  }
})

test_that("im_regression() checks its inputs", {
  expect_error(im_regression(fixture_fm, fixture_data(), 4, 4, clusters = 1), "at least 2")
  expect_error(im_regression(fixture_fm, fixture_data(), 4, 4, clusters = 4, weights = TRUE),
               "no variable called weights")
})

test_that("im_table() makes coefficient and cluster tables", {
  im <- fit_im()
  coef_tab <- im_table(im)
  expect_identical(names(coef_tab@data),
    c("Variable", "Estimate", "Pseudo SE", "Pseudo t", "p", "CI", "Clusters"))
  clus_tab <- im_table(im, type = "clusters")
  expect_identical(names(clus_tab@data), c("Cluster", fixture_terms, "Observations"))
  expect_identical(clus_tab@data$Cluster, paste("Cluster", 1:4))
  expect_equal(sum(clus_tab@data$Observations), nrow(fixture_data()))
  expect_error(im_table(list()), "im_regression")
})

test_that("scpc_regression() matches scpcR::scpc() for every coefficient", {
  sc <- fit_scpc()
  expect_identical(sc$Coefs$term, fixture_terms)
  expect_error(scpc_regression(fixture_fm, fixture_data(), 4, 4, avc = c(0.02, 0.04)), "single number")
  skip_if_not_installed("scpcR")
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, FALSE)
  eq <- build_formulas(prep$rhs, prep$pc, "explan_var")$eq_est
  fit <- quietly(fixest::feols(eq, data = prep$df, weights = ~wts))
  direct <- scpcR::scpc(fit, data = prep$df, lon = "X", lat = "Y", avc = 0.04)$scpcstats
  expect_equal(unname(as.matrix(sc$Coefs[, -1])), unname(direct[1:7, ]))
})

test_that("regression_table() combines BCH, IM and SCPC with notes", {
  skip_if_not_installed("modelsummary")
  bch <- quietly(basis_regression(fixture_fm, fixture_data(), 4, 4, clusters = 4))
  tab <- regression_table(bch = bch, im = fit_im(), scpc = fit_scpc(), output = "data.frame")
  expect_identical(names(tab)[4:6], c("BCH", "IM", "SCPC"))
  expect_false(any(grepl("^PC|Intercept", tab$term)))

  tt <- regression_table(bch = bch, im = fit_im(), scpc = fit_scpc())
  expect_match(unlist(tt@notes), "BCH: 4 clusters. IM: 4 clusters. SCPC: avc = 0.04.", fixed = TRUE, all = FALSE)

  expect_identical(names(regression_table(im = fit_im(), output = "data.frame"))[4], "IM")
  orig <- quietly(fixest::feols(fixture_fm, fixture_data(), vcov = "hetero"))
  tab <- regression_table(orig = orig, bch = bch, im = fit_im(), output = "data.frame")
  expect_identical(names(tab)[4:6], c("Original", "BCH", "IM"))
  expect_identical(names(regression_table(orig = orig, output = "data.frame"))[4], "Original")
  expect_error(regression_table(orig = fit_im()), "feols")

  expect_error(regression_table(), "at least one")
  expect_error(regression_table(im = bch), "im_regression")
})

test_that("regression_table() shows confidence intervals, not standard errors", {
  skip_if_not_installed("modelsummary")
  tab <- regression_table(im = fit_im(), scpc = fit_scpc(), output = "data.frame")
  expect_true(all(tab$statistic[tab$part == "estimates"] %in% c("estimate", "conf.int")))
  expect_error(regression_table(im = fit_im(), statistic = "std.error", output = "data.frame"))
})

test_that("regression_plot() plots BCH, IM and SCPC with notes", {
  skip_if_not_installed("modelsummary")
  bch <- quietly(basis_regression(fixture_fm, fixture_data(), 4, 4, clusters = 4))
  p <- regression_plot(bch = bch, im = fit_im(), scpc = fit_scpc())
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$caption, "BCH: 4 clusters. IM: 4 clusters. SCPC: avc = 0.04.")
  dat <- regression_plot(bch = bch, im = fit_im(), scpc = fit_scpc(), draw = FALSE)
  expect_setequal(as.character(dat$model), c("BCH", "IM", "SCPC"))
  expect_false(any(grepl("^PC|Intercept", dat$term)))
  expect_error(regression_plot(), "at least one")
})
