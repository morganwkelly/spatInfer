# Clusters shared by placebos and final regressions, treatment found by name, SCPC failures
# for each avc, the range search warning and the synth trend.

fixture_coords <- function() as.matrix(fixture_data()[, c("X", "Y")])

test_that("final regressions use the placebo's clusters, including Clara", {
  hold_clus <- generate_clusters(fixture_data(), k_medoids = FALSE, max_clus = 4)
  clara <- hold_clus$clust_4
  expect_false(isTRUE(all.equal(as.integer(clara),
                                as.integer(make_clusters(fixture_coords(), 4, TRUE)))))

  bch <- quietly(basis_regression(fixture_fm, fixture_data(), 4, 4, clusters = 4, k_medoids = FALSE))
  expect_equal(fixest::se(bch), fixest::se(summary(bch, cluster = clara)))

  im <- quietly(im_regression(fixture_fm, fixture_data(), 4, 4, clusters = 4, k_medoids = FALSE))
  expect_equal(unname(im$Cluster_N), unname(c(table(clara))))

  # Clara uses a fixed seed and leaves the user's random numbers unchanged.
  set.seed(42)
  before <- .Random.seed
  expect_identical(unname(make_clusters(fixture_coords(), 4, FALSE)), clara)
  expect_identical(.Random.seed, before)
})

test_that("BCH placebo p values are those of a separate fit for each number of clusters", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  hold_clus <- generate_clusters(prep$df, TRUE, 4)
  set.seed(1)
  Sim <- matrix(rnorm(2 * nrow(prep$df)), ncol = 2)
  out <- bch_p_values(Sim, eqs$eq_sim, prep$df, hold_clus, 2, 4, Parallel = FALSE)
  expect_named(out, c("clus_2", "clus_3", "clus_4"))
  for (j in 1:2) for (k in 1:3) {
    df1 <- cbind.data.frame(sim = Sim[, j], prep$df, clust = hold_clus[, k])
    fit <- fixest::feols(eqs$eq_sim, data = df1, weights = ~wts, vcov = fixest::vcov_cluster(cluster = ~clust))
    expect_equal(out[j, k], unname(fixest::pvalue(fit)["sim"]))
  }
})

test_that("a treatment dropped as collinear is an error, not the next coefficient", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  Sim <- matrix(0, nrow = nrow(prep$df), ncol = 1)
  hold_clus <- generate_clusters(prep$df, TRUE, 3)
  expect_error(quietly(hc_sim(1, Sim, eqs$eq_sim, prep$df)), "dropped from the regression")
  expect_error(quietly(bch_sim(1, Sim, eqs$eq_sim, prep$df, hold_clus)), "dropped from the regression")
})

test_that("SCPC failures are left out separately for each avc", {
  prep <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)
  eqs <- build_formulas(prep$rhs, prep$pc, "explan_var")
  avc <- c(0.02, 0.1)
  scpc_out <- data.frame(avc_0.02 = c(0.01, NA, 0.5, 0.5), avc_0.1 = c(0.01, 0.01, 0.5, 0.9))
  hc_out <- data.frame(hc_p = rep(0.5, 4))
  warnings <- capture_warnings(summ <- quietly_notes(summary_scpc(prep$df, eqs$eq_est, avc, scpc_out, hc_out)))
  expect_match(warnings, "1 of 4 simulations for avc = 0.02", all = FALSE)
  expect_false(any(grepl("avc = 0.1", warnings)))
  expect_equal(summ$Results$sim_05[-1], round(c(1 / 3, 2 / 4), 3))
})

test_that("Noise_Sim() warns when the range is at the upper limit of its search", {
  d <- fixture_data()
  expect_warning(Noise_Sim(d, lm(X ~ 1, data = d), 2, TRUE, FALSE), "upper limit of the search")
  expect_no_warning(Noise_Sim(d, lm(single_mothers ~ X + Y, data = d), 2, TRUE, FALSE))
})

test_that("synth trend adds the product of X and Y only when asked", {
  d <- prepare_spatial_data(fixture_fm, fixture_data(), 4, 4, NULL)$df
  expect_false("X:Y" %in% names(coef(synth_trend(d, FALSE))))
  expect_true("X:Y" %in% names(coef(synth_trend(d, TRUE))))
  expect_gte(summary(synth_trend(d, TRUE))$r.squared, summary(synth_trend(d, FALSE))$r.squared)

  syn <- run_sim(synth, xy_interaction = TRUE)
  expect_named(syn, c("Results", "Spatial_Params"))
  expect_equal(syn$Spatial_Params$R2, round(summary(synth_trend(d, TRUE))$r.squared, 3))
})
