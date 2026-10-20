test_that("missing or incomplete coordinates are rejected", {
  d <- fixture_data()
  no_x <- d[, setdiff(names(d), "X")]
  na_y <- d
  na_y$Y[1] <- NA

  for (fun in list(placebo, placebo_im, synth, synth_im)) {
    expect_error(fun(fixture_fm, no_x, 4, 4), "named X and Y")
    expect_error(fun(fixture_fm, na_y, 4, 4), "missing values")
  }
  expect_error(basis_regression(fixture_fm, no_x, 4, 4, 4), "named X and Y")
  expect_error(basis_regression(fixture_fm, na_y, 4, 4, 4), "missing values")
})

test_that("max_clus below 3 is rejected", {
  for (fun in list(placebo, placebo_im, synth, synth_im)) {
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, max_clus = 2), "greater than 2")
  }
})

test_that("weights must name a numeric variable without missing values", {
  d <- fixture_data(weighted = TRUE)
  d$wt_text <- as.character(d$wt)
  d$wt_na <- replace(d$wt, 1, NA)
  for (fun in list(placebo, placebo_im, synth, synth_im)) {
    expect_error(fun(fixture_fm, d, 4, 4, weights = "pop"), "no variable called pop")
    expect_error(fun(fixture_fm, d, 4, 4, weights = TRUE), "name of the weighting variable")
    expect_error(fun(fixture_fm, d, 4, 4, weights = "wt_text"), "must be numeric")
    expect_error(fun(fixture_fm, d, 4, 4, weights = "wt_na"), "missing values")
  }
  for (w in list("pop", TRUE, c("wt", "wt_na")))
    expect_error(basis_regression(fixture_fm, d, 4, 4, 4, weights = w))
})

test_that("optimal_basis() rejects max_splines outside 3 to 12", {
  expect_error(optimal_basis(fixture_fm, fixture_data(), max_splines = 13), "12")
  expect_error(optimal_basis(fixture_fm, fixture_data(), max_splines = 2), "minimum")
})

test_that("exact_cholesky = FALSE without BRISC gives an informative error", {
  local_mocked_bindings(has_brisc = function() FALSE)
  expect_error(run_sim(placebo, Parallel = FALSE, exact_cholesky = FALSE),
    "needs the BRISC package")
})
