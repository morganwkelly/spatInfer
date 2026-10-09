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

test_that("weights = TRUE requires a weights column", {
  for (fun in list(placebo, placebo_im, synth, synth_im)) {
    expect_error(fun(fixture_fm, fixture_data(), 4, 4, weights = TRUE), "no variable called weights")
  }
})

test_that("optimal_basis() rejects max_splines above 12", {
  expect_error(optimal_basis(fixture_fm, fixture_data(), max_splines = 13), "12")
})
