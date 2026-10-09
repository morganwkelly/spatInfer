# Characterization tests: outputs must match those recorded from commit ed80437.

test_that("optimal_basis() BIC/R2 curves and chosen basis are unchanged", {
  expect_golden("optimal_basis")
})

test_that("basis_regression() estimates are unchanged", {
  expect_golden("basis_regression_bch")
  expect_golden("basis_regression_hc")
  expect_golden("basis_regression_weighted")
})

test_that("placebo() and placebo_im() results are unchanged", {
  expect_golden("placebo")
  expect_golden("placebo_im")
})

test_that("synth() and synth_im() results are unchanged", {
  expect_golden("synth")
  expect_golden("synth_im")
})

test_that("weighted simulation results are unchanged", {
  expect_golden("placebo_weighted")
  expect_golden("synth_im_weighted")
})

test_that("clara clustering (k_medoids = FALSE) results are unchanged", {
  expect_golden("placebo_clara")
})

test_that("BRISC simulation (exact_cholesky = FALSE) results are unchanged", {
  skip_if_not_installed("BRISC")
  expect_golden("placebo_brisc")
})
