# The package must not change global state that belongs to the user.

test_that("simulation functions leave the user's random number state unchanged", {
  for (fun in list(placebo, placebo_im, synth, synth_im)) {
    set.seed(42)
    expected <- runif(3)
    set.seed(42)
    run_sim(fun, Parallel = FALSE)
    expect_identical(runif(3), expected)
  }
})

test_that("preserving the random number state does not change results", {
  set.seed(1)
  first <- run_sim(placebo, Parallel = FALSE)
  set.seed(2)
  second <- run_sim(placebo, Parallel = FALSE)
  expect_identical(first, second)
  expect_identical(first, readRDS(golden_path("placebo")))
})

test_that("plot_basis() restores graphics parameters", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  graphics::par(mai = c(1, 1, 1, 1))
  quietly(plot_basis(fixture_fm, fixture_data(), splines = 4))
  expect_equal(graphics::par("mai"), c(1, 1, 1, 1))
})
