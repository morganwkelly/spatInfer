# Parallel execution must reproduce serial results exactly: the simulated noise
# is generated serially before any parallel step.

skip_parallel <- function() {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(isTRUE(parallel::detectCores() < 4), "Fewer than 4 cores available")
}

test_that("placebo() gives identical results in parallel", {
  skip_parallel()
  expect_identical(run_sim(placebo, Parallel = TRUE), readRDS(golden_path("placebo")))
})

test_that("placebo_im() gives identical results in parallel", {
  skip_parallel()
  expect_identical(run_sim(placebo_im, Parallel = TRUE), readRDS(golden_path("placebo_im")))
})
