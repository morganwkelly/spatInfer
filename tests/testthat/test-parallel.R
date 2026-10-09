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

test_that("parallel runs restore the fixest thread setting", {
  skip_parallel()
  before <- fixest::getFixest_nthreads()
  run_sim(placebo, Parallel = TRUE)
  expect_identical(fixest::getFixest_nthreads(), before)
})

test_that("errors in parallel workers are raised", {
  skip_parallel()
  expect_error(run_sims(4, function(j) if (j == 3) stop("worker failed") else j, TRUE, FALSE),
    "worker failed")
})

test_that("n_workers() uses all cores but two, and at least one", {
  for (cores in list(NA_integer_, 1L, 2L, 3L, 10L)) {
    local_mocked_bindings(detectCores = function(...) cores, .package = "parallel")
    expected <- if (is.na(cores)) 1L else max(1L, cores - 2L)
    expect_identical(n_workers(), expected)
  }
})
