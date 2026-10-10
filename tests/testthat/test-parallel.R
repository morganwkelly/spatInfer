# Parallel execution must reproduce serial results exactly: the simulated noise
# is generated serially before any parallel step. Each parallel run is compared
# with a serial run in the same session, so the comparison can be exact. Both
# parallel back ends are tested on every platform that can fork: separate worker
# processes (used on macOS and Windows) and forked workers (used on Linux).

skip_parallel <- function() {
  skip_on_cran()
  skip_if(isTRUE(parallel::detectCores() < 4), "Fewer than 4 cores available")
}

local_backend <- function(fork, env = parent.frame()) {
  if (fork) skip_on_os("windows")
  local_mocked_bindings(use_fork = function() fork, .env = env)
}

for (fork in c(FALSE, TRUE)) {
  backend <- if (fork) "forked workers" else "worker processes"

  test_that(paste("placebo() gives identical results with", backend), {
    skip_parallel()
    local_backend(fork)
    expect_identical(run_sim(placebo, Parallel = TRUE), run_sim(placebo, Parallel = FALSE))
  })

  test_that(paste("placebo_im() gives identical results with", backend), {
    skip_parallel()
    local_backend(fork)
    expect_identical(run_sim(placebo_im, Parallel = TRUE), run_sim(placebo_im, Parallel = FALSE))
  })

  test_that(paste("placebo_scpc() gives identical results with", backend), {
    skip_parallel()
    local_backend(fork)
    expect_identical(run_scpc(placebo_scpc, Parallel = TRUE), run_scpc(placebo_scpc, Parallel = FALSE))
  })

  test_that(paste("errors are raised from", backend), {
    skip_parallel()
    local_backend(fork)
    expect_error(run_sims(4, function(j) if (j == 3) stop("worker failed") else j, TRUE, FALSE),
      "worker failed")
  })
}

test_that("forked workers are not used on macOS or Windows", {
  if (Sys.info()[["sysname"]] == "Darwin" || .Platform$OS.type == "windows") {
    expect_false(use_fork())
  } else {
    expect_true(use_fork())
  }
})

test_that("parallel runs restore the fixest thread setting", {
  skip_parallel()
  before <- fixest::getFixest_nthreads()
  run_sim(placebo, Parallel = TRUE)
  expect_identical(fixest::getFixest_nthreads(), before)
})

test_that("n_workers() uses all cores but two, and at least one", {
  for (cores in list(NA_integer_, 1L, 2L, 3L, 10L)) {
    local_mocked_bindings(detectCores = function(...) cores, .package = "parallel")
    expected <- if (is.na(cores)) 1L else max(1L, cores - 2L)
    expect_identical(n_workers(), expected)
  }
})
