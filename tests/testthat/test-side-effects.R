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
  expect_matches_golden(first, "placebo")
})

test_that("plot_basis() restores graphics parameters", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  graphics::par(mai = c(1, 1, 1, 1))
  quietly(plot_basis(fixture_fm, fixture_data(), splines = 4))
  expect_equal(graphics::par("mai"), c(1, 1, 1, 1))
})

test_that("parallel workers use one BLAS thread and the user's setting is restored", {
  old <- RhpcBLASctl::blas_get_num_procs()
  on.exit(RhpcBLASctl::blas_set_num_threads(old))
  RhpcBLASctl::blas_set_num_threads(2)
  # R's reference BLAS is single threaded and ignores the setting.
  if (RhpcBLASctl::blas_get_num_procs() != 2L) skip("BLAS thread count cannot be set")
  for (fork in c(FALSE, if (.Platform$OS.type == "unix") TRUE)) {
    local_mocked_bindings(use_fork = function() fork)
    threads <- run_sims(2, function(j) RhpcBLASctl::blas_get_num_procs(), TRUE, FALSE)
    expect_identical(unlist(threads), c(1L, 1L))
    expect_identical(RhpcBLASctl::blas_get_num_procs(), 2L)
  }
})
