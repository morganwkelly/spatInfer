# Known issues, pinned so that fixing one is a deliberate, visible change.
# When an issue is fixed, replace its test with one for the corrected behaviour.

test_that("KNOWN ISSUE: jitter_coords has no effect", {
  with_jitter <- readRDS(golden_path("placebo"))
  without_jitter <- run_sim(placebo, Parallel = FALSE, jitter_coords = FALSE)
  expect_identical(without_jitter, with_jitter)
})
