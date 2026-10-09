# Known issues, pinned so that fixing one is a deliberate, visible change.
# When an issue is fixed, replace its test with one for the corrected behaviour.

test_that("KNOWN ISSUE: jitter_coords has no effect", {
  with_jitter <- readRDS(golden_path("placebo"))
  without_jitter <- run_sim(placebo, Parallel = FALSE, jitter_coords = FALSE)
  expect_identical(without_jitter, with_jitter)
})

test_that("KNOWN ISSUE: tables pad NA rows when max_clus < 5", {
  plbo <- run_sim(placebo, Parallel = FALSE, max_clus = 3)
  expect_equal(nrow(plbo$Results), 2)
  plac_tab <- placebo_table(plbo)
  expect_equal(nrow(plac_tab@data), 5)
  expect_true(all(is.na(plac_tab@data$Clusters[3:5])))
})
