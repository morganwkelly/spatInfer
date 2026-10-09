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

test_that("KNOWN ISSUE: optimal_basis(max_splines = 3) adds 4x4 and duplicates 3x3", {
  expect_golden("optimal_basis_max3")
  ob3 <- readRDS(golden_path("optimal_basis_max3"))
  expect_setequal(unique(ob3$bic$name), c("BIC_3.x", "BIC_4", "BIC_3.y"))
})
