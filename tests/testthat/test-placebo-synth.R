# placebo_synth() and iv_placebo_synth() share the data, clusters and SCPC spatial weights
# between the two tests, and must give exactly the results of the separate functions.

test_that("placebo_synth() matches placebo() and synth()", {
  ps <- run_sim(placebo_synth, Parallel = FALSE, avc = fixture_avc)
  expect_named(ps, c("Placebo", "Synth"))
  expect_matches_golden(ps$Placebo, "placebo")
  expect_matches_golden(ps$Synth, "synth")

  # Every option reaches both tests.
  args <- list(Parallel = FALSE, avc = fixture_avc, weights = "wt", data = fixture_data(TRUE),
               max_clus = 4, xy_interaction = TRUE, k_medoids = FALSE, scpc_weighting = "sqrt",
               progress = FALSE)
  ps <- do.call(run_sim, c(list(placebo_synth), args))
  args$progress <- NULL
  expect_identical(ps$Placebo, do.call(run_sim, c(list(placebo), args[names(args) != "xy_interaction"])))
  expect_identical(ps$Synth, do.call(run_sim, c(list(synth), args)))
})

test_that("iv_placebo_synth() matches iv_placebo() and iv_synth()", {
  ps <- run_iv(iv_placebo_synth, Parallel = FALSE, avc = iv_avc)
  expect_named(ps, c("Placebo", "Synth"))
  expect_matches_golden(ps$Placebo, "iv_placebo")
  expect_matches_golden(ps$Synth, "iv_synth")

  ps <- run_iv(iv_placebo_synth, Parallel = FALSE, weights = "timepwt48", avc = iv_avc,
               xy_interaction = TRUE, progress = FALSE)
  expect_matches_golden(ps$Placebo, "iv_placebo_weighted_scpc")
  expect_identical(ps$Synth, run_iv(iv_synth, Parallel = FALSE, weights = "timepwt48", avc = iv_avc,
                                    xy_interaction = TRUE))
})

test_that("the progress bar counts every simulated regression", {
  steps <- character()
  counts <- integer()
  local_mocked_bindings(sim_progress = function(total, .envir) {
    expect_identical(total, 12)
    list(status = function(msg) steps <<- c(steps, msg), add = function(k) counts <<- c(counts, k))
  })
  run_sim(placebo_synth, Parallel = FALSE, nSim = 6)
  expect_identical(steps, c("Placebo: fitting spatial noise", "Placebo: simulations",
                            "Synth: fitting spatial noise", "Synth: simulations"))
  expect_identical(sum(counts), 12)

  steps <- character()
  counts <- integer()
  run_iv(iv_placebo_synth, Parallel = FALSE, nSim = 6)
  expect_identical(steps, c("IV placebo: fitting spatial noise", "IV placebo: simulations",
                            "IV synth: fitting spatial noise", "IV synth: simulations"))
  expect_identical(sum(counts), 12)

  steps <- character()
  run_sim(placebo_synth, Parallel = FALSE, nSim = 6, progress = FALSE)
  expect_identical(steps, character())
})

test_that("placebo_synth() leaves the user's random number state unchanged", {
  set.seed(42)
  before <- .Random.seed
  run_sim(placebo_synth, Parallel = FALSE, nSim = 5)
  expect_identical(.Random.seed, before)
})

test_that("placebo_synth_table() adds the synth p values to every placebo row", {
  ps <- run_sim(placebo_synth, Parallel = FALSE, avc = fixture_avc, progress = FALSE)
  tab <- placebo_synth_table(ps)
  plac <- placebo_table(ps$Placebo)@data
  expect_identical(names(tab@data), append(names(plac), "Synth p", after = which(names(plac) == "Plac 5%")))
  expect_identical(tab@data[names(plac)], plac)
  expect_equal(tab@data$`Synth p`, round(ps$Synth$Results$`Synth p`, 3))
  expect_match(unlist(tab@notes)[2], "Treatment: Structure=.*Outcome: Structure=")

  # With select it is the combined table.
  sel <- c(BCH = 4, SCPC = 0.04)
  expect_identical(placebo_synth_table(ps, select = sel)@data,
                   combined_table(ps$Placebo, ps$Synth, select = sel)@data)
  # combined_table() also takes the placebo_synth() output itself.
  expect_identical(combined_table(ps, select = sel)@data,
                   combined_table(ps$Placebo, ps$Synth, select = sel)@data)
  expect_identical(combined_table(ps, sel)@data, combined_table(ps, select = sel)@data)
  expect_error(combined_table(ps, ps$Synth, sel), "must be left out")
  expect_error(combined_table(ps), "named vector")

  expect_error(placebo_synth_table(ps$Placebo), "not the output of placebo_synth")
  bad <- ps
  bad$Synth$Results <- bad$Synth$Results[-2, ]
  expect_error(placebo_synth_table(bad), "same inference methods")
})

test_that("placebo_synth_table() shows IV results in 2SLS and Anderson-Rubin blocks", {
  ps <- run_iv(iv_placebo_synth, Parallel = FALSE, avc = iv_avc, progress = FALSE)
  tab <- placebo_synth_table(ps)
  plac <- placebo_table(ps$Placebo)@data
  expect_identical(tab@data[names(plac)], plac)
  expect_equal(tab@data$`Synth p`, round(c(ps$Synth$Results$`Synth p`), 3))
  expect_match(paste(unlist(tab@notes), collapse = " "), "Noise instrument")

  sel <- c(BCH = 4, IM = 5)
  expect_identical(placebo_synth_table(ps, select = sel, ar_select = c(IM = 4))@data,
                   iv_combined_table(ps$Placebo, ps$Synth, select = sel, ar_select = c(IM = 4))@data)
  expect_identical(iv_combined_table(ps, select = sel, ar_select = c(IM = 4))@data,
                   iv_combined_table(ps$Placebo, ps$Synth, select = sel, ar_select = c(IM = 4))@data)
  expect_identical(iv_combined_table(ps, sel)@data,
                   iv_combined_table(ps$Placebo, ps$Synth, select = sel)@data)
  expect_error(iv_combined_table(ps, ps$Synth, sel), "must be left out")
})
