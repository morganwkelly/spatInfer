# spatInfer 0.1.0.9000

## New features

* `placebo_scpc()` and `synth_scpc()` run the placebo and synthetic outcome
  tests with SCPC inference (Mueller and Watson 2022, 2023) from the scpcR
  package. Instead of numbers of clusters they examine values of the average
  correlation bound `avc`, by default 0.02, 0.04, 0.06, 0.08 and 0.1. Results
  include a `Pseudo SE` column, one quarter of the width of the 95% confidence
  interval in the units of the coefficient, and an `Estimates` component with
  the full SCPC output for each `avc`. `placebo_table()` and `synth_table()`
  accept their output.

* scpcR (installed from GitHub with `remotes::install_github("spatial-spur/scpcR")`)
  and geodist are new dependencies. scpcR requires fixest 0.14.0 or later.

## Behaviour changes

* `basis_regression()` now applies regression weights when `weights = TRUE`.
  Previously the argument was ignored and the regression was unweighted. It
  now also stops with an error if there is no `weights` variable, as the other
  functions do. Unweighted results are unchanged.

* `placebo()`, `placebo_im()`, `synth()` and `synth_im()` no longer change the
  user's random number state. They still use the same fixed internal seeds, so
  their results are unchanged, but the caller's `.Random.seed` is restored when
  they return.

* Parallel runs no longer leave global settings changed. The fixest thread
  setting is restored afterwards, and no foreach backend is registered:
  simulations now run with `parallel::mclapply()` (Linux) or a temporary
  PSOCK cluster that is always stopped (macOS, Windows). An error in a worker is
  now raised instead of being lost. `foreach` and `doParallel` are no longer
  dependencies; `withr` is new.

* On macOS, parallel runs no longer fork. R on macOS uses Apple's Accelerate
  BLAS, which is not fork-safe, so forked workers could crash inside
  `fields::Krig()`, especially under RStudio. Previously the results of
  crashed workers were silently dropped: the Matern search could skip some
  ranges, and placebo p values could be computed from fewer than `nSim`
  simulations. Results from runs where no worker crashed are unchanged.

* `optimal_basis(max_splines = 3)` now examines only the 3x3 tensor.
  Previously its loop ran backwards, adding a 4x4 tensor and duplicating the
  3x3 results. `max_splines` below 3 is now rejected with an error. Results for
  `max_splines` of 4 or more are unchanged.

* `optimal_basis()` with `max_splines` of 11 or 12 now gives every tensor its
  own colour. The palette had eight colours, so the 11x11 and 12x12 curves
  were drawn in the default grey for missing colours. Two colours (black and
  wine) are added; plots with up to eight tensors are unchanged.

* `placebo_table()` and `synth_table()` no longer add rows of `NA` when
  `max_clus` is below 5, and show the standard error label (HC, BCH or IM) only
  in the first row of each block for any `max_clus`. Previously rows 3 to 5
  were blanked, which also left the label showing from row 6 when
  `max_clus` is 7 or more. Tables for the default `max_clus = 6` are unchanged.

* `exact_cholesky = FALSE` now stops with an informative error, before any
  simulation work, when the suggested `BRISC` package is not installed.

* Parallel runs use at least one worker. Previously `detectCores() - 2`
  workers were requested, which fails on machines with two or fewer cores or
  where the number of cores cannot be detected. Machines with three or more
  cores use the same number of workers as before.

* `jitter_coords = FALSE` now turns off jittering of identical coordinates
  in the Moran test of `placebo()`, `placebo_im()`, `synth()` and `synth_im()`.
  Previously the argument was ignored and identical coordinates were always
  jittered. The default, `TRUE`, is unchanged. The documentation now gives the
  jitter as about 1 km (0.01 degrees), not 10 km.

* `placebo_im()` and `synth_im()` now warn when the treatment coefficient
  cannot be estimated in some clusters of the regression on the actual data.
  Those clusters are left out of the IM t-test, as before; previously this
  happened without notice. The IM t-test no longer passes an
  `na.action = na.fail()` argument, which `t.test()` ignored. Results are
  unchanged.

* `plot_basis()` now uses its `theta` and `phi` arguments. Previously the view
  was fixed at the default values.

* `plot_basis()` restores the graphics parameters it changes.
