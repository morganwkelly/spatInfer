# spatInfer 0.1.0.9000

## New features

* `combined_table()` summarises several placebo and synthetic outcome tests
  in one table: an HC row, then for each placebo the row with a chosen number
  of clusters or SCPC `avc` value, with the synthetic outcome p value of the
  matching synth. Its notes give the structure, effective range and R2 of the
  treatment and of the outcome separately.

* `iv_placebo()` and `iv_synth()` run the placebo and synthetic outcome tests
  for instrumental variable regressions with one endogenous treatment and one
  instrument, written as in fixest: `outcome ~ controls | treatment ~
  instrument`. Each gives estimated and simulated p values and confidence
  intervals for 2SLS with HC, BCH and IM inference, and for the
  Anderson-Rubin (AR) test with HC and BCH inference (IM is not valid for
  AR), for 3 up to `max_clus` clusters. AR confidence sets are
  computed exactly and can be unbounded when the instrument is weak. In
  `iv_placebo()` the placebo instrument is generated from the placebo treatment
  through the estimated relationship between the instrument and the treatment,
  plus its own spatial noise. Each HC and BCH row of the results gives the
  first stage F statistic with its own inference and number of clusters; IM
  rows give none.
  `placebo_table()` and `synth_table()` show the results in separate BCH, IM
  and AR tables, chosen with their new `type` argument, and `optimal_basis()`
  and `plot_basis()` accept IV formulas.

* `iv_regression()` estimates the spatial basis IV regression by 2SLS with
  BCH or HC standard errors (`cov`), and `iv_im_regression()` with IM
  inference, using the same clusters as `iv_placebo()`. `iv_regression()`
  also gives the Anderson-Rubin p value and confidence set of the treatment
  effect and the first stage F statistic; `iv_im_regression()` gives neither. `im_table()` shows
  `iv_im_regression()` results.

* `iv_regression_table()` and `iv_regression_plot()` compare an original IV
  regression with `iv_regression()` and `iv_im_regression()` results. The
  table adds rows for the first stage F statistic and the AR p value and 95%
  set of each column except IM. The plot shows each coefficient's 2SLS
  confidence interval and, for the treatment, the AR set of each column except IM, with arrows for unbounded
  parts.

* New dataset `autor_china`: Autor, Dorn and Hanson (2013) data on Chinese
  import competition and US commuting zones.

* `im_regression()` runs the spatial basis regression with IM inference on
  every coefficient: each is estimated separately in k-medoids clusters (the
  same clusters `basis_regression()` uses for BCH) and tested with a t-test of
  its cluster estimates. Standard errors are pseudo standard errors, one
  quarter of the width of the 95% confidence interval, and t statistics are
  the estimates divided by these. `Cluster_Coefs` has one row per cluster and
  one column per coefficient. `im_table()` gives a coefficient table of the results
  (`type = "coefs"`) or the estimate of each coefficient in each cluster
  (`type = "clusters"`, rounded to `decimals` places, 2 by default).

* `scpc_regression()` runs the spatial basis regression with SCPC inference on
  every coefficient for a chosen value of `avc`.

* `regression_table()` makes a modelsummary table comparing BCH, IM and SCPC
  coefficients and standard errors, with notes giving the number of BCH and IM
  clusters and the SCPC `avc`. An optional `orig` argument adds a first
  column, "Original", for a `fixest::feols()` regression such as the original
  regression without a spatial basis.

* `placebo_scpc()` and `synth_scpc()` run the placebo and synthetic outcome
  tests with SCPC inference (Mueller and Watson 2022, 2023), ported from the
  scpcR package. Instead of numbers of clusters they examine values of the average
  correlation bound `avc`, by default 0.02, 0.04, 0.06, 0.08 and 0.1. Results
  include a `Pseudo SE` column, one quarter of the width of the 95% confidence
  interval in the units of the coefficient, and an `Estimates` component with
  the full SCPC output for each `avc`. `placebo_table()` and `synth_table()`
  accept their output; the tables leave out `Pseudo SE`.

* SCPC is computed within spatInfer, from code ported from the scpcR package
  (MIT licence), so scpcR is not needed. Tests check the results against saved
  output of scpcR 0.1.3. The spatial weights, which depend only on the
  locations and `avc`, are computed once and reused in every simulation, making
  `placebo_scpc()` and `synth_scpc()` about five times faster. geodist is a new
  dependency.

## Behaviour changes

* `weights` now takes the name of the weighting variable, for example
  `weights = "pop"`, in every function that has it: the placebo, synth and
  regression functions, including the IV ones. The variable no longer has to
  be called `weights`, and `weights = TRUE` is an error. Leave `weights` out
  for an unweighted regression. The weighting variable must be numeric with no
  missing values. Results are unchanged.

* Confidence intervals and their widths in the results and tables of all
  placebo and synthetic outcome functions are now in the units of the
  treatment coefficient. Previously they were divided by the absolute value of
  the coefficient. They are given to 3 significant digits instead of 2
  decimal places, so that intervals of small coefficients do not round to
  zero. p values are unchanged.

* `synth()` now runs in parallel by default (`Parallel = TRUE`), like the
  other placebo and synthetic outcome functions. Results are unchanged.

* Parallel runs no longer stall on Linux systems with a multithreaded BLAS
  such as OpenBLAS, the Ubuntu default. Each worker started a BLAS thread per
  core, and the oversubscribed workers ran so slowly that `placebo()`,
  `placebo_im()`, `synth_im()` and the SCPC functions appeared to hang. Workers
  now use one BLAS thread, and the user's setting is restored afterwards.
  RhpcBLASctl is a new dependency.

* The `jitter_coords` argument of `placebo()`, `placebo_im()`, `placebo_scpc()`,
  `synth()`, `synth_im()` and `synth_scpc()` is removed. The Moran test no
  longer jitters identical coordinates, which it does not need, so with
  duplicated sites its value differs slightly from before.

* `basis_regression()` now applies regression weights. Previously the
  argument was ignored and the regression was unweighted. Unweighted results
  are unchanged.

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

* `placebo_im()` and `synth_im()` now warn when the treatment coefficient
  cannot be estimated in some clusters of the regression on the actual data.
  Those clusters are left out of the IM t-test, as before; previously this
  happened without notice. The IM t-test no longer passes an
  `na.action = na.fail()` argument, which `t.test()` ignored. Results are
  unchanged.

* `plot_basis()` now uses its `theta` and `phi` arguments. Previously the view
  was fixed at the default values.

* `plot_basis()` restores the graphics parameters it changes.
