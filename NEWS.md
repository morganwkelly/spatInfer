# spatInfer 0.1.0.9000

## New features

* `basis_regression()` is renamed `bch_regression()`, and `iv_regression()`
  is renamed `iv_bch_regression()`, to match `im_regression()`,
  `scpc_regression()`, `iv_im_regression()` and `iv_scpc_regression()`. The
  old names are removed.

* `placebo_synth()` runs `placebo()` and `synth()` in one command, and
  `iv_placebo_synth()` runs `iv_placebo()` and `iv_synth()`. They show a
  progress bar with the step being run and the number of simulated regressions
  done (turn it off with `progress = FALSE`). The data, clusters and SCPC
  spatial weights are prepared once for both tests, and the results, returned
  as `Placebo` and `Synth`, are the same as those of the separate functions.
  `placebo_synth_table()` shows both in one table: every row of the placebo
  table with its synthetic outcome p value added or, with `select`, only the
  chosen rows, as `combined_table()` and `iv_combined_table()` do.
  `combined_table()` and `iv_combined_table()` also take the output of
  `placebo_synth()` or `iv_placebo_synth()` in place of a placebo and a synth,
  for example `combined_table(ps, select = c(BCH = 7, IM = 5))`.

* `placebo()` and `synth()` now run the placebo and synthetic outcome tests
  with HC, BCH, IM and SCPC inference in one call, as `iv_placebo()` and
  `iv_synth()` do. Each simulated regression is estimated once for every
  inference method. `Results` has an HC row, BCH and IM rows for 3 up to
  `max_clus` clusters and, for each SCPC bound given in the new `avc`
  argument, an SCPC row, whose `Clusters` holds `avc`. The default,
  `avc = NULL`, gives no SCPC rows. `placebo_im()`, `placebo_scpc()`,
  `synth_im()` and `synth_scpc()` are removed: their results are the IM and
  SCPC rows of `placebo()` and `synth()`, which agree with them exactly. Their
  `Coefs`, `Estimates` and `Pseudo SE` outputs are dropped, and the
  `Spatial_Params` of `placebo()` now include `N`.

* `placebo_table()` and `synth_table()` show every inference method in a single
  table, with a `Clusters/avc` column when there are SCPC rows. For
  `iv_placebo()` and `iv_synth()` results the table has a block of 2SLS rows
  followed by a block of Anderson-Rubin rows; the `type` argument is removed.

* `combined_table()` now takes one `placebo()` and the matching `synth()`, and
  a named `select` like that of `iv_combined_table()`, for example
  `select = c(BCH = 7, IM = 5, SCPC = 0.04)`, instead of lists of placebos and
  synths.

* `combined_table()` summarises several placebo and synthetic outcome tests
  in one table: an HC row, then for each placebo the row with a chosen number
  of clusters or SCPC `avc` value, with the synthetic outcome p value of the
  matching synth. Its notes give the structure, effective range and R2 of the
  treatment and of the outcome separately.

* `iv_placebo()` and `iv_synth()` run the placebo and synthetic outcome tests
  for instrumental variable regressions with one endogenous treatment and one
  instrument, written as in fixest: `outcome ~ controls | treatment ~
  instrument`. Each gives estimated and simulated p values and confidence
  intervals for 2SLS and for the Anderson-Rubin (AR) test, each with HC,
  BCH, IM and, for each SCPC bound given in `avc`, SCPC inference, for 3 up
  to `max_clus` clusters. AR tests use the
  reduced form approach of Chernozhukov and Hansen (2008): the treatment
  effect `b` is rejected when the instrument is significant in the regression
  of `outcome - b*treatment` on the instrument, controls and spatial basis,
  under the chosen inference. The IM version estimates this regression in
  each cluster, and needs no strong instrument within clusters. AR confidence
  sets are computed exactly from the reduced forms of the outcome and
  treatment, and can be unbounded when the instrument is weak. In
  `iv_placebo()` the placebo instrument is generated from the placebo treatment
  through the estimated relationship between the instrument and the treatment,
  plus its own spatial noise. Each HC and BCH row of the results gives the
  first stage F statistic with its own inference and number of clusters; IM
  and SCPC rows give none.
  `placebo_table()` and `synth_table()` show the results, and `optimal_basis()`
  and `plot_basis()` accept IV formulas.

* `iv_regression()` estimates the spatial basis IV regression by 2SLS with
  BCH or HC standard errors (`cov`), and `iv_im_regression()` with IM
  inference, using the same clusters as `iv_placebo()`. `iv_regression()`
  also gives the Anderson-Rubin p value and confidence set of the treatment
  effect and the first stage F statistic; `iv_im_regression()` gives the
  Anderson-Rubin results with IM inference, and no F statistic. `im_table()` shows
  `iv_im_regression()` results. `iv_scpc_regression()` gives 2SLS and
  Anderson-Rubin results with SCPC inference: 2SLS is treated as least
  squares on the instrumented regressors, with structural residuals.

* `iv_combined_table()` is the IV version of `combined_table()`: from one
  `iv_placebo()` and the matching `iv_synth()` it shows the HC row and the
  chosen BCH, IM and SCPC rows, with their synthetic outcome p values, in a
  2SLS block and an Anderson-Rubin block (`select`, e.g.
  `c(BCH=7, IM=5, SCPC=0.04)`). `ar_select` gives the AR rows their own
  inference methods, numbers of clusters or `avc`, and defaults to `select`.

* `iv_regression_table()` and `iv_regression_plot()` compare an original IV
  regression with `iv_regression()`, `iv_im_regression()` and
  `iv_scpc_regression()` results. The table adds rows for the first stage F
  statistic and the AR p value and 95% set of each column (no F for IM or
  SCPC). The plot shows each coefficient's 2SLS
  confidence interval and, for the treatment, the AR set of each column, with arrows for unbounded
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

* `synth()`, `synth_im()`, `synth_scpc()` and `iv_synth()` have an
  `xy_interaction` argument. Set to `TRUE`, it adds the product of longitude
  and latitude to the quadratic trend that generates the synthetic outcome.
  The default, `FALSE`, keeps the previous trend.

## Behaviour changes

* `basis_regression()`, `im_regression()`, `iv_regression()` and
  `iv_im_regression()` have a `k_medoids` argument. Set it to `FALSE` when the
  placebo test used `k_medoids = FALSE`, so that the final regression uses the
  same Clara clusters. Previously they always used PAM clusters, which could
  differ from the clusters the placebo test chose between. Results with the
  default are unchanged.

* The treatment's p values and confidence intervals are now found by name, not
  by their position in the coefficient table. If the treatment is dropped as
  collinear, an error says so, instead of results being reported for another
  variable. Results are otherwise unchanged.

* `placebo_scpc()` and `synth_scpc()` leave out simulations where SCPC fails
  separately for each `avc`, instead of leaving a simulation out for every
  `avc` when it fails for one.

* The placebo and synthetic outcome functions warn when the estimated range of
  spatial correlation is at the upper limit of its search, the 95th percentile
  of distances between points, since the true range may then be longer.

* The BCH simulations of `placebo()` and `synth()` estimate each simulated
  regression once and compute its clustered p value for every number of
  clusters from that fit, instead of re-estimating it for each. Results are
  unchanged.

* Weighted SCPC: every function with SCPC inference has a `scpc_weighting`
  argument. The default, `"scpcR"`, conditions on the unweighted regressors,
  as scpcR does, so weighted results match scpcR. `scpc_weighting = "sqrt"`
  instead computes conditional SCPC for the regression rescaled by the square
  root of the weights, the least squares regression equivalent to weighted
  least squares. The choice changes only the SCPC p values and confidence
  intervals of weighted regressions.

* SCPC inference on a weighted regression in which fixest drops observations
  with missing values no longer fails: the weights of the dropped observations
  are left out.

* The help pages of `placebo_scpc()`, `synth_scpc()` and `scpc_regression()`
  note that with 4,500 or more observations SCPC, as in scpcR, uses distances
  in degrees where great-circle distances are intended. The calculation is
  unchanged.

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
