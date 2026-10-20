# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`spatInfer` is an R package implementing Conley and Kelly (2024): spatial basis regressions with large-cluster inference, plus spatial-noise placebo and synthetic-outcome diagnostics. It supports three estimators: BCH (Bester-Conley-Hansen), IM (Ibragimov-Mueller) and SCPC (Mueller-Watson).

The user workflow is `optimal_basis()` → `placebo*()` → `synth*()` → `basis_regression()` / `im_regression()` / `scpc_regression()`, with `placebo_table()`, `synth_table()`, `combined_table()` (chosen rows of several placebos with their synth p values), `im_table()` and `regression_table()` for output. Instrumental variable regressions use `iv_placebo()` and `iv_synth()`, which report 2SLS and Anderson-Rubin tests, then `iv_regression()`.

## Commands

Run from the package root:

```r
devtools::document()                 # regenerate man/ and NAMESPACE (roxygen; never edit those by hand)
devtools::test()                     # full test suite
devtools::test(filter = "scpc")      # one file: tests/testthat/test-scpc.R
testthat::test_file("tests/testthat/test-golden.R")
devtools::check()
devtools::build_readme()             # README.md is generated from README.Rmd; edit the .Rmd
```

To regenerate golden fixtures, run `Rscript tests/testthat/fixtures/make-golden.R [case ...]`. Do this only when a change in results is intended and explained. Case names are the keys of `golden_cases` in `helper-golden.R`, which match the `.rds` basenames. `scpcR_reference.rds` is not one of them. It was recorded from scpcR itself, which is no longer a dependency, so it can't be regenerated and must never be overwritten.

## Architecture

**Placebo vs synth.** Both kinds of test have the same pipeline, implemented in `R/spatial_helpers.R`:
1. `prepare_spatial_data()` validates input (coordinates must be named `X`/`Y`). Every function's `weights` argument names the weighting column as a string, or is `NULL`. `regression_weights()` checks it and copies it to `wts`, which all fits use. `prepare_spatial_data()` then renames the outcome to `dep_var` and the first RHS variable (the treatment) to `explan_var`, then appends the principal components of an `mgcv` tensor-spline basis (`prin_comp()`).
2. `build_formulas()` builds `eq_est` and `eq_sim`. In `eq_sim`, simulated noise `sim` replaces the treatment (placebo, `sim_replaces="explan_var"`) or the outcome (synth, `"dep_var"`).
3. `Noise_Sim()` fits a Matérn model to the residuals by an MLE range search with `fields::Krig`. It then draws `nSim` noise columns using an exact Cholesky decomposition, or `BRISC` when `exact_cholesky=FALSE`. The placebo takes residuals from regressing the treatment on the PCs. The synth takes residuals from regressing the outcome on a quadratic in X/Y.
4. Estimator-specific p-values come from `bch_helpers.R`, `im_helpers.R` or `scpc_helpers.R`, with HC p-values from `hc_p_values()`. They are summarised into a `Results` data frame by `summary_bch` / `summary_im` / `summary_scpc`.

The BCH and IM variants loop over k-medoids clusters from `generate_clusters()`: `cluster::pam`, or `clara` when `k_medoids=FALSE`, for 2..`max_clus` clusters. The SCPC variants loop over `avc` values instead of cluster counts.

**IV** (`R/iv_helpers.R`, `iv_placebo.R`, `iv_synth.R`) handles one endogenous treatment and one instrument, with fixest formulas `y ~ controls | x ~ z`, renamed to `dep_var`, `explan_var` and `inst_var`. Each simulated draw (`iv_sim()`) tests 2SLS (`fit_*` term) with HC, BCH and IM, and Anderson-Rubin (the instrument in the reduced form) with HC and BCH only, for every cluster count. IM is not valid for AR: there are no `AR IM` rows, and `iv_im_regression()` has no `AR`/`AR_Set` components, so don't add them back. BCH reuses one fit through `fixest::pvalue(fit, cluster=)`. IV uses clusters 3..`max_clus`; there is no 2-cluster row. Every HC and BCH Results row has a first-stage F with its own inference and cluster count, so there is no single F: the squared first-stage t (`first_stage_F()`). IM rows, and `iv_im_regression()`, deliberately have no F (`NA`, shown blank in tables). AR rows (HC and BCH) share the F of their inference. `placebo_table()`/`synth_table()` split IV results by `type` ("BCH", "IM", "AR") using `iv_table_rows()` (in `R/placebo_table.R`). AR confidence sets are analytic (`ar_conf_set()`), not a grid search; a grid search gives the same sets. The variance of `y - b*x` is quadratic in b, so three fits (`ar_fits()`) recover it exactly, and the set may be two rays or the whole line. In `iv_placebo()` the placebo instrument is `a + b*x_sim + e_sim`. Here `e_sim` comes from a second `Noise_Sim()` call with its own seeds (`seeds=c(5678,567), trend=FALSE`). Keep those seeds distinct from the treatment draw. `iv_regression()` (BCH/HC) and `iv_im_regression()` (IM; class `c("iv_im_regression","iv_regression")`) share `iv_basis_regression()`. They must agree exactly with the matching `iv_placebo()` rows, and `test-iv.R` checks this. IV results go to `iv_regression_table()`/`iv_regression_plot()` (`R/iv_regression_table.R`). The generic `regression_table()` rejects them with an error pointing to those functions. It fits on the renamed data from `prepare_iv_data()`, then refits with the original treatment and instrument names so its terms read `fit_<treatment>`. The original standalone prototypes are in `../IV/`.

**SCPC** (`R/scpc_core.R`) is ported from scpcR 0.1.3 (MIT; its licence is kept in the file header). It must reproduce scpcR's results: `test-scpc.R` compares against `fixtures/golden/scpcR_reference.rds`. scpcR is no longer a dependency.

**Parallelism.** Every simulation loop goes through `run_sims()`. It uses forked `mclapply` on Linux only. macOS and Windows use a PSOCK cluster, because Accelerate BLAS is not fork-safe. Workers are pinned to one BLAS thread (`RhpcBLASctl`) and optionally one fixest thread, and the user's settings are restored on exit. Noise is generated serially before any parallel step, so parallel results must be *identical* to serial ones.

**Side-effect contract.** Simulation functions call `withr::local_preserve_seed()`, which leaves the user's RNG state untouched, and set fixed seeds internally. Plotting restores `par()`. `test-side-effects.R` enforces this.

## Tests

- `test-golden.R` holds characterization tests. Outputs must match the `.rds` files in `tests/testthat/fixtures/golden/` to a relative tolerance of 1e-10. The shared fixture is a deterministic 99-row subset of `opportunity` with small settings, defined in `helper-golden.R` (`run_sim()`, `run_scpc()`). Any numerical change in a method will break these tests. Treat that as a regression unless the change is intended.
- `optimal_basis` golden checks use a 1% tolerance for nearly saturated rows (5×5 tensor, ≥23 PCs), because those rows are BLAS-dependent.
- A golden test is skipped, not failed, when its `.rds` is missing (`expect_matches_golden()`). A case newly added to `golden_cases` tests nothing until you run `make-golden.R <case>`.
- The `placebo_brisc` golden case is skipped when BRISC (Suggests only) is not installed, so a passing suite may not have exercised `exact_cholesky=FALSE`.
- `test-parallel.R` mocks `use_fork()` to test both back ends. It is skipped on CRAN and with fewer than 4 cores. PSOCK workers (used on macOS) load the *installed* spatInfer, not the `load_all()` copy. Tests of new internal functions therefore fail with "could not find function" until the package is installed. Install into a temporary library and put it first in `R_LIBS` to avoid replacing the user's copy.
- The IV fixture (`iv_fixture_data()`, every 6th row of `autor_china` for 2000–07) has a weak instrument (F ≈ 2). Its golden files include unbounded AR sets.

## Code style

The package uses `=` for assignment, base-R pipes `|>`, and `pkg::fun()` calls (only `stats` is imported wholesale). Internal comments use `######` banners. Simulation functions return plain lists such as `Results`, `Spatial_Params` and `Estimates`, without S3 classes. Record user-visible changes in `NEWS.md`. Datasets (`opportunity`, `autor_china`) live in `data/*.rda` and are documented by roxygen in `R/opportunity.R` and `R/autor_china.R`; keep the two in step.
