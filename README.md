
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatInfer

<!-- badges: start -->

<!-- badges: end -->

Strong autocorrelation and directional trends mean that regressions
using spatial observations are prone to return inflated t-statistics
regardless of the validity of their identification procedure. The
purpose of `spatInfer` is to estimate reliable results from regressions
using spatial observations by implementing the diagnostic and regression
procedures introduced by Conley and Kelly (2024). A simple workflow is
used to carry out spatial noise placebo and synthetic outcome tests, and
estimate regressions.

These are applied to three estimators: Bester-Conley-Hansen,
Ibragimov-Mueller, and Mueller-Watson. If these very different
procedures give similar results you can be reasonably confident that
your results are reliable.

## Workflow.

`spatInfer` is based on a sequence of four steps.

1.  `optimal_basis()` gives the spatial basis to include in regressions.

2.  `placebo()`, `placebo_im()` and `placebo_scpc()` estimate treatment
    significance from spatial noise placebos and select the tuning
    parameter of each estimator: the number of clusters for
    Bester-Conley-Hansen and Ibragimov-Mueller, and the bound on average
    correlation `avc` for Mueller-Watson (SCPC). `placebo_table()`
    displays the results.

3.  `synth()`, `synth_im()` and `synth_scpc()` test the null hypothesis
    that the outcome is trending spatial noise, and therefore
    independent of the treatment. `synth_table()` displays the results,
    and `combined_table()` shows the chosen row of each placebo test
    with its synthetic outcome p value.

4.  `basis_regression()`, `im_regression()` and `scpc_regression()`
    estimate the regression with each procedure, and
    `regression_table()` and `regression_plot()` display them along with
    the original, uncorrected estimate.

The final regressions should use the same settings as the placebo tests.
In particular, if the placebo used `k_medoids = FALSE` (Clara clusters,
for large datasets), pass `k_medoids = FALSE` to the regression so that
its clusters are the same. By default the synthetic outcome is generated
by a quadratic in longitude and latitude; `xy_interaction = TRUE` adds
their product.

## Instrumental variables.

Regressions with one endogenous treatment and one instrument, written as
in fixest, `outcome ~ controls | treatment ~ instrument`, follow the
same workflow with `iv_placebo()`, `iv_synth()`, and then
`iv_regression()`, `iv_im_regression()` and `iv_scpc_regression()`.
These report 2SLS estimates and Anderson-Rubin tests, computed by the
reduced form approach of Chernozhukov and Hansen (2008), with BCH, IM or
SCPC inference, along with the first stage F statistic. Display the
results with `iv_regression_table()` and `iv_regression_plot()`.

## Tutorial

A step by step tutorial can be found
[here](https://morganwkelly.github.io/spatInfer_tutor/)

## Installation

You can install the development version of spatInfer from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("morganwkelly/spatInfer")
```
