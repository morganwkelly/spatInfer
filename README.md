
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatInfer

<!-- badges: start -->
<!-- badges: end -->

Strong autocorrelation and directional trends mean that regressions
using spatial observations are prone to return inflated t-statistics
regardless of the validity of their identification procedure.
`spatInfer` introduces a regression procedure to obtain reliable results
in these circumstances by adding a spatial basis and using large cluster
standard errors. It generates a spatial noise placebo to ensure the
validity of the inference procedure, and a synthetic outcome test to
ascertain whether the outcome is spatial noise independent of the
treatment.

## Workflow.

`spatInfer` is based on a sequence of four commands.

1.  `optimal_basis()` gives the spatial basis to include in regressions.

2.  `placebo()` estimates treatment significance from spatial noise
    placebos and chooses the optimal number of large clusters for
    standard error estimation.

3.  `synth()` tests the null hypothesis that the outcome is trending
    spatial noise, and therefore independent of the treatment.

4.  `basis_regression()` runs the regression with the spatial basis and
    large cluster standard errors.

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
