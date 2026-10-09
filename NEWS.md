# spatInfer (development version)

## Behaviour changes

* `basis_regression()` now applies regression weights when `weights = TRUE`.
  Previously the argument was ignored and the regression was unweighted. It
  now also stops with an error if there is no `weights` variable, as the other
  functions do. Unweighted results are unchanged.
