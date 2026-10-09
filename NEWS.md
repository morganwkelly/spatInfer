# spatInfer 0.1.0.9000

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
  simulations now run with `parallel::mclapply()` (macOS, Linux) or a temporary
  PSOCK cluster that is always stopped (Windows). An error in a worker is now
  raised instead of being lost. `foreach` and `doParallel` are no longer
  dependencies; `withr` is new.

* `plot_basis()` restores the graphics parameters it changes.
