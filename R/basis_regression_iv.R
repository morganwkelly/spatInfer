#' Carry out spatial basis IV (2SLS) regressions with BCH standard errors.
#'
#' @description
#' `basis_regression_iv()` runs two-stage least squares (2SLS) with a spatial
#' basis added as exogenous controls, and large-cluster BCH or heteroskedasticity-
#' consistent standard errors. It mirrors [basis_regression()] in interface and
#' output but accepts an IV formula specifying an endogenous variable and one or
#' more instruments.
#'
#' The spatial basis principal components are added to the exogenous part of the
#' first-stage and reduced-form regressions, de-trending the data and reducing
#' residual spatial correlation before inference.
#'
#' @param fm The IV formula in the form `y ~ exog_controls | endog ~ instruments`.
#'   `y` is the outcome, `exog_controls` are additional exogenous regressors
#'   (use `1` if none), `endog` is the endogenous variable of interest, and
#'   `instruments` are the excluded instruments. Multiple instruments or controls
#'   are separated by `+` in the usual way, e.g.
#'   `y ~ ctrl1 + ctrl2 | endog ~ z1 + z2`.
#' @param df The dataset. Longitude and latitude must be named `X` and `Y` and
#'   have no missing values.
#' @param splines The dimension of the linear tensor spline basis, from
#'   [optimal_basis()].
#' @param pc_num The number of principal components to include, from
#'   [optimal_basis()].
#' @param clusters Number of k-medoids clusters for BCH standard errors, chosen
#'   using the placebo test from [placebo()].
#' @param weights Set `weights=TRUE` if the regression is weighted. The
#'   weighting variable in the dataset must be named `wts`.
#' @param cov Defaults to `"BCH"` for large-cluster bias-corrected
#'   heteroskedasticity-consistent standard errors. Any other value gives plain
#'   heteroskedasticity-consistent standard errors.
#'
#' @return A `fixest` object (from [fixest::feols()]) that can be printed and
#'   exported using the `modelsummary` package. The first-stage F-statistic is
#'   accessible via `fixest::fitstat(result, "ivf")`. Because t-statistics from
#'   large-cluster regressions have low degrees of freedom, reporting confidence
#'   intervals and p-values (rather than standard errors) is recommended in
#'   `modelsummary` output.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' set.seed(123)
#' opportunity <- opportunity |> dplyr::slice_sample(n=250)
#' # Instrument single_mothers with dropout_rate, controlling for gini.
#' # Use splines/PCs from optimal_basis and clusters from placebo.
#' ck_iv <- basis_regression_iv(
#'   mobility ~ gini | single_mothers ~ dropout_rate,
#'   opportunity,
#'   splines = 4, pc_num = 3, clusters = 4
#' )
#'
#' # First-stage F-statistic
#' fixest::fitstat(ck_iv, "ivf")
#'
#' # Regression table (omit basis PCs and intercept)
#' modelsummary::modelsummary(list(`IV-BCH` = ck_iv),
#'   statistic = c("conf.int", "p = {p.value}"),
#'   coef_omit = "Intercept|PC",
#'   gof_map = c("nobs", "r.squared"), fmt = 2)

basis_regression_iv <- function(fm, df, splines, pc_num, clusters,
                                 weights = FALSE, cov = "BCH") {

  if (is.null(df$X) | is.null(df$Y))
    stop("You must have longitude and latitude variables named X and Y")
  if (sum(is.na(df$X)) > 0 | sum(is.na(df$Y)) > 0)
    stop("You cannot have missing values in longitude and latitude.")

  # Parse IV formula: y ~ exog_controls | endog ~ instruments
  iv_names <- set_names_iv(fm)
  var_dep  <- iv_names$var_dep   # outcome variable name
  rhs      <- iv_names$rhs       # exogenous controls (raw string, no leading "+")
  var_expl <- iv_names$var_expl  # endogenous variable name
  var_inst <- iv_names$var_inst  # excluded instrument(s) name

  # Rename outcome to dep_var so that prin_comp() can call bam(dep_var ~ te(X,Y,...))
  df1 <- df |> dplyr::rename(dep_var = iv_names$var_dep)

  # Generate spatial basis principal components and add to dataset
  pc     <- prin_comp(df1, splines, pc_num)
  df1    <- cbind.data.frame(df1, pc)
  pc_str <- paste(names(pc), collapse = "+")

  # Build the exogenous RHS: user-supplied controls (if any) plus spatial PCs
  exog_part <- if (is.na(rhs) || rhs == "") {
    pc_str
  } else {
    paste0(rhs, "+", pc_str)
  }

  # Full IV formula passed to feols:  dep_var ~ exog_controls + PCs | endog ~ instruments
  eq_iv <- as.formula(paste0("dep_var~", exog_part, "|", var_expl, "~", var_inst))

  wts_formula <- if (weights) ~wts else NULL

  if (cov == "BCH") {
    Coords    <- as.matrix(df |> dplyr::select(X, Y))
    clust_bch <- factor(cluster::pam(Coords, k = clusters)$clustering)
    CK <- fixest::feols(eq_iv,
                        data     = df1,
                        weights  = wts_formula,
                        cluster  = clust_bch,
                        data.save = TRUE)
  } else {
    CK <- fixest::feols(eq_iv,
                        data     = df1,
                        weights  = wts_formula,
                        vcov     = "hetero",
                        data.save = TRUE)
  }

  return(CK)
}
