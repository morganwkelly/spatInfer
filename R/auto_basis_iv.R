#' Automated spatial IV regression with optimal basis and cluster selection
#'
#' @description
#' `auto_basis_iv()` runs the full three-step spatial IV pipeline automatically:
#'
#' 1. **Optimal basis** — sweeps tensor spline dimensions from 3×3 to
#'    `max_splines`×`max_splines`, extracts principal components at each
#'    dimension, and chooses the (spline, PC count) combination that minimises
#'    BIC for the outcome variable. This replicates the logic of
#'    [optimal_basis()] but returns values instead of a plot.
#'
#' 2. **Optimal clusters** — runs [placebo()] using the endogenous variable as
#'    the "treatment" (so that its spatial structure drives the noise
#'    simulations) and automatically selects the cluster count whose placebo
#'    rejection rate at the 5 % level is closest to 0.05.
#'
#' 3. **IV regression** — calls [basis_regression_iv()] with the chosen
#'    parameters and returns its output together with the selected values and
#'    the full placebo diagnostics.
#'
#' @param fm The IV formula in the form `y ~ exog_controls | endog ~ instruments`.
#'   See [basis_regression_iv()] for details.
#' @param df The dataset. Must contain columns `X` (longitude) and `Y`
#'   (latitude) with no missing values.
#' @param max_splines Maximum tensor dimension to examine (3 to 12). Defaults
#'   to 8.
#' @param nSim Number of placebo simulations for cluster selection. Defaults to
#'   1000; use a smaller value (e.g. 100) for exploratory runs.
#' @param weights Set `weights=TRUE` if the regression is weighted. The
#'   weighting variable must be named `wts` in `df` for
#'   [basis_regression_iv()] and `weights` in `df` for [placebo()].
#' @param max_clus Maximum number of clusters to examine in the placebo test.
#'   Defaults to 6.
#' @param Parallel Run simulations in parallel. Set to `FALSE` if there are
#'   memory problems.
#' @param exact_cholesky For very large datasets, set to `FALSE` to use the
#'   BRISC Cholesky approximation.
#' @param k_medoids For large datasets, set to `FALSE` to use CLARA instead of
#'   PAM for medoid computation.
#' @param jitter_coords If some sites share identical coordinates, jitter by
#'   Gaussian noise (sd = 0.01) so the Moran test can be computed.
#' @param max_pc Maximum number of principal components to consider when
#'   searching for the BIC-optimal basis. Defaults to 6. Increase if you
#'   suspect the spatial trend requires more components; decrease to keep the
#'   regression parsimonious.
#' @param cov Covariance estimator passed to [basis_regression_iv()]. Defaults
#'   to `"BCH"`.
#' @param verbose If `TRUE` (default), prints the chosen splines, PC count, and
#'   cluster number to the console.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`result`}{A `fixest` object from [basis_regression_iv()]. Can be
#'     passed directly to `modelsummary::modelsummary()` and
#'     `fixest::fitstat()`.}
#'   \item{`splines`}{Chosen tensor spline dimension.}
#'   \item{`pc_num`}{Chosen number of principal components.}
#'   \item{`clusters`}{Chosen number of BCH clusters.}
#'   \item{`placebo`}{Full list returned by [placebo()], containing `Results`
#'     and `Spatial_Params`.}
#' }
#'
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' set.seed(123)
#' opp <- opportunity |> dplyr::slice_sample(n = 150)
#'
#' res <- auto_basis_iv(
#'   mobility ~ gini | single_mothers ~ dropout_rate,
#'   opp,
#'   max_splines = 6, nSim = 100, Parallel = FALSE
#' )
#'
#' # Chosen parameters
#' res$splines; res$pc_num; res$clusters
#'
#' # First-stage F-statistic
#' fixest::fitstat(res$result, "ivf")
#'
#' # Regression table
#' modelsummary::modelsummary(list(`IV-BCH` = res$result),
#'   statistic = c("conf.int", "p = {p.value}"),
#'   coef_omit = "Intercept|PC", fmt = 2)

auto_basis_iv <- function(fm, df, max_splines = 8, nSim = 1000,
                           weights = FALSE, max_clus = 6, max_pc = 6,
                           Parallel = TRUE, exact_cholesky = TRUE, k_medoids = TRUE,
                           jitter_coords = TRUE, cov = "BCH", verbose = TRUE) {

  if (is.null(df$X) | is.null(df$Y))
    stop("You must have longitude and latitude variables named X and Y")
  if (sum(is.na(df$X)) > 0 | sum(is.na(df$Y)) > 0)
    stop("You cannot have missing values in longitude and latitude.")
  if (max_splines > 12)
    stop("The maximum number of splines you can use is 12.")
  if (max_clus < 3)
    stop("Your maximum number of clusters max_clus must be greater than 2.")

  # ── Step 1: parse IV formula ─────────────────────────────────────────────────
  iv_names <- set_names_iv(fm)
  var_dep  <- iv_names$var_dep   # outcome
  rhs      <- iv_names$rhs       # exogenous controls string (may be "1" or "")
  var_expl <- iv_names$var_expl  # endogenous variable
  var_inst <- iv_names$var_inst  # excluded instrument(s)

  # ── Step 2: BIC loop over spline dimensions ───────────────────────────────────
  # Rename outcome to dep_var so mgcv::bam() can find it (same convention as
  # optimal_basis() and prin_comp()).
  df_tmp <- df |> dplyr::rename(dep_var = !!var_dep)
  mx <- max_splines - 2   # number of dimensions above 3x3

  bic_rows <- list()
  idx <- 1L

  for (spl in 3:max_splines) {
    gm <- mgcv::bam(dep_var ~ te(X, Y, bs = "bs", k = spl, m = 1),
                    data = df_tmp, discrete = TRUE)
    pc_all <- prcomp(model.matrix(gm))
    pc_df  <- cbind.data.frame(dep_var = df_tmp$dep_var, pc_all$x)
    n_pc   <- ncol(pc_df) - 1L

    for (j in seq_len(min(n_pc, max_pc))) {
      ll <- lm(dep_var ~ ., pc_df[, 1:(j + 1L)])
      bic_rows[[idx]] <- data.frame(BIC = BIC(ll), spline = spl, pc = j)
      idx <- idx + 1L
    }
  }

  all_bic      <- purrr::list_rbind(bic_rows)
  best_row     <- all_bic[which.min(all_bic$BIC), ]
  best_splines <- best_row$spline
  best_pc_num  <- best_row$pc

  # ── Step 3: placebo test for optimal cluster count ────────────────────────────
  # Use the endogenous variable as the "treatment" so that the spatial structure
  # of the variable being instrumented drives the noise simulations.
  if (rhs == "1" || rhs == "") {
    fm_ols <- as.formula(paste0(var_dep, "~", var_expl))
  } else {
    fm_ols <- as.formula(paste0(var_dep, "~", var_expl, "+", rhs))
  }

  plbo <- placebo(fm_ols, df,
                  splines        = best_splines,
                  pc_num         = best_pc_num,
                  nSim           = nSim,
                  weights        = weights,
                  max_clus       = max_clus,
                  Parallel       = Parallel,
                  exact_cholesky = exact_cholesky,
                  k_medoids      = k_medoids,
                  jitter_coords  = jitter_coords)

  # Auto-select: BCH row with placebo rejection rate closest to 0.05
  bch_rows      <- plbo$Results |> dplyr::filter(SE == "BCH")
  best_clusters <- as.integer(
    bch_rows$Clusters[which.min(abs(bch_rows$`Plac 5%` - 0.05))]
  )

  # ── Step 4: IV regression ────────────────────────────────────────────────────
  if (verbose) {
    plac_rate <- bch_rows$`Plac 5%`[bch_rows$Clusters == as.character(best_clusters)]
    message(sprintf("Optimal basis  : %dx%d tensor with %d PCs",
                    best_splines, best_splines, best_pc_num))
    message(sprintf("Optimal clusters: %d  (placebo rejection rate at 5%%: %.1f%%)",
                    best_clusters, 100 * plac_rate))
  }

  result <- basis_regression_iv(fm, df,
                                splines  = best_splines,
                                pc_num   = best_pc_num,
                                clusters = best_clusters,
                                weights  = weights,
                                cov      = cov)

  # ── Return ───────────────────────────────────────────────────────────────────
  list(
    result   = result,
    splines  = best_splines,
    pc_num   = best_pc_num,
    clusters = best_clusters,
    placebo  = plbo
  )
}
