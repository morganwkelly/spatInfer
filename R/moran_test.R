#' Moran's I test on the residuals of a fixest model
#'
#' Computes Moran's I for spatial autocorrelation in the residuals of a
#' [fixest::feols()] object, including those returned by [basis_regression()]
#' and [basis_regression_iv()]. Use this as a post-estimation diagnostic: after
#' adding a spatial basis the residuals should show little remaining
#' autocorrelation.
#'
#' @param model A `fixest` object (from [fixest::feols()], [basis_regression()],
#'   or [basis_regression_iv()]).
#' @param df A data frame with columns `X` (longitude) and `Y` (latitude), one
#'   row per observation used in `model`. If `NULL` (default), the function
#'   looks for data saved in `model$data`; this is automatically available for
#'   objects returned by [basis_regression()] and [basis_regression_iv()] because
#'   they call `feols(..., data.save = TRUE)` internally.
#' @param near_neigh Number of nearest neighbours used to build the spatial
#'   weights matrix. Default `5`.
#'
#' @return An object of class `moran_test` with components:
#'   \describe{
#'     \item{`statistic`}{The observed Moran's I statistic.}
#'     \item{`p.value`}{One-sided p-value (H1: positive spatial autocorrelation).}
#'   }
#'
#' @details
#' The weights matrix is built from the `k` nearest neighbours of each
#' observation (row-standardised), matching the approach used internally by
#' [placebo()] and [synth()]. Duplicate coordinates are jittered by a small
#' random amount before neighbour search.
#'
#' Moran's I is computed via [spdep::moran.test()] on the (weighted) residuals
#' extracted directly from the `fixest` model object. If the model was estimated
#' with observation weights, the residuals are multiplied by the square root of
#' those weights before the test, matching the convention of
#' [spdep::lm.morantest()].
#'
#' When `model$data` is available (i.e. the model was estimated with
#' `data.save = TRUE`, as [basis_regression()] and [basis_regression_iv()] do
#' automatically), coordinates are inferred from `model$data` and no `df`
#' argument is needed.
#'
#' For `modelsummary` integration, `moran_i` and `moran_p` are added
#' automatically to the goodness-of-fit section when `spatInfer` is loaded and
#' the model has spatial data (i.e. `model$data` contains `X` and `Y`). Control
#' the display with `gof_map`:
#'
#' ```r
#' modelsummary::modelsummary(
#'   list(CK = ck),
#'   gof_map = list(
#'     list(raw = "nobs",     clean = "N",       fmt = 0),
#'     list(raw = "r.squared", clean = "R\u00b2", fmt = 2),
#'     list(raw = "moran_i",  clean = "Moran I", fmt = 3),
#'     list(raw = "moran_p",  clean = "Moran p", fmt = 3)
#'   )
#' )
#' ```
#'
#' For a plain [fixest::feols()] call (without `data.save = TRUE`), supply
#' coordinates explicitly and use the `gof_function` argument:
#'
#' ```r
#' modelsummary::modelsummary(
#'   list(OLS = mod),
#'   gof_function = function(m, ...) {
#'     mt <- moran_test(m, df = mydata)
#'     data.frame(moran_i = mt$statistic, moran_p = mt$p.value)
#'   }
#' )
#' ```
#'
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' set.seed(123)
#' opportunity <- opportunity |> dplyr::slice_sample(n = 250)
#'
#' ck <- basis_regression(mobility ~ single_mothers + gini, opportunity,
#'                        splines = 4, pc_num = 3, clusters = 5)
#'
#' # Standalone test (coordinates inferred from model$data)
#' moran_test(ck)
#'
#' # Passing coordinates explicitly
#' moran_test(ck, df = opportunity)
#'
#' # In modelsummary — moran_i / moran_p appear automatically via glance_custom
#' modelsummary::modelsummary(
#'   list(CK = ck),
#'   statistic  = c("conf.int", "p = {p.value}"),
#'   coef_omit  = "Intercept|PC",
#'   gof_map    = list(
#'     list(raw = "nobs",      clean = "N",       fmt = 0),
#'     list(raw = "r.squared", clean = "R2",       fmt = 2),
#'     list(raw = "moran_i",   clean = "Moran I",  fmt = 3),
#'     list(raw = "moran_p",   clean = "Moran p",  fmt = 3)
#'   )
#' )
moran_test <- function(model, df = NULL, near_neigh = 5) {

  model_data <- model$data  # NULL when data.save = FALSE

  # ── Resolve coordinates ───────────────────────────────────────────────────────
  if (!is.null(df)) {
    if (!all(c("X", "Y") %in% names(df)))
      stop("The data frame must contain columns named X (longitude) and Y (latitude).")
    Coords <- as.matrix(df[, c("X", "Y")])

  } else if (!is.null(model_data)) {
    if (!all(c("X", "Y") %in% names(model_data)))
      stop("model$data does not contain X and Y columns. Supply `df` explicitly.")
    # feols with data.save=TRUE stores all rows; use complete cases on formula
    # variables to recover exactly the rows used in estimation.
    fml_vars <- all.vars(formula(model))
    keep     <- stats::complete.cases(
      model_data[, intersect(fml_vars, names(model_data)), drop = FALSE]
    )
    Coords   <- as.matrix(model_data[keep, c("X", "Y")])

  } else {
    stop(
      "No data found. Either supply `df` with X and Y columns, ",
      "or re-estimate the model with `data.save = TRUE` (basis_regression ",
      "and basis_regression_iv do this automatically)."
    )
  }

  if (anyDuplicated(Coords) > 0) {
    set.seed(123)
    Coords <- Coords + matrix(stats::rnorm(2L * nrow(Coords), 0, 0.01), ncol = 2L)
  }

  # ── Spatial weights matrix ────────────────────────────────────────────────────
  nearest <- spdep::knn2nb(spdep::knearneigh(Coords, k = near_neigh, longlat = FALSE))
  listw   <- spdep::nb2listw(nearest, style = "W")

  # ── Residuals ─────────────────────────────────────────────────────────────────
  # Apply sqrt(weights) when the model was estimated with observation weights,
  # matching the weighted.residuals() convention used by lm.morantest().
  resids <- stats::residuals(model)
  wts    <- model$weights
  if (!is.null(wts) && length(wts) == length(resids))
    resids <- resids * sqrt(wts)

  result <- spdep::moran.test(resids, listw = listw)
  stat   <- unname(result$estimate["Moran I statistic"])
  pval   <- result$p.value

  structure(list(statistic = stat, p.value = pval), class = "moran_test")
}

#' @export
print.moran_test <- function(x, ...) {
  cat("Moran's I test on residuals\n")
  cat(sprintf("  I = %.4f,  p = %.4f\n", x$statistic, x$p.value))
  invisible(x)
}

#' Add Moran's I to modelsummary goodness-of-fit for fixest models
#'
#' Automatically called by [modelsummary::modelsummary()] when `spatInfer` is
#' loaded. Adds `moran_i` and `moran_p` rows to the GOF section for any
#' `fixest` model that has spatial coordinates saved (i.e. `model$data`
#' contains `X` and `Y`). Returns an empty data frame silently for other
#' `fixest` models.
#'
#' @param x A `fixest` object.
#' @param ... Ignored.
#' @return A one-row data frame with columns `moran_i` and `moran_p`, or an
#'   empty data frame if spatial data is unavailable.
#' @exportS3Method modelsummary::glance_custom
glance_custom.fixest <- function(x, ...) {
  dat <- x$data
  if (is.null(dat) || !all(c("X", "Y") %in% names(dat)))
    return(data.frame())

  mt <- tryCatch(moran_test(x), error = function(e) NULL)
  if (is.null(mt)) return(data.frame())
  data.frame(moran_i = mt$statistic, moran_p = mt$p.value)
}
