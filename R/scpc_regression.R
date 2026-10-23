#' Spatial basis regression with SCPC inference on every coefficient.
#'
#' @description
#' `scpc_regression()` runs the spatial basis regression of [basis_regression()]
#' with SCPC inference (Mueller and Watson 2022, 2023) on every coefficient. The
#' SCPC code is ported from the scpcR package, and gives the same results for
#' weighted and unweighted regressions; see Details.
#'
#' @param fm  The equation formula. The variable of interest is the first one on the right hand side.
#' @param df  The name of the dataset use. Longitude and latitude must be named X and Y and have no missing values.
#' @param splines The dimension of the linear tensor used, from optimal_basis function.
#' @param pc_num The number of principal component to include, again from optimal_basis.
#' @param avc The SCPC upper bound on the average pairwise correlation, based on
#'   the placebo test of [placebo()]. It must lie between 0.001 and 0.99.
#' @param weights The name of the weighting variable in the dataset, as a
#'   string, for example `weights="pop"`. Leave it out for an unweighted
#'   regression.
#' @inheritParams placebo
#'
#' @details With 4,500 or more observations, SCPC uses scpcR's large sample approximation, which
#'   computes the spatial weights from random subsets of 1,000 locations and extends them
#'   to the rest. As in scpcR, this extension measures the distance between locations in
#'   degrees of longitude and latitude, whereas the rate at which correlation decays with
#'   distance is calibrated on great-circle distances in fractions of the earth's
#'   circumference. The weights therefore decay faster with distance than intended, and
#'   results can differ from the exact calculation used for smaller samples. The
#'   calculation is left unchanged so that results match scpcR.
#'
#'   With `weights`, conditional SCPC leaves the regressors unweighted by
#'   default, as scpcR does. With `scpc_weighting = "sqrt"` it is computed for
#'   the regression rescaled by the square root of the weights, the least squares
#'   regression equivalent to weighted least squares.
#'
#' @return An object of class `scpc_regression`, a list containing `Coefs`, a
#'   data frame with one row for each coefficient giving the estimate, SCPC
#'   standard error, t statistic, p value and 95% confidence interval, and `avc`,
#'   `N`, `Splines` and `PCs`. The principal components of the spatial basis are
#'   left out of `Coefs`. Use [regression_table()] to compare the results with
#'   BCH and IM results.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=250)
#' sc=scpc_regression(mobility~single_mothers+short_commute+gini+dropout_rate+social_cap,
#' opportunity,splines=4,pc_num=3,avc=0.04)
#' sc

scpc_regression=function(fm,df,splines,pc_num,avc,weights=NULL,scpc_weighting="scpcR"){
  if(length(avc)!=1)
    stop("avc must be a single number between 0.001 and 0.99.",call.=FALSE)
  avc=check_avc(avc)
  scpc_weighting=check_scpc_weighting(scpc_weighting)
  basis=basis_data(fm,df,splines,pc_num,weights)
  fit=fixest::feols(basis$eq_est,data=basis$df,weights=~wts)
  obs=fixest::obs(fit)
  stats=scpc_coef_stats(fit,scpc_spatial(as.matrix(basis$df[obs,c("X","Y")]),avc,
                                          cond_weights=scpc_weighting),
                        seq_along(stats::coef(fit)))
  stats=stats[!rownames(stats)%in%names(basis$pc),,drop=FALSE]
  Coefs=data.frame(term=rownames(stats),
                   estimate=stats[,1],
                   std.error=stats[,2],
                   statistic=stats[,3],
                   p.value=stats[,4],
                   conf.low=stats[,5],
                   conf.high=stats[,6],
                   row.names=NULL)
  obj=list(Coefs=Coefs,avc=avc,N=stats::nobs(fit),Splines=splines,PCs=pc_num)
  class(obj)="scpc_regression"
  return(obj)
}

#' @export
print.scpc_regression=function(x,...){
  cat("SCPC inference with avc =",x$avc,",",x$N,"observations.\n\n")
  print(x$Coefs,digits=3,row.names=FALSE)
  invisible(x)
}
