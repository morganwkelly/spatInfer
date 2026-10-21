#' Regression table comparing BCH, IM and SCPC inference.
#'
#' @description
#' `regression_table()` makes a [modelsummary::modelsummary()] table with a
#' column for each of the BCH, IM and SCPC results supplied, giving coefficients
#' and 95% confidence intervals, optionally preceded by a column for the original
#' regression. Notes give the number of clusters used for BCH and IM and
#' the value of the SCPC average correlation bound `avc`.
#'
#' @param orig An optional [fixest::feols()] regression, typically the original
#'   regression without a spatial basis, shown in a first column labelled "Original"
#'   with the confidence intervals of the standard errors it was estimated with.
#' @param bch Output of [basis_regression()] with BCH standard errors.
#' @param im Output of [im_regression()].
#' @param scpc Output of [scpc_regression()].
#' @param statistic Passed to [modelsummary::modelsummary()]. Defaults to 95%
#'   confidence intervals; use for example `c("conf.int","p = {p.value}")` to add
#'   p values. Standard errors are not available, since IM gives only pseudo
#'   standard errors.
#' @param coef_omit Coefficients to leave out, as a regular expression. Defaults to
#'   the intercept and the principal components of the spatial basis.
#' @param gof_map Goodness of fit statistics shown. Defaults to the number of observations.
#' @param notes Further notes, added after the notes on clusters and `avc`.
#' @param ... Other arguments passed to [modelsummary::modelsummary()], such as
#'   `fmt`, `stars`, `title` or `output`.
#'
#' @details At least one of `orig`, `bch`, `im` and `scpc` must be given. The IM estimate
#'   is the mean of the coefficients estimated in each cluster, so it differs from
#'   the BCH and SCPC estimates, which come from the full sample. Use
#'   [regression_plot()] to plot the same results. Requires the modelsummary package.
#'
#' @return A modelsummary table, by default a tinytable object.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=250)
#' fm=mobility~single_mothers+short_commute+gini+dropout_rate+social_cap
#' bch=basis_regression(fm,opportunity,splines=4,pc_num=3,clusters=5)
#' im=im_regression(fm,opportunity,splines=4,pc_num=3,clusters=5)
#' sc=scpc_regression(fm,opportunity,splines=4,pc_num=3,avc=0.04)
#' orig=fixest::feols(fm,opportunity,vcov="hetero")
#' regression_table(orig=orig,bch=bch,im=im,scpc=sc,fmt=2)

regression_table=function(orig=NULL,bch=NULL,im=NULL,scpc=NULL,statistic="conf.int",
                          coef_omit="Intercept|^PC[0-9]+$",gof_map="nobs",notes=NULL,...){
  res=regression_models(orig,bch,im,scpc,"regression_table()")

  return(modelsummary::modelsummary(res$models,
                                    statistic=statistic,
                                    coef_omit=coef_omit,
                                    gof_map=gof_map,
                                    notes=c(list(paste(res$notes,collapse=" ")),as.list(notes)),
                                    ...))
}

#' Plot of regression coefficients comparing BCH, IM and SCPC inference.
#'
#' @description
#' `regression_plot()` makes a [modelsummary::modelplot()] plot of the
#' coefficients and 95% confidence intervals of each of the BCH, IM and SCPC
#' results supplied, optionally with the original regression. It is the plot
#' analogue of [regression_table()].
#'
#' @inheritParams regression_table
#' @param ... Other arguments passed to [modelsummary::modelplot()], such as
#'   `coef_map`, `coef_rename`, `facet` or `draw`.
#'
#' @details At least one of `orig`, `bch`, `im` and `scpc` must be given. The
#'   caption gives the number of clusters used for BCH and IM and the value of
#'   the SCPC average correlation bound `avc`. Requires the modelsummary package.
#'
#' @return A ggplot object, which can be modified further with ggplot2.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=250)
#' fm=mobility~single_mothers+short_commute+gini+dropout_rate+social_cap
#' bch=basis_regression(fm,opportunity,splines=4,pc_num=3,clusters=5)
#' im=im_regression(fm,opportunity,splines=4,pc_num=3,clusters=5)
#' sc=scpc_regression(fm,opportunity,splines=4,pc_num=3,avc=0.04)
#' orig=fixest::feols(fm,opportunity,vcov="hetero")
#' regression_plot(orig=orig,bch=bch,im=im,scpc=sc)

regression_plot=function(orig=NULL,bch=NULL,im=NULL,scpc=NULL,
                         coef_omit="Intercept|^PC[0-9]+$",...){
  res=regression_models(orig,bch,im,scpc,"regression_plot()")

  p=modelsummary::modelplot(res$models,coef_omit=coef_omit,...)
  if(inherits(p,"ggplot"))
    p=p+ggplot2::geom_vline(xintercept=0,linetype="dashed",colour="grey50")+
      ggplot2::labs(caption=paste(res$notes,collapse=" "))
  return(p)
}

######Check the regressions supplied and put them in a named list modelsummary accepts,
######with notes on clusters and avc.
regression_models=function(orig,bch,im,scpc,caller){
  if(!requireNamespace("modelsummary",quietly=TRUE))
    stop(caller," needs the modelsummary package. Install it with install.packages(\"modelsummary\").",
         call.=FALSE)
  if(is.null(orig)&&is.null(bch)&&is.null(im)&&is.null(scpc))
    stop("Give at least one of orig, bch, im and scpc.",call.=FALSE)
  if(inherits(bch,"iv_regression")||inherits(im,"iv_regression")||inherits(scpc,"iv_regression"))
    stop("For iv_regression(), iv_im_regression() and iv_scpc_regression() results use iv_regression_table() and iv_regression_plot().",
         call.=FALSE)

  models=list()
  model_notes=character(0)
  if(!is.null(orig)){
    if(!inherits(orig,"fixest"))
      stop("orig must be a fixest::feols() regression.",call.=FALSE)
    models$Original=as_modelsummary_list(broom::tidy(orig,conf.int=TRUE),stats::nobs(orig))
  }
  if(!is.null(bch)){
    if(!inherits(bch,"fixest"))
      stop("bch must be the output of basis_regression().",call.=FALSE)
    G=attr(bch$cov.scaled,"G")
    models$BCH=as_modelsummary_list(broom::tidy(bch,conf.int=TRUE),stats::nobs(bch))
    model_notes=c(model_notes,if(is.null(G)) "BCH: heteroskedasticity robust standard errors."
                  else paste0("BCH: ",G[1]," clusters."))
  }
  if(!is.null(im)){
    if(!inherits(im,"im_regression"))
      stop("im must be the output of im_regression().",call.=FALSE)
    models$IM=as_modelsummary_list(im$Coefs,im$N)
    model_notes=c(model_notes,paste0("IM: ",im$Clusters," clusters."))
  }
  if(!is.null(scpc)){
    if(!inherits(scpc,"scpc_regression"))
      stop("scpc must be the output of scpc_regression().",call.=FALSE)
    models$SCPC=as_modelsummary_list(scpc$Coefs,scpc$N)
    model_notes=c(model_notes,paste0("SCPC: avc = ",scpc$avc,"."))
  }
  return(list(models=models,notes=model_notes))
}

######Results in the form modelsummary accepts for a model it has no tidier for.
######Standard errors are left out, since IM gives only pseudo standard errors.
as_modelsummary_list=function(coefs,nobs){
  tidy=coefs |>
    dplyr::select(term,estimate,statistic,p.value,conf.low,conf.high)
  out=list(tidy=as.data.frame(tidy),glance=data.frame(nobs=nobs))
  class(out)="modelsummary_list"
  return(out)
}
