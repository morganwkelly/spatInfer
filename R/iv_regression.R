#' Spatial basis instrumental variable regressions with BCH standard errors.
#'@description
#'
#'  `iv_regression()` estimates an instrumental variable regression by 2SLS
#'  with the spatial basis added, using BCH (clustered) or HC standard errors,
#'  and gives the Anderson-Rubin (AR) test and confidence set of the treatment
#'  effect and the first stage F statistic. Choose the
#'  number of clusters with [iv_placebo()]: the clusters are the same k-medoids
#'  clusters it uses. For IM inference use [iv_im_regression()].
#'
#'@inheritParams iv_placebo
#'@param clusters Number of k-medoids clusters to use for BCH, based on the
#'  placebo test. Not used with HC.
#'@param cov `"BCH"` (the default) for clustered standard errors, or `"HC"`
#'  for heteroskedasticity robust standard errors.
#'@param k_medoids Use k-medoids clustering (PAM), the default. Set to `FALSE`
#'  to use Clara, as for large datasets: use the same setting as in
#'  [iv_placebo()] so that the clusters are the same.
#'@details The AR confidence set contains the treatment effects `b` where the
#'  instrument is not significant at 5% in a regression of
#'  `outcome - b*treatment`, and can be unbounded when the instrument is weak:
#'  see [iv_placebo()].
#'
#'@return An object of class `iv_regression`: a list containing `Coefs`, the
#'  coefficient table of the treatment (named `fit_` followed by its name, as in
#'  fixest) and the controls, with the spatial basis left out; `AR`, the AR p
#'  value and confidence set of the treatment effect, with `AR_Set` giving the
#'  set as numbers (an interval, or `c(-Inf, a, b, Inf)` for two rays); the
#'  first stage F statistic `F`; the inference used, number of clusters,
#'  observations, splines and principal components; and the [fixest::feols()] regression `Fit`. Show
#'  it with [iv_regression_table()] or [iv_regression_plot()].
#'@export
#'
#' @examples
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' fm=d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag
#' # Use the number of splines and PCs from optimal_basis() and clusters from iv_placebo().
#' iv_bch=iv_regression(fm,china,splines=6,pc_num=8,clusters=6,weights="timepwt48")
#' iv_bch

iv_regression=function(fm,df,splines,pc_num,clusters=NULL,weights=NULL,cov=c("BCH","HC"),k_medoids=TRUE){
  cov=match.arg(cov)
  return(iv_basis_regression(fm,df,splines,pc_num,clusters,weights,cov,k_medoids))
}

#' Spatial basis instrumental variable regressions with IM inference.
#'@description
#'
#'  `iv_im_regression()` runs the spatial basis instrumental variable
#'  regression with IM inference on every coefficient: each is estimated by
#'  2SLS separately in k-medoids clusters (the same clusters [iv_regression()]
#'  and [iv_placebo()] use) and tested with a t-test of its cluster estimates.
#'  No Anderson-Rubin (AR) test is given, since IM is not valid for it: use
#'  [iv_regression()] for AR. No first stage F statistic is given.
#'
#'@inheritParams iv_regression
#'@param clusters Number of k-medoids clusters, based on the placebo test.
#'@details Estimates are means of the coefficients estimated in each cluster.
#'  Standard errors are pseudo standard errors, one quarter of the width of the
#'  95% confidence interval, and t statistics are the estimates divided by
#'  these.
#'
#'@return An object of class `iv_im_regression`, which is also an
#'  `iv_regression`: as for [iv_regression()], but without `AR` and `AR_Set`, with `F` set to `NA`, and in place of `Fit` it
#'  contains `Cluster_Coefs`, with one row per cluster and one column per
#'  coefficient, and `Cluster_N`, the number of observations in each cluster.
#'  [im_table()] shows its coefficients or its estimates in each cluster.
#'@export
#'
#' @examples
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' fm=d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag
#' iv_im=iv_im_regression(fm,china,splines=6,pc_num=8,clusters=4,weights="timepwt48")
#' iv_im
#' im_table(iv_im,type="clusters")

iv_im_regression=function(fm,df,splines,pc_num,clusters,weights=NULL,k_medoids=TRUE){
  obj=iv_basis_regression(fm,df,splines,pc_num,clusters,weights,"IM",k_medoids)
  class(obj)=c("iv_im_regression","iv_regression")
  return(obj)
}

######2SLS with the spatial basis and BCH, HC or IM inference, for iv_regression() and iv_im_regression().
iv_basis_regression=function(fm,df,splines,pc_num,clusters,weights,cov,k_medoids=TRUE){
  if(cov!="HC"&&(!is.numeric(clusters)||length(clusters)!=1||clusters<2||clusters!=round(clusters)))
    stop("clusters must be a whole number of at least 2 for BCH and IM.",call.=FALSE)
  if(cov=="HC") clusters=NULL
  ivf=parse_iv_formula(fm)
  prep=prepare_iv_data(fm,df,splines,pc_num,weights,max_clus=NULL)
  df1=prep$df
  eqs=build_iv_formulas(prep$rhs,prep$pc)

  #the regression reported keeps the names of the treatment and instrument
  df_user=df1 |> dplyr::rename(dplyr::all_of(stats::setNames(c("explan_var","inst_var"),c(ivf$treatment,ivf$instrument))))
  eq_user=as.formula(paste0("dep_var~",paste(names(prep$pc),collapse="+"),prep$rhs,
                            "|",ivf$treatment,"~",ivf$instrument))
  full=fixest::feols(eq_user,data=df_user,weights=~wts,vcov="hetero")
  terms=setdiff(names(stats::coef(full)),names(prep$pc))

  clust=NULL
  if(!is.null(clusters)){
    Coords=as.matrix(df1 |> dplyr::select(X,Y))
    clust=make_clusters(Coords,clusters,k_medoids)   #same clusters as iv_placebo()
  }

  obj=list()
  if(cov=="IM"){
    cluster_coefs=vapply(terms,function(term) iv_cluster_coefs(eq_user,df_user,clust,term),
                         numeric(nlevels(clust)))
    cluster_coefs=matrix(cluster_coefs,nrow=nlevels(clust),dimnames=list(levels(clust),terms))
    obj$Coefs=purrr::list_rbind(lapply(terms,function(term) im_coef_row(term,cluster_coefs[,term])))
    first_F=NA_real_   #no first stage F or AR test for IM
  }else{
    fit=if(cov=="BCH") summary(full,cluster=clust) else full
    obj$Coefs=broom::tidy(fit,conf.int=TRUE) |>
      dplyr::filter(term %in% terms) |>
      as.data.frame()
    fits=ar_fits(df1,eqs)
    ar_p=if(cov=="BCH") fixest::pvalue(fits$y,cluster=clust)[["inst_var"]] else fixest::pvalue(fits$y)[["inst_var"]]
    ar_set=ar_set_fixest(fits,clust)
    first_F=first_stage_F(fits$x,clust)
    obj$AR=data.frame(term=paste0("fit_",ivf$treatment),p.value=ar_p,CI=format_conf_set(ar_set))
    obj$AR_Set=ar_set
  }
  obj$F=first_F
  obj$Inference=cov
  obj$Clusters=clusters
  obj$N=stats::nobs(full)
  obj$Splines=splines
  obj$PCs=pc_num
  if(cov=="IM"){
    obj$Cluster_Coefs=data.frame(cluster=levels(clust),cluster_coefs,check.names=FALSE,row.names=NULL)
    obj$Cluster_N=c(table(clust))
  }else{
    obj$Fit=fit
  }
  class(obj)="iv_regression"
  return(obj)
}

######Description of the inference of an iv_regression object, for printing and table notes.
iv_inference_label=function(x){
  if(x$Inference=="HC") return("HC inference")
  paste(x$Inference,"inference with",x$Clusters,"clusters")
}

#' @export
print.iv_regression=function(x,...){
  cat("2SLS with",paste0(iv_inference_label(x),","),x$N,"observations.\n\n")
  print(x$Coefs,digits=3,row.names=FALSE)
  if(!is.null(x$AR))
    cat("\nAnderson-Rubin test of",x$AR$term,": p =",signif(x$AR$p.value,3),
        ", 95% confidence set",x$AR$CI,"\n")
  if(!is.na(x$F)) cat("First stage F (",iv_inference_label(x),") = ",signif(x$F,4),"\n",sep="")
  if(!is.null(x$Cluster_Coefs)){
    cat("\nCoefficients in each cluster:\n")
    print(x$Cluster_Coefs,digits=3,row.names=FALSE)
  }
  invisible(x)
}
