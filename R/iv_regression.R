#' Spatial basis instrumental variable regressions with BCH standard errors.
#'@description
#'
#'  `iv_regression()` estimates an instrumental variable regression by 2SLS
#'  with the spatial basis added, using BCH (clustered) or HC standard errors,
#'  and gives the Anderson-Rubin (AR) test and confidence set of the treatment
#'  effect, with the same inference, and the first stage F statistic. Choose the
#'  number of clusters with [iv_placebo()]: the clusters are the same k-medoids
#'  clusters it uses. For IM inference use [iv_im_regression()], and for SCPC
#'  [iv_scpc_regression()].
#'
#'@inheritParams iv_placebo
#'@param clusters Number of k-medoids clusters to use for BCH, based on the
#'  placebo test. Not used with HC.
#'@param cov `"BCH"` (the default) for clustered standard errors, or `"HC"`
#'  for heteroskedasticity robust standard errors.
#'@param k_medoids Use k-medoids clustering (PAM), the default. Set to `FALSE`
#'  to use Clara, as for large datasets: use the same setting as in
#'  [iv_placebo()] so that the clusters are the same.
#'@details The AR test and confidence set use the reduced form approach of
#'  Chernozhukov and Hansen (2008): the set contains the treatment effects `b`
#'  where the instrument is not significant at 5% in a regression of
#'  `outcome - b*treatment` on the instrument, controls and spatial basis, with
#'  BCH or HC inference. It is valid however weak the instrument, and can be
#'  unbounded: see [iv_placebo()].
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
#'  The Anderson-Rubin (AR) test and confidence set use IM inference too: the
#'  reduced forms of the outcome and treatment on the instrument are estimated
#'  in each cluster, and `b` is in the set when the t-test does not reject a
#'  zero mean of the cluster coefficients of `outcome - b*treatment` (see
#'  [iv_placebo()]). No first stage F statistic is given.
#'
#'@inheritParams iv_regression
#'@param clusters Number of k-medoids clusters, based on the placebo test.
#'@details Estimates are means of the 2SLS coefficients estimated in each
#'  cluster. Unlike these, the AR test needs no strong instrument within each
#'  cluster.
#'  Standard errors are pseudo standard errors, one quarter of the width of the
#'  95% confidence interval, and t statistics are the estimates divided by
#'  these.
#'
#'@return An object of class `iv_im_regression`, which is also an
#'  `iv_regression`: as for [iv_regression()], with IM inference for `AR` and `AR_Set`, with `F` set to `NA`, and in place of `Fit` it
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

#' Spatial basis instrumental variable regressions with SCPC inference.
#'@description
#'
#'  `iv_scpc_regression()` runs the spatial basis instrumental variable
#'  regression by 2SLS with SCPC inference (Mueller and Watson 2022, 2023) on
#'  every coefficient, and gives the Anderson-Rubin (AR) test and confidence set
#'  of the treatment effect with SCPC inference. Choose `avc` with the SCPC rows
#'  of [iv_placebo()].
#'
#'@inheritParams iv_regression
#'@param avc The SCPC upper bound on the average pairwise correlation, based on
#'  the placebo test of [iv_placebo()]. It must lie between 0.001 and 0.99.
#'@details 2SLS is least squares of the outcome on the instrumented treatment
#'  (its first stage fitted values) and the controls, with residuals from the
#'  structural equation. SCPC inference on the 2SLS coefficients treats it in
#'  this way, as [scpc_regression()] treats least squares, with the conditional
#'  critical value computed for the instrumented regressors. Like other 2SLS
#'  inference it is unreliable when the instrument is weak, unlike the AR test:
#'  for SCPC this is the test of the instrument in the regression of
#'  `outcome - b*treatment` on the instrument, controls and spatial basis,
#'  whose conditional critical value does not depend on `b` (see
#'  [iv_placebo()]). See [scpc_regression()] for samples of 4,500 or more
#'  observations and weighted regressions.
#'
#'@return An object of class `iv_scpc_regression`, which is also an
#'  `iv_regression`: as for [iv_regression()], with SCPC standard errors,
#'  t statistics, p values and confidence intervals in `Coefs`, SCPC inference
#'  for `AR` and `AR_Set`, `F` set to `NA`, and `avc` in place of `Clusters`.
#'  Show it with [iv_regression_table()] or [iv_regression_plot()].
#'@export
#'
#' @examples
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' fm=d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag
#' iv_sc=iv_scpc_regression(fm,china,splines=6,pc_num=8,avc=0.05,weights="timepwt48")
#' iv_sc

iv_scpc_regression=function(fm,df,splines,pc_num,avc,weights=NULL){
  if(length(avc)!=1)
    stop("avc must be a single number between 0.001 and 0.99.",call.=FALSE)
  avc=check_avc(avc)
  obj=iv_basis_regression(fm,df,splines,pc_num,NULL,weights,"SCPC",avc=avc)
  class(obj)=c("iv_scpc_regression","iv_regression")
  return(obj)
}

######2SLS with the spatial basis and BCH, HC, IM or SCPC inference, for iv_regression(),
######iv_im_regression() and iv_scpc_regression().
iv_basis_regression=function(fm,df,splines,pc_num,clusters,weights,cov,k_medoids=TRUE,avc=NULL){
  if(cov%in%c("BCH","IM")&&(!is.numeric(clusters)||length(clusters)!=1||clusters<2||clusters!=round(clusters)))
    stop("clusters must be a whole number of at least 2 for BCH and IM.",call.=FALSE)
  if(cov%in%c("HC","SCPC")) clusters=NULL
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
    ar=ar_im(eqs$eq_rf,df1,"inst_var",cbind(y=df1$dep_var,x=df1$explan_var),clust)
    if(ar$n_missing>0)
      warning("AR IM with ",clusters," clusters: the reduced form could not be estimated in ",
              ar$n_missing," cluster(s), which are left out of the t-test.",call.=FALSE)
    first_F=NA_real_   #no first stage F for IM
  }else if(cov=="SCPC"){
    fit=full
    spc=scpc_spatial(as.matrix(df1[fixest::obs(full),c("X","Y")]),avc)
    stats=scpc_coef_stats(full,spc,seq_along(stats::coef(full)))
    stats=stats[rownames(stats)%in%terms,,drop=FALSE]
    obj$Coefs=data.frame(term=rownames(stats),estimate=stats[,"Coef"],std.error=stats[,"Std_Err"],
                         statistic=stats[,"t"],p.value=stats[,"p"],conf.low=stats[,"lower"],
                         conf.high=stats[,"upper"],row.names=NULL)
    ar=ar_scpc(ar_data_rf(df1,eqs),spc)
    first_F=NA_real_   #no first stage F for SCPC
  }else{
    fit=if(cov=="BCH") summary(full,cluster=clust) else full
    obj$Coefs=broom::tidy(fit,conf.int=TRUE) |>
      dplyr::filter(term %in% terms) |>
      as.data.frame()
    rf=ar_data_rf(df1,eqs)
    ar=if(cov=="BCH") ar_bch(rf,clust) else ar_hc(rf)
    first_F=ar_F(ar)
  }
  conf_set=ar_set(ar)   #Anderson-Rubin test by the reduced form approach, with the same inference
  obj$AR=data.frame(term=paste0("fit_",ivf$treatment),p.value=ar$p,CI=format_conf_set(conf_set))
  obj$AR_Set=conf_set
  obj$F=first_F
  obj$Inference=cov
  obj$Clusters=clusters
  obj$avc=avc
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
  if(x$Inference=="SCPC") return(paste("SCPC inference with avc =",x$avc))
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
