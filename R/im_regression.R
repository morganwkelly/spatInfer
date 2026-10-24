#' Spatial basis regression with IM inference on every coefficient.
#'
#' @description
#' `im_regression()` runs the spatial basis regression of [bch_regression()]
#' with Ibragimov-Mueller (IM) inference on every coefficient instead of BCH
#' standard errors. The sample is divided into k-medoids clusters, the regression
#' is estimated separately in each cluster, and each coefficient is tested with a
#' t-test of its cluster estimates.
#'
#' @param fm  The equation formula. The variable of interest is the first one on the right hand side.
#' @param df  The name of the dataset use. Longitude and latitude must be named X and Y and have no missing values.
#' @param splines The dimension of the linear tensor used, from optimal_basis function.
#' @param pc_num The number of principal component to include, again from optimal_basis.
#' @param clusters Number of k-medoids clusters to use, based on the placebo test of [placebo()].
#' @param weights The name of the weighting variable in the dataset, as a
#'   string, for example `weights="pop"`. Leave it out for an unweighted
#'   regression.
#' @param k_medoids Use k-medoids clustering (PAM), the default. Set to `FALSE`
#'   to use Clara, as for large datasets: use the same setting as in the placebo
#'   test so that the clusters are the same.
#'
#' @details The IM estimate of a coefficient is the mean of its estimates in the
#'   clusters. p values and 95% confidence intervals come from a t-test of these
#'   with one degree of freedom fewer than the number of clusters. The reported
#'   standard error is a pseudo standard error, one quarter of the width of the
#'   95% confidence interval, and the t statistic is the estimate divided by this
#'   pseudo standard error, so they are comparable with BCH and SCPC results. Clusters where a coefficient cannot be estimated, for example
#'   because its variable does not vary within the cluster, are left out of its
#'   t-test with a warning. The clusters are the same as those [bch_regression()]
#'   uses for BCH standard errors with the same number of clusters. The principal
#'   components of the spatial basis are left out of the results.
#'
#' @return An object of class `im_regression`, a list containing
#'   * `Coefs`: a data frame with one row for each coefficient giving the IM
#'     estimate, pseudo standard error, pseudo t statistic, p value, 95%
#'     confidence interval and the number of clusters used.
#'   * `Cluster_Coefs`: a data frame giving the estimate of each coefficient in
#'     each cluster, one row per cluster and one column per coefficient.
#'   * `Cluster_N`: the number of observations in each cluster.
#'   * `Clusters`, `N`, `Splines` and `PCs`: the number of clusters and observations,
#'     and the spatial basis used.
#'
#'   Use [im_table()] to make tables of the results, and [regression_table()] to
#'   compare them with BCH and SCPC results.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=250)
#' # Use the number of splines and PCs indicated by optimal_basis and
#' # the number of clusters from the IM placebo test.
#' im=im_regression(mobility~single_mothers+short_commute+gini+dropout_rate+social_cap,
#' opportunity,splines=4,pc_num=3,clusters=5)
#' im_table(im)
#' im_table(im,type="clusters")

im_regression=function(fm,df,splines,pc_num,clusters,weights=NULL,k_medoids=TRUE){
  if(!is.numeric(clusters)||length(clusters)!=1||clusters<2||clusters!=round(clusters))
    stop("clusters must be a whole number of at least 2.",call.=FALSE)
  basis=basis_data(fm,df,splines,pc_num,weights)
  df1=basis$df

  Coords=as.matrix(df1 |> dplyr::select(X,Y))
  clust=make_clusters(Coords,clusters,k_medoids)   #same clusters as BCH in bch_regression()

  #coefficients estimated in each cluster: one column per cluster, NA where not estimated
  fits=lapply(split(df1,clust),function(d) fixest::feols(basis$eq_est,data=d,weights=~wts))
  full=fixest::feols(basis$eq_est,data=df1,weights=~wts)
  terms=setdiff(names(stats::coef(full)),names(basis$pc))
  cluster_coefs=vapply(fits,function(f) unname(stats::coef(f)[terms]),numeric(length(terms)))
  cluster_coefs=matrix(cluster_coefs,nrow=length(terms),dimnames=list(terms,levels(clust)))

  Coefs=purrr::list_rbind(lapply(terms,function(term) im_coef_row(term,cluster_coefs[term,])))

  Cluster_Coefs=data.frame(cluster=levels(clust),t(cluster_coefs),check.names=FALSE,row.names=NULL)
  obj=list(Coefs=Coefs,
           Cluster_Coefs=Cluster_Coefs,
           Cluster_N=c(table(clust)),
           Clusters=clusters,
           N=stats::nobs(full),
           Splines=splines,
           PCs=pc_num)
  class(obj)="im_regression"
  return(obj)
}

######IM t-test of one coefficient from its cluster estimates. Clusters where it was not
######estimated are left out, with a warning.
im_coef_row=function(term,estimates){
  ok=is.finite(estimates)
  if(sum(!ok)>0)
    warning("IM: the coefficient of ",term," could not be estimated in ",sum(!ok),
            " cluster(s), which are left out of its t-test.",call.=FALSE)
  if(sum(ok)<2){
    warning("IM: the coefficient of ",term," is estimated in fewer than 2 clusters, so it cannot be tested.",
            call.=FALSE)
    return(data.frame(term=term,estimate=NA_real_,std.error=NA_real_,statistic=NA_real_,
                      p.value=NA_real_,conf.low=NA_real_,conf.high=NA_real_,clusters=sum(ok)))
  }
  tt=stats::t.test(estimates[ok])
  pseudo_se=diff(tt$conf.int)/4   #one quarter of the width of the 95% confidence interval
  return(data.frame(term=term,
                    estimate=unname(tt$estimate),
                    std.error=pseudo_se,
                    statistic=unname(tt$estimate)/pseudo_se,
                    p.value=tt$p.value,
                    conf.low=tt$conf.int[1],
                    conf.high=tt$conf.int[2],
                    clusters=sum(ok)))
}

#' @export
print.im_regression=function(x,...){
  cat("IM inference with",x$Clusters,"clusters,",x$N,"observations.\n\n")
  print(x$Coefs,digits=3,row.names=FALSE)
  cat("\nCoefficients in each cluster:\n")
  print(x$Cluster_Coefs,digits=3,row.names=FALSE)
  invisible(x)
}

#' Tables of IM regression results.
#'
#' @description
#' `im_table()` makes a table of the results of [im_regression()]: either a
#' standard coefficient table or the estimate of each coefficient in each cluster.
#'
#' @param im An object created by [im_regression()] or [iv_im_regression()].
#' @param type `"coefs"` (the default) for a table giving the IM estimate,
#'   pseudo standard error, p value and 95% confidence interval of each
#'   coefficient. `"clusters"` for a table giving each coefficient's estimate in
#'   each cluster, one row per cluster, with the number of observations in each
#'   cluster in the last column.
#' @param caption An optional caption for the table.
#' @param digits Number of significant digits shown in the `"coefs"` table.
#'   Defaults to 3.
#' @param decimals Number of decimal places the cluster estimates are rounded
#'   to in the `"clusters"` table. Defaults to 2.
#'
#' @return A tinytable object.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=250)
#' im=im_regression(mobility~single_mothers+short_commute+gini+dropout_rate+social_cap,
#' opportunity,splines=4,pc_num=3,clusters=5)
#' im_table(im,caption="IM estimates.")
#' im_table(im,type="clusters",caption="Coefficients estimated in each cluster.")

im_table=function(im,type=c("coefs","clusters"),caption="",digits=3,decimals=2){
  if(!inherits(im,c("im_regression","iv_im_regression")))
    stop("im must be the output of im_regression() or iv_im_regression().",call.=FALSE)
  type=match.arg(type)
  basis_note=paste0("Splines=",im$Splines,", PCs=",im$PCs,", N=",im$N,".")
  if(inherits(im,"iv_im_regression"))
    basis_note=paste0("2SLS in each cluster. ",basis_note)

  if(type=="coefs"){
    tab=im$Coefs |>
      dplyr::mutate(CI=paste0("[",signif(conf.low,digits),", ",signif(conf.high,digits),"]")) |>
      dplyr::select(Variable=term,Estimate=estimate,`Pseudo SE`=std.error,
                    p=p.value,CI,Clusters=clusters) |>
      dplyr::mutate(dplyr::across(c(Estimate,`Pseudo SE`),\(x) signif(x,digits)),
                    p=round(p,3))
    notes=list(paste0("IM inference with ",im$Clusters," k-medoids clusters. Estimates are means of the
                      coefficients estimated in each cluster; p values and 95% confidence intervals come
                      from a t-test of these with one degree of freedom fewer than the number of clusters used.
                      Pseudo SE is one quarter of the width of the confidence interval."),
               basis_note)
  }else{
    tab=im$Cluster_Coefs |>
      dplyr::mutate(dplyr::across(-cluster,\(x) round(x,decimals)),
                    Observations=unname(im$Cluster_N[cluster]),
                    cluster=paste("Cluster",cluster)) |>
      dplyr::rename(Cluster=cluster)
    notes=list("Coefficient values in each k-medoids cluster.",basis_note)
  }

  out=tinytable::tt(tab,theme="striped",align=paste0("l",strrep("r",ncol(tab)-1)),
                    notes=notes,caption=caption) |>
    tinytable::style_tt(i=0,bold=TRUE) |>
    tinytable::style_tt(j=1,bold=TRUE)
  return(out)
}

######Data for a spatial basis regression, as in bch_regression(): the outcome is renamed
######dep_var, the regression weights are in wts, and the principal components of the
######spatial basis are added. Explanatory variables keep their names.
basis_data=function(fm,df,splines,pc_num,weights){
  if(is.null(df$X)|is.null(df$Y))
    stop("You must have longitude and latitude variables named X and Y")
  if(sum(is.na(df$X))>0|sum(is.na(df$Y))>0)
    stop("You cannot have missing values in longitude and latitude.")
  df$wts=regression_weights(df,weights)

  new_names=set_names(fm)
  df=df |> dplyr::rename(dep_var=new_names$gg[1])
  pc=prin_comp(df,splines,pc_num)
  df=cbind.data.frame(df,pc)
  eq_est=paste0("dep_var~",new_names$orig_name,new_names$rhs)
  eq_est=as.formula(paste(eq_est,paste(names(pc),collapse="+"),sep="+"),env=parent.frame())
  return(list(df=df,eq_est=eq_est,pc=pc))
}
