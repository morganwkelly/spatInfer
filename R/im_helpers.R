######Coefficient of the treatment (or placebo) estimated separately within each cluster.
im_cluster_coefs=function(eq,df,clust){
  split(df,clust) |>
    purrr::map(~fixest::feols(eq, data = .x,weights=~wts)) |>
    purrr::map_df(broom::tidy) |>
    dplyr::filter(term == 'sim'|term == 'explan_var') |>
    dplyr::select(estimate)
}

######t-test of the cluster coefficients. Clusters where the coefficient could not be
######estimated (missing or not finite) are left out.
im_t_test=function(im_coefs){
  estimate=im_coefs$estimate
  return(t.test(estimate[is.finite(estimate)]))
}
