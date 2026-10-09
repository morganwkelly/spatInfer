im_p_values=function(Sim,eq_sim,df,hold_clus,nSim,max_clus,Parallel){
  cluster_p_values(Sim,eq_sim,df,hold_clus,nSim,Parallel,sim_fun=im_sim,clust_name="clust_im")
}

im_sim=function(j,Sim,eq_sim,df2){
  ##get simulated regression results with IM
  sim=Sim[,j]
  df1=cbind.data.frame(sim,df2)
  im_coefs=im_cluster_coefs(eq_sim,df1,df1$clust_im)
  im=t.test(im_coefs,na.action=na.fail())
  sim_res1=data.frame(
    im$p.value
  )
  return(sim_res1)
}

######Coefficient of the treatment (or placebo) estimated separately within each cluster.
im_cluster_coefs=function(eq,df,clust){
  split(df,clust) |>
    purrr::map(~fixest::feols(eq, data = .x,weights=~wts)) |>
    purrr::map_df(broom::tidy) |>
    dplyr::filter(term == 'sim'|term == 'explan_var') |>
    dplyr::select(estimate)
}

summary_im=function(df,eq_est,hold_clus,max_clus,im_out,hc_out){
  sim_summ=list()
  for(k in 1:(max_clus-1)){
    im=t.test(im_cluster_coefs(eq_est,df,hold_clus[,k]),na.action=na.fail())
    ci=im$conf.int[1:2]/abs(im$estimate)    #normalize by coef
    sim_summ[[k]]=ci_summary_row("IM",k+1,im$p.value,ci,im_out[,k])
  }
  sim_summ[[max_clus]]=hc_summary_row(df,eq_est,hc_out)
  return(finish_summary(sim_summ))
}


coef_im=function(df,eq_est,hold_clus,max_clus){
  
  estimate=fixest::feols(eq_est,data=df,             #baseline estimate with real variables
                         weights = ~wts,
                         vcov="hetero")
  coef_orig=estimate$coeftable[2,1]
  
  im_coef=list()
  for(k in 1:(max_clus-1)){
    im_coef[[k]]=im_cluster_coefs(eq_est,df,hold_clus[,k])/coef_orig
  }
  return(im_coef)
}
