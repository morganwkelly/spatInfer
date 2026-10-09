im_p_values=function(Sim,eq_sim,df,hold_clus,nSim,max_clus,Parallel){
  cluster_p_values(Sim,eq_sim,df,hold_clus,nSim,Parallel,sim_fun=im_sim,clust_name="clust_im")
}

im_sim=function(j,Sim,eq_sim,df2){
  ##get simulated regression results with IM
  sim=Sim[,j]
  df1=cbind.data.frame(sim,df2)
  im_coefs=split(df1,df1$clust_im) |>
    purrr::map(~fixest::feols(eq_sim, data = .x,weights=~wts)) |>
    purrr::map_df(broom::tidy) |>
    dplyr::filter(term == 'sim'|term == 'explan_var') |>
    dplyr::select(estimate)
  im=t.test(im_coefs,na.action=na.fail())
  sim_res1=data.frame(
    im$p.value
  )
  return(sim_res1)
}

summary_im=function(df,eq_est,hold_clus,max_clus,im_out,hc_out){
  sim_summ=list()
  for(k in 1:(max_clus-1)){
    df1=df
    df1$clust=hold_clus[,k]
    im=split(df1,df1$clust) |>
      purrr::map(~fixest::feols(eq_est, data = .x,weights=~wts)) |>
      purrr::map_df(broom::tidy) |>
      dplyr::filter(term == 'sim'|term == 'explan_var') |>
      dplyr::select(estimate)
    im=t.test(im,na.action=na.fail())
    est_p=im$p.value
    ci=im$conf.int[1:2]/abs(im$estimate)    #normalize by coef
    width_ci=round(ci[2]-ci[1],2)
    conf_int=paste0("[",round(ci[1],2),", ",round(ci[2],2),"]")

    sim_p=mean(im_out[,k]<est_p)
    sim_05=mean(im_out[,k]<0.05)
    sim_summ[[k]]=data.frame(SE="IM",Clusters=k+1,est_p,sim_p,sim_05,
                             width_ci,CI=conf_int)
  }
  ##HC
  estimate=fixest::feols(eq_est,data=df1,             #baseline estimate with real variables
                         weights = ~wts,
                         vcov="hetero")
  est_p=estimate$coeftable[2,4]
  ci=confint(estimate)/abs(estimate$coeftable[2,1])     #normalize by coef
  width_ci=round(ci[2,2]-ci[2,1],2)
  conf_int=paste0("[",round(ci[2,1],2),", ",round(ci[2,2],2),"]")
  sim_p=mean(hc_out[,1]<est_p)
  sim_05=mean(hc_out[,1]<0.05)
  sim_summ[[max_clus]]=data.frame(SE="HC",Clusters=0,est_p,sim_p,sim_05,width_ci,CI=conf_int)
  ###############

  sim_summ=purrr::list_rbind(sim_summ) |>
    dplyr::filter(Clusters!=2) |>
    dplyr::arrange(Clusters) |>
    dplyr::mutate(Clusters=ifelse(Clusters==0,".",Clusters)) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 3)))
  return(sim_summ)
}


coef_im=function(df,eq_est,hold_clus,max_clus){
  
  estimate=fixest::feols(eq_est,data=df,             #baseline estimate with real variables
                         weights = ~wts,
                         vcov="hetero")
  coef_orig=estimate$coeftable[2,1]
  
  im_coef=list()
  for(k in 1:(max_clus-1)){
    df1=df
    df1$clust=hold_clus[,k]
    im=split(df1,df1$clust) |>
      purrr::map(~fixest::feols(eq_est, data = .x,weights=~wts)) |>
      purrr::map_df(broom::tidy) |>
      dplyr::filter(term == 'sim'|term == 'explan_var') |>
      dplyr::select(estimate)
    im_coef[[k]]=im/coef_orig
  }
  return(im_coef)
}