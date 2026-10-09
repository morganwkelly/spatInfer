bch_sim=function(j,Sim,eq_sim,df2){
  ##get simulated regression results with clustered standard errors
  sim=Sim[,j]
  df1=cbind.data.frame(sim,df2)
  rob_bch=fixest::feols(eq_sim,data=df1,
                        weights = ~wts,
                        vcov=fixest::vcov_cluster(cluster=~clust_bch))
  sim_res1=data.frame(
    rob_bch$coeftable[2,4]
  )
  return(sim_res1)
}



bch_p_values=function(Sim,eq_sim,df,hold_clus,nSim,max_clus,Parallel){
  cluster_p_values(Sim,eq_sim,df,hold_clus,nSim,Parallel,sim_fun=bch_sim,clust_name="clust_bch")
}

summary_bch=function(df,eq_est,hold_clus,max_clus,bch_out,hc_out){
  sim_summ=list()
  for(k in 1:(max_clus-1)){
    df1=df
    df1$clust=hold_clus[,k]
    estimate=fixest::feols(eq_est,data=df1,             #baseline estimate with real variables
                           weights = ~wts,
                           vcov=fixest::vcov_cluster(cluster=~clust))
    est_p=estimate$coeftable[2,4]
    ci=confint(estimate)/abs(estimate$coeftable[2,1])     #normalize by coef
    width_ci=round(ci[2,2]-ci[2,1],2)
    conf_int=paste0("[",round(ci[2,1],2),", ",round(ci[2,2],2),"]")
    # vr=data.frame(res=df1$explan_var,cl=hold_clus[,k]) |>
    #   dplyr::group_by(cl) |>
    #   dplyr::summarize(vr=var(res),2)

    sim_p=mean(bch_out[,k]<est_p)
    sim_05=mean(bch_out[,k]<0.05)
    sim_summ[[k]]=data.frame(SE="BCH",Clusters=k+1,est_p,sim_p,sim_05,
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
