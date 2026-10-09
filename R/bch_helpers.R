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
    sim_summ[[k]]=ci_summary_row("BCH",k+1,estimate$coeftable[2,4],normalized_ci(estimate),bch_out[,k])
  }
  sim_summ[[max_clus]]=hc_summary_row(df,eq_est,hc_out)
  return(finish_summary(sim_summ))
}
