bch_sim=function(j,Sim,eq_sim,df,hold_clus){
  ##simulated regression j, estimated once: the coefficients do not depend on the clusters,
  ##so only its clustered p value is recomputed for each number of clusters
  df1=cbind.data.frame(sim=Sim[,j],df)
  fit=fixest::feols(eq_sim,data=df1,weights=~wts)
  term=treatment_term(fit)
  p=vapply(hold_clus,function(clust) fixest::pvalue(fit,cluster=clust)[[term]],numeric(1))
  return(as.data.frame(as.list(stats::setNames(p,paste0("clus_",seq_along(p)+1)))))
}



bch_p_values=function(Sim,eq_sim,df,hold_clus,nSim,max_clus,Parallel){
  bch_out=run_sims(nSim,function(j) bch_sim(j,Sim,eq_sim,df,hold_clus),Parallel,fixest_single_thread=TRUE)
  return(purrr::list_rbind(bch_out))
}

summary_bch=function(df,eq_est,hold_clus,max_clus,bch_out,hc_out){
  sim_summ=list()
  for(k in 1:(max_clus-1)){
    df1=df
    df1$clust=hold_clus[,k]
    estimate=fixest::feols(eq_est,data=df1,             #baseline estimate with real variables
                           weights = ~wts,
                           vcov=fixest::vcov_cluster(cluster=~clust))
    sim_summ[[k]]=ci_summary_row("BCH",k+1,estimate$coeftable[treatment_term(estimate),4],coef_ci(estimate),bch_out[,k])
  }
  sim_summ[[max_clus]]=hc_summary_row(df,eq_est,hc_out)
  return(finish_summary(sim_summ))
}
