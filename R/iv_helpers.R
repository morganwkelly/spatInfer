set_names_iv=function(fm){
  fm=deparse1(fm)
  fm=stringr::str_replace_all(fm," ","")
  gg=unlist(stringr::str_split_1(fm,"~"))
  var_dep=gg[1]
  var_inst=gg[3]
  rh_all=unlist(stringr::str_split(gg[2],"\\|"))
  rhs=rh_all[1]
  var_expl=rh_all[2]
  return(list(rhs=rhs,var_dep=var_dep,var_expl=var_expl,var_inst=var_inst,orig_name=var_expl))
}


#####noise instr variable. Placebo for explan var in first stage regression
bch_sim_iv1_placebo=function(j,Sim,eq_sim,df2){
  ##get simulated regression results with clustered standard errors
  sim=Sim[,j]
  df1=cbind.data.frame(sim,df2)
  rob_bch=fixest::feols(eq_sim,data=df1,
                        weights = ~wts,
                        vcov=fixest::vcov_cluster(cluster=~clust_bch))
  sim_res1=data.frame(
    rob_bch$iv_first_stage$explan_var$coeftable[2,4]
  )
  return(sim_res1)
}



bch_p_values_iv1_placebo=function(Sim,eq_sim,df,hold_clus,nSim,max_clus,Parallel){
  bch_out=list()
  for (l in 1:ncol(hold_clus)){
    df2=df
    df2$clust_bch=hold_clus[,l]
    clus_out=list()
    if(Parallel){
      n_cores=parallel::detectCores()-2  #number of cores to use
      fixest::setFixest_nthreads(nthreads=1)
      `%dopar%` <- foreach::`%dopar%`
      cl_k <- parallel::makeForkCluster(n_cores)
      doParallel::registerDoParallel(cl_k)
      clus_out=foreach::foreach(j=1:nSim) %dopar% {bch_sim_iv1_placebo(j,Sim,eq_sim,df2)}
      parallel::stopCluster(cl_k)
    }else{
      for (j in 1:nSim){
        clus_out[[j]]=bch_sim(j,Sim,eq_sim,df2)
      }
    }
    clus_out=purrr::list_rbind(clus_out)
    bch_out[[l]]=clus_out
    
  }
  names(bch_out)=paste0("clus_",2:(length(bch_out)+1))
  bch_out=purrr::list_cbind(bch_out)
  return(bch_out)
}


########synth outcome: noise explan var on instr in stage 1
bch_sim_iv1_synth=function(j,Sim,eq_sim,df2){
  ##get simulated regression results with clustered standard errors
  sim=Sim[,j]
  df1=cbind.data.frame(sim,df2)
  rob_bch=fixest::feols(eq_sim,data=df1,
                        weights = ~wts,
                        vcov=fixest::vcov_cluster(cluster=~clust_bch))
  sim_res1=data.frame(
    rob_bch$iv_first_stage$sim$coeftable[2,4]
  )
  return(sim_res1)
}



bch_p_values_iv1_synth=function(Sim,eq_sim,df,hold_clus,nSim,max_clus,Parallel){
  bch_out=list()
  for (l in 1:ncol(hold_clus)){
    df2=df
    df2$clust_bch=hold_clus[,l]
    clus_out=list()
    if(Parallel){
      n_cores=parallel::detectCores()-2  #number of cores to use
      fixest::setFixest_nthreads(nthreads=1)
      `%dopar%` <- foreach::`%dopar%`
      cl_k <- parallel::makeForkCluster(n_cores)
      doParallel::registerDoParallel(cl_k)
      clus_out=foreach::foreach(j=1:nSim) %dopar% {bch_sim_iv1_synth(j,Sim,eq_sim,df2)}
      parallel::stopCluster(cl_k)
    }else{
      for (j in 1:nSim){
        clus_out[[j]]=bch_sim(j,Sim,eq_sim,df2)
      }
    }
    clus_out=purrr::list_rbind(clus_out)
    bch_out[[l]]=clus_out
    
  }
  names(bch_out)=paste0("clus_",2:(length(bch_out)+1))
  bch_out=purrr::list_cbind(bch_out)
  return(bch_out)
}


hc_sim_iv1=function(j,Sim,eq_sim,df){
  ##get simulated regression results with hetero standard errors
  sim=Sim[,j]
  df1=cbind.data.frame(sim,df)
  rob_hc=fixest::feols(eq_sim,data=df1,
                       weights = ~wts,
                       vcov="hetero")
  
  sim_res2=data.frame(
    hc_p=rob_hc$coeftable[2,4]
  )
  return(sim_res2)
}

hc_p_values_iv1=function(Sim,eq_sim,df,nSim,Parallel){
  hc_out=list()
  if(Parallel){
    n_cores=parallel::detectCores()-2  #number of cores to use
    fixest::setFixest_nthreads(nthreads=1)
    `%dopar%` <- foreach::`%dopar%`
    cl_k <- parallel::makeForkCluster(n_cores)
    doParallel::registerDoParallel(cl_k)
    hc_out=foreach::foreach(j=1:nSim) %dopar% {hc_sim(j,Sim,eq_sim,df)}
    parallel::stopCluster(cl_k)
  }else{
    for (j in 1:nSim){
      hc_out[[j]]=hc_sim_iv1(j,Sim,eq_sim,df)
    }
  }
  hc_out=purrr::list_rbind(hc_out)
  
  return(hc_out)
}
