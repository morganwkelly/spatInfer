
set_names=function(fm){
  fm=deparse1(fm)
  fm=stringr::str_replace_all(fm," ","")
  gg=stringr::str_split_1(fm,"~")
  rhs=stringr::str_split(gg[2],"\\+",n=2)
  rhs=unlist(rhs)[2]
  if(is.na(rhs)){
    rhs=""}else{
      rhs=paste0("+",rhs)
    }
  orig_name=stringr::str_split_1(gg[2],"\\+")[1]
  return(list(rhs=rhs,gg=gg,orig_name=orig_name))
}

generate_clusters=function(df,k_medoids,max_clus){
  Coords=as.matrix(df |> dplyr::select(X,Y))
  hold_clus=matrix(NA,nrow=nrow(Coords),ncol=(max_clus-1))
  hold_clus=as.data.frame(hold_clus)
  for(m in 1:(max_clus-1)){
    if(k_medoids){
      hold_clus[,m]=factor(cluster::pam(Coords,k=m+1)$clustering ) #k-medoids clusters
    }else{
      set.seed(321)
      hold_clus[,m]=factor(cluster::clara(Coords,k=m+1,samples=50)$clustering ) #use clara for large datasets
     # hold_clus[,m]=factor(kmeans(Coords,m+1)$cluster)
    }
  }

  names(hold_clus)=paste0("clust_",2:(max_clus))
  return(hold_clus)
}

prin_comp=function(df,splines,pc_num){
gm_2=mgcv::bam(dep_var~
           te(X,Y,bs=c("bs"),
              k=splines,
              m=1),
         data=df,
         discrete=T)

df_p=prcomp(model.matrix(gm_2))
pc=as.data.frame(df_p$x)
pc=pc[,1:pc_num]   #keep number than minimize BIC
return(pc)
}


hc_sim=function(j,Sim,eq_sim,df){
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

#############search for mle ests of matern params
Noise_Sim=function(df,lm_res,nSim,exact_cholesky,Parallel){
  Residuals=lm_res$residuals
  Coords=as.matrix(df |> dplyr::select(X,Y))
  rng_search=seq(0.025,1,by=0.025)*quantile(fields::rdist(x1=Coords),probs=0.95)   #fraction of 95th distance bw point
  kriging_search=function(j){
    ##find range, structure
    ##rng_search is set of ranges to try
    ##range measured in degrees
    hold_search=data.frame(range_search=rng_search[j],lambda=NA,loglik=NA,converge=NA,eff_df=NA)
    fit_search = fields::Krig(  x = Coords,
                                Y = scale(as.vector(Residuals)),
                                Covariance = "Matern",
                                smoothness = 0.5,
                                aRange = rng_search[j],
                                m=1,
                                na.rm=T,
                                give.warnings = F
    )

    hold_search[1,2:4]=fit_search$lambda.est[6,c(1,5,6)]
    hold_search[1,5]=fit_search$eff.df
    cov_par=hold_search |> as.numeric()
    Range=cov_par[1]
    Effective_Range=2*Range/quantile(fields::rdist(x1=Coords),probs=0.95)
    results=data.frame(Effective_Range,Structure=1/(1+cov_par[2]),Range,
                       Noise_to_Signal=cov_par[2],Likelihood=cov_par[3],
                       df=cov_par[5]/nrow(Coords),sigma_2=fit_search$sigma.MLE,tau_2=fit_search$tauHat.MLE)

    return(results)

  }
  if(Parallel){
    n_cores=parallel::detectCores()-2  #number of cores to use
    cl_k <- parallel::makeForkCluster(n_cores)
    doParallel::registerDoParallel(cl_k)
    `%dopar%` <- foreach::`%dopar%`
    sim_krig=foreach::foreach(j=1:length(rng_search)) %dopar% {kriging_search(j)}
    parallel::stopCluster(cl_k)
  }else{
    sim_krig=list()
    environment(kriging_search)=environment()
    for (j in 1:length(rng_search)){
      sim_krig[[j]]=kriging_search(j)
    }
  }
  sim_krig=purrr::list_rbind(sim_krig)
  matern_params=sim_krig|> dplyr::arrange(Likelihood) |> dplyr::slice(1)
  #return(matern_params)
  Range=matern_params$Range
  Structure=matern_params$Structure
  Effective_Range=2*Range/quantile(fields::rdist(x1=Coords),probs=0.95)   #fraction of 95th distance bw points


  if(exact_cholesky){                  #choose method to produce simulations:exact or approx Cholesky decomp
    # ################# spatial correlation matrix
    KL=Structure*fields::Matern(fields::rdist(x1=Coords),
                                range=Range,
                                smoothness=0.5   #exponential falloff
    )+
      diag(nrow(Coords))*(1-Structure)
    KL=t(chol(KL))

    ###################Generate simulated variables.
    set.seed(1234)
    Sim=KL%*%matrix(rnorm(nSim*nrow(Coords)),ncol=nSim)   #sims without original trend
  }else{
    set.seed(123)
    beg_seed=round(1e4*runif(nSim))
    Sim=BRISC::BRISC_simulation(
      Coords,
      sim_number = nSim,
      seeds=beg_seed,
      tau.sq = 1-matern_params$Structure,
      sigma.sq = matern_params$Structure,
      phi=1/matern_params$Range)


    Sim=Sim$output.data
    #Sim=lm_res$fitted+sd(Residuals)*Sim
  }
  Sim=lm_res$fitted+sd(Residuals)*Sim         #Adding trend leaves results unchanged: resids orthogonal to spatial basis.
  return(list(Sim=Sim,matern_params=matern_params))
}


hc_p_values=function(Sim,eq_sim,df,nSim,Parallel){
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
      hc_out[[j]]=hc_sim(j,Sim,eq_sim,df)
    }
  }
  hc_out=purrr::list_rbind(hc_out)

  return(hc_out)
}



#Moran z test for autocorr in residuals. Uses five nearest neighbours
moran=function(fm,df){
  Coords=as.matrix(df |> dplyr::select(X,Y))
  lm_1=lm(fm,df,weights=wts)
  if(anyDuplicated(Coords)>0){
    set.seed(123)
    Coords=Coords+matrix(rnorm(2*nrow(Coords),0,0.01),ncol=2)    #jitter by 1km to remove potential duplication
  }
  nearest=spdep::knn2nb(spdep::knearneigh(Coords,k=5,longlat = F))   #k nearest neighbours for Moran
  nearest=spdep::nb2listw(nearest,style="W")
  moran=spdep::lm.morantest(lm_1,listw=nearest)$statistic[1,1]
  return(moran)
}

