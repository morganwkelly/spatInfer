
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

#####Shared setup for the placebo and synth functions: check inputs, set weights,
#####rename outcome and treatment to dep_var and explan_var, and add the spatial basis principal components.
prepare_spatial_data=function(fm,df,splines,pc_num,weights,max_clus=NULL){
  if(is.null(df$X)|is.null(df$Y))
    stop("You must have longitude and latitude variables named X and Y")
  if(sum(is.na(df$X))>0|sum(is.na(df$Y))>0)
    stop("You cannot have missing values in longitude and latitude.")
  if(!is.null(max_clus)&&max_clus<3)
    stop("Your maximum number of clusters max_clus must be greater than 2.")

  if(!weights){
    df$wts=1
  }else{
    if(is.null(df$weights)){
      stop("There is no variable called weights in your data.")
    }else{
    df$wts=df$weights}
  }

#rename dependent and explanatory variables as dep_var and explan_var and list all other variables in a string called rhs
  new_names=set_names(fm)
  df=df |> dplyr::rename(dep_var=new_names$gg[1],
                   explan_var=new_names$orig_name)

#get the principal components that minimise BIC and add them to the dataset.
  pc=prin_comp(df,splines,pc_num)
  df=cbind.data.frame(df,pc)

  return(list(df=df,rhs=new_names$rhs,pc=pc))
}

#####Estimated and simulated equations with the spatial basis principal components added.
#####The simulated noise `sim` replaces either the treatment (placebo) or the outcome (synth).
#####Formulas keep the caller's environment, as when built inline with as.formula().
build_formulas=function(rhs,pc,sim_replaces=c("explan_var","dep_var")){
  sim_replaces=match.arg(sim_replaces)
  env=parent.frame()
  basis=paste(names(pc),collapse="+")
  eq_est=paste0("dep_var~explan_var",rhs)
  eq_sim=if(sim_replaces=="explan_var") paste0("dep_var~sim",rhs) else paste0("sim~explan_var",rhs)
  return(list(eq_est=as.formula(paste(eq_est,basis,sep="+"),env=env),
              eq_sim=as.formula(paste(eq_sim,basis,sep="+"),env=env)))
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

#####Get BIC minimizing set of principal components of tensor spline.
prin_comp=function(df,splines,pc_num){
gm_2=mgcv::bam(dep_var~
           te(X,Y,bs=c("bs"),
              k=splines,
              m=1),
         data=df,
         discrete=TRUE)

df_p=prcomp(model.matrix(gm_2))
pc=as.data.frame(df_p$x)
pc=pc[1:pc_num]   #keep number than minimize BIC

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

#############search for mle ests of matern params and generate simulated noise with these parameters.
#############For large datasets use exact_cholesky=FALSE which uses BRISC.
has_brisc=function() requireNamespace("BRISC",quietly=TRUE)

Noise_Sim=function(df,lm_res,nSim,exact_cholesky,Parallel){
  if(!exact_cholesky&&!has_brisc())
    stop("exact_cholesky = FALSE needs the BRISC package. Install it with install.packages(\"BRISC\").",
         call.=FALSE)
  Residuals=lm_res$residuals
  Coords=as.matrix(df |> dplyr::select(X,Y))
  rng_search=seq(0.025,1,by=0.025)*                 #search in increments of 0.025: proportions of
    quantile(fields::rdist(x1=Coords),probs=0.95)   #95th percentile distance bw point
  kriging_search=function(j){                       #find MLE of structure at each range
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
                                na.rm=TRUE,
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
  sim_krig=run_sims(length(rng_search),kriging_search,Parallel,fixest_single_thread=FALSE)
  sim_krig=purrr::list_rbind(sim_krig)
  matern_params=sim_krig|> dplyr::arrange(Likelihood) |> dplyr::slice(1)   #Choose MLE params
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
  }else{    #Use BRISC approx Cholesky for large data
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
    }
  Sim=lm_res$fitted+sd(Residuals)*Sim         #Adding trend leaves results unchanged: resids orthogonal to spatial basis.
  return(list(Sim=Sim,matern_params=matern_params))
}

######Calculate p_values for simulated data, using HC standard errors.
hc_p_values=function(Sim,eq_sim,df,nSim,Parallel){
  hc_out=run_sims(nSim,function(j) hc_sim(j,Sim,eq_sim,df),Parallel,fixest_single_thread=TRUE)
  hc_out=purrr::list_rbind(hc_out)

  return(hc_out)
}

######Number of parallel workers: all cores but two, and at least one.
######detectCores() can return NA, in which case one worker is used.
n_workers=function(){
  cores=parallel::detectCores()
  if(is.na(cores)) return(1L)
  return(as.integer(max(1,cores-2)))
}

######Forked workers (mclapply) are used only on Linux and other non-macOS Unix systems.
######On macOS, the Accelerate BLAS that R uses is not fork-safe: forked workers crash in
######linear algebra calls such as eigen() inside fields::Krig, especially under RStudio.
######Windows cannot fork. Both use separate worker processes (a PSOCK cluster) instead.
use_fork=function(){
  .Platform$OS.type=="unix"&&Sys.info()[["sysname"]]!="Darwin"
}

######Run fun(1), ..., fun(n) serially or in parallel and return the results as a list.
######fixest_single_thread=TRUE runs fixest on one thread in each worker; the user's
######thread setting is restored afterwards. No parallel backend is left registered.
run_sims=function(n,fun,Parallel,fixest_single_thread){
  if(!Parallel){
    out=list()
    for (j in 1:n){
      out[[j]]=fun(j)
    }
    return(out)
  }

  n_cores=n_workers()
  if(fixest_single_thread){
    old_threads=fixest::getFixest_nthreads()
    fixest::setFixest_nthreads(nthreads=1)
    on.exit(fixest::setFixest_nthreads(nthreads=old_threads),add=TRUE)
  }

  if(!use_fork()){
    cl=parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    if(fixest_single_thread) parallel::clusterEvalQ(cl,fixest::setFixest_nthreads(nthreads=1))
    out=parallel::parLapply(cl,1:n,fun)
  }else{
    #mclapply only warns about failed workers; those are raised as errors below
    out=suppressWarnings(parallel::mclapply(1:n,fun,mc.cores=n_cores))
    failed=vapply(out,function(x) is.null(x)||inherits(x,"try-error"),logical(1))
    if(any(failed)){
      first=out[[which(failed)[1]]]
      if(inherits(first,"try-error")) stop(attr(first,"condition"))
      stop("A parallel worker returned no result, possibly from running out of memory. Try Parallel = FALSE.")
    }
  }
  return(out)
}

######p values of simulated regressions for each set of clusters in hold_clus.
######The clusters are stored in column clust_name, which sim_fun uses.
cluster_p_values=function(Sim,eq_sim,df,hold_clus,nSim,Parallel,sim_fun,clust_name){
  clus_p=list()
  for (l in 1:ncol(hold_clus)){
    df2=df
    df2[[clust_name]]=hold_clus[,l]
    clus_out=run_sims(nSim,function(j) sim_fun(j,Sim,eq_sim,df2),Parallel,fixest_single_thread=TRUE)
    clus_out=purrr::list_rbind(clus_out)
    clus_p[[l]]=clus_out
  }
  names(clus_p)=paste0("clus_",2:(length(clus_p)+1))
  clus_p=purrr::list_cbind(clus_p)
  return(clus_p)
}


######Summary rows shared by summary_bch() and summary_im().
######Confidence interval of the treatment coefficient, normalized by the coefficient.
normalized_ci=function(estimate){
  ci=confint(estimate)/abs(estimate$coeftable[2,1])     #normalize by coef
  return(c(ci[2,1],ci[2,2]))
}

######One row of results: estimated p value, share of simulated p values below it and below 0.05,
######and the normalized confidence interval.
ci_summary_row=function(SE,Clusters,est_p,ci,sim_p_values){
  width_ci=round(ci[2]-ci[1],2)
  conf_int=paste0("[",round(ci[1],2),", ",round(ci[2],2),"]")
  sim_p=mean(sim_p_values<est_p)
  sim_05=mean(sim_p_values<0.05)
  return(data.frame(SE=SE,Clusters=Clusters,est_p,sim_p,sim_05,width_ci,CI=conf_int))
}

######HC row: baseline estimate with real variables and heteroskedasticity robust standard errors.
hc_summary_row=function(df,eq_est,hc_out){
  estimate=fixest::feols(eq_est,data=df,
                         weights = ~wts,
                         vcov="hetero")
  return(ci_summary_row("HC",0,estimate$coeftable[2,4],normalized_ci(estimate),hc_out[,1]))
}

######Drop the two cluster row, put HC first, and round.
finish_summary=function(sim_summ){
  sim_summ=purrr::list_rbind(sim_summ) |>
    dplyr::filter(Clusters!=2) |>
    dplyr::arrange(Clusters) |>
    dplyr::mutate(Clusters=ifelse(Clusters==0,".",Clusters)) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 3)))
  return(sim_summ)
}

#Moran z test for autocorr in residuals. Uses near_neigh nearest neighbours.
#Identical coordinates are jittered by about 1km when jitter_coords=TRUE.
moran=function(fm,df,near_neigh=5,jitter_coords=TRUE){
  Coords=as.matrix(df |> dplyr::select(X,Y))
  lm_1=lm(fm,df,weights=wts)
  if(jitter_coords&&anyDuplicated(Coords)>0){
    set.seed(123)
    Coords=Coords+matrix(rnorm(2*nrow(Coords),0,0.01),ncol=2)    #jitter by 1km to remove potential duplication
  }
  nearest=spdep::knn2nb(spdep::knearneigh(Coords,k=near_neigh,longlat = FALSE))   #k nearest neighbours for Moran
  nearest=spdep::nb2listw(nearest,style="W")
  moran=spdep::lm.morantest(lm_1,listw=nearest)$statistic[1,1]
  return(moran)
}

