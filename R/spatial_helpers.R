
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

#####Regression weights: the variable of df named by weights, or 1 for every observation
#####when weights is NULL (unweighted).
regression_weights=function(df,weights){
  if(is.null(weights)) return(rep(1,nrow(df)))
  if(is.logical(weights))
    stop("weights takes the name of the weighting variable, for example weights=\"pop\", ",
         "not TRUE or FALSE. Leave it out for an unweighted regression.",call.=FALSE)
  if(!is.character(weights)||length(weights)!=1)
    stop("weights must be the name of one variable in your data.",call.=FALSE)
  if(!weights %in% names(df))
    stop("There is no variable called ",weights," in your data.",call.=FALSE)
  wts=df[[weights]]
  if(!is.numeric(wts))
    stop("The weighting variable ",weights," must be numeric.",call.=FALSE)
  if(anyNA(wts))
    stop("The weighting variable ",weights," cannot have missing values.",call.=FALSE)
  return(wts)
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

  df$wts=regression_weights(df,weights)

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
    hold_clus[,m]=make_clusters(Coords,m+1,k_medoids)
  }

  names(hold_clus)=paste0("clust_",2:(max_clus))
  return(hold_clus)
}

######Clusters used by the placebo and synth tests: 3 up to max_clus k-medoids clusters, as a list of factors.
sim_clusters=function(df,k_medoids,max_clus){
  hold_clus=generate_clusters(df,k_medoids,max_clus)
  return(as.list(hold_clus)[-1])
}

######Spatial trend of the outcome for synthetic outcome tests: a quadratic in longitude and
######latitude, with their product added when xy_interaction is TRUE.
synth_trend=function(df,xy_interaction){
  fm=if(xy_interaction) dep_var~poly(X,2)+poly(Y,2)+X:Y else dep_var~poly(X,2)+poly(Y,2)
  return(lm(fm,data=df))
}

######k clusters of the coordinates: k-medoids (pam), or clara for large datasets. The placebo
######tests and the final regressions both use this, so their clusters are the same. clara
######uses a fixed seed and leaves the user's random numbers unchanged.
make_clusters=function(Coords,k,k_medoids){
  if(k_medoids) return(factor(cluster::pam(Coords,k=k)$clustering))
  withr::local_preserve_seed()
  set.seed(321)
  return(factor(cluster::clara(Coords,k=k,samples=50)$clustering))
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


#############search for mle ests of matern params and generate simulated noise with these parameters.
#############For large datasets use exact_cholesky=FALSE which uses BRISC.
#############seeds are the seeds of the exact and BRISC simulations. trend=FALSE returns the noise
#############without the fitted values of lm_res added.
has_brisc=function() requireNamespace("BRISC",quietly=TRUE)

Noise_Sim=function(df,lm_res,nSim,exact_cholesky,Parallel,seeds=c(1234,123),trend=TRUE){
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
  if(Range>=max(rng_search))
    warning("The estimated range of spatial correlation is at the upper limit of the search, ",
            "the 95th percentile of distances between points, so the true range may be longer.",
            call.=FALSE)
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
    set.seed(seeds[1])
    Sim=KL%*%matrix(rnorm(nSim*nrow(Coords)),ncol=nSim)   #sims without original trend
  }else{    #Use BRISC approx Cholesky for large data
    set.seed(seeds[2])
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
  Sim=sd(Residuals)*Sim
  if(trend) Sim=lm_res$fitted+Sim         #Adding trend leaves results unchanged: resids orthogonal to spatial basis.
  return(list(Sim=Sim,matern_params=matern_params))
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
######Parallel workers use one BLAS thread each. A multithreaded BLAS such as OpenBLAS
######(the default on Ubuntu) otherwise starts a thread per core in every worker, and the
######oversubscribed workers run so slowly that the simulations appear to hang.
######progress, if given, is called with the number of runs finished after each one, or in parallel
######after each batch of runs: the runs are then split into about 20 batches.
run_sims=function(n,fun,Parallel,fixest_single_thread,progress=NULL){
  if(!Parallel){
    out=list()
    for (j in 1:n){
      out[[j]]=fun(j)
      if(!is.null(progress)) progress(1)
    }
    return(out)
  }

  n_cores=n_workers()
  old_blas=RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(1)    #inherited by forked workers
  on.exit(RhpcBLASctl::blas_set_num_threads(old_blas),add=TRUE)
  if(fixest_single_thread){
    old_threads=fixest::getFixest_nthreads()
    fixest::setFixest_nthreads(nthreads=1)
    on.exit(fixest::setFixest_nthreads(nthreads=old_threads),add=TRUE)
  }

  batch=if(is.null(progress)) n else n_cores*ceiling(n/(20*n_cores))
  batches=split(seq_len(n),ceiling(seq_len(n)/batch))
  run_batches=function(run){
    out=list()
    for(b in batches){
      out=c(out,run(b))
      if(!is.null(progress)) progress(length(b))
    }
    return(out)
  }

  if(!use_fork()){
    cl=parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterEvalQ(cl,RhpcBLASctl::blas_set_num_threads(1))
    if(fixest_single_thread) parallel::clusterEvalQ(cl,fixest::setFixest_nthreads(nthreads=1))
    #fun, which may hold large simulated data, is sent to the workers once, not with every batch
    parallel::clusterCall(cl,assign_sim_fun,fun)
    out=run_batches(function(b) parallel::parLapply(cl,b,call_sim_fun))
  }else{
    #mclapply only warns about failed workers; those are raised as errors below
    out=run_batches(function(b) suppressWarnings(parallel::mclapply(b,fun,mc.cores=n_cores)))
    failed=vapply(out,function(x) is.null(x)||inherits(x,"try-error"),logical(1))
    if(any(failed)){
      first=out[[which(failed)[1]]]
      if(inherits(first,"try-error")) stop(attr(first,"condition"))
      stop("A parallel worker returned no result, possibly from running out of memory. Try Parallel = FALSE.")
    }
  }
  return(out)
}

######Store the simulation function in a PSOCK worker, and call it. They are defined at the top
######level of the package, so sending them to a worker does not send the data of run_sims().
assign_sim_fun=function(fun){
  assign(".spatInfer_sim_fun",fun,envir=globalenv())
  return(NULL)
}
call_sim_fun=function(j) get(".spatInfer_sim_fun",envir=globalenv())(j)

######Progress bar for placebo_synth() and iv_placebo_synth(), counting the simulated regressions,
######with a status giving the step being run. It is created in the frame .envir of the calling
######function, and closes when that function returns. add(k) counts k more regressions.
sim_progress=function(total,.envir=parent.frame()){
  id=cli::cli_progress_bar(total=total,clear=FALSE,.envir=.envir,
                           format="{cli::pb_status} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}")
  done=0
  list(status=function(msg) cli::cli_progress_update(id=id,status=msg,set=done,force=TRUE,.envir=.envir),
       add=function(k){
         done<<-done+k
         cli::cli_progress_update(id=id,set=done,.envir=.envir)
       })
}

######Show the step being run on a progress bar, if there is one.
progress_status=function(progress,msg){
  if(!is.null(progress)) progress$status(msg)
}

######Inputs shared by the placebo and synthetic outcome tests: the data with the treatment and
######outcome renamed and the spatial basis added, the clusters for BCH and IM, and the SCPC spatial
######weights for each avc, which depend only on the locations.
sim_inputs=function(fm,df,splines,pc_num,weights,max_clus,k_medoids,avc,scpc_weighting){
  if(!is.null(avc)) avc=check_avc(avc)
  scpc_weighting=check_scpc_weighting(scpc_weighting)
  prep=prepare_spatial_data(fm,df,splines,pc_num,weights,max_clus)
  return(sim_inputs_finish(prep,splines,pc_num,max_clus,k_medoids,avc,scpc_weighting))
}

######Clusters from 3 up to max_clus and SCPC spatial weights added to the prepared data prep.
sim_inputs_finish=function(prep,splines,pc_num,max_clus,k_medoids,avc,scpc_weighting){
  prep$hold_clus=sim_clusters(prep$df,k_medoids,max_clus)
  prep$avc=avc
  prep$spcs=if(is.null(avc)) NULL else scpc_setup(prep$df,avc,scpc_weighting)
  prep$splines=splines
  prep$pc_num=pc_num
  return(prep)
}

######Confidence interval of the treatment coefficient, in the units of the coefficient.
coef_ci=function(estimate){
  ci=confint(estimate)
  term=treatment_term(estimate)
  return(c(ci[term,1],ci[term,2]))
}

######Name of the treatment in a fit: explan_var, or the placebo sim. Coefficients are found by
######name, not position, because feols drops collinear variables.
treatment_term=function(fit){
  term=intersect(c("sim","explan_var"),names(stats::coef(fit)))
  if(length(term)!=1)
    stop("The treatment variable was dropped from the regression because of collinearity.",call.=FALSE)
  return(term)
}

######One row of results: estimated p value, share of simulated p values below it and below 0.05,
######and the confidence interval. The interval is in the units of the coefficient, so it and
######its width are kept to 3 significant digits rather than rounded to a fixed number of decimals.
ci_summary_row=function(SE,Clusters,est_p,ci,sim_p_values){
  width_ci=signif(ci[2]-ci[1],3)
  conf_int=paste0("[",signif(ci[1],3),", ",signif(ci[2],3),"]")
  sim_p=mean(sim_p_values<est_p)
  sim_05=mean(sim_p_values<0.05)
  return(data.frame(SE=SE,Clusters=Clusters,est_p,sim_p,sim_05,width_ci,CI=conf_int))
}

######Show "." as the number of clusters of HC rows, and round. CI widths keep their significant digits.
finish_summary=function(rows){
  out=purrr::list_rbind(rows) |>
    dplyr::mutate(Clusters=ifelse(Clusters==0,".",as.character(Clusters))) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric)&!dplyr::any_of("width_ci"), \(x) round(x, 3)))
  rownames(out)=NULL
  return(out)
}

######p values of the treatment (or placebo) in simulated regression j with HC standard errors,
######with BCH and IM for each set of clusters in hold_clus, and with SCPC for each avc with spatial
######weights spcs. The regression is estimated once: its coefficients do not depend on the inference.
######A draw where SCPC fails gives NA SCPC p values. Names are hc, bch_<clusters>, im_<clusters>
######and scpc_<avc>.
ols_sim=function(j,Sim,eq_sim,df,hold_clus,avc=NULL,spcs=scpc_setup(df,avc)){
  df1=cbind.data.frame(sim=Sim[,j],df)
  fit=fixest::feols(eq_sim,data=df1,weights=~wts,vcov="hetero")
  term=treatment_term(fit)
  p=c(hc=unname(fixest::pvalue(fit)[term]))
  for(clust in hold_clus){
    G=nlevels(clust)
    p[paste0("bch_",G)]=fixest::pvalue(fit,cluster=clust)[[term]]
    p[paste0("im_",G)]=im_t_test(im_cluster_coefs(eq_sim,df1,clust))$p.value
  }
  if(length(avc)>0)
    p[paste0("scpc_",avc)]=tryCatch(scpc_stats(fit,df1,avc,spcs)[,"p"],
                                    error=function(e) rep(NA_real_,length(avc)))
  return(as.data.frame(as.list(p)))
}

######Simulated p values for every draw in Sim, run in parallel if asked.
######progress is passed to run_sims().
ols_sim_p_values=function(df,Sim,eq_sim,hold_clus,nSim,Parallel,avc=NULL,spcs=NULL,progress=NULL){
  out=run_sims(nSim,function(j) ols_sim(j,Sim,eq_sim,df,hold_clus,avc,spcs),Parallel,
               fixest_single_thread=TRUE,progress=progress)
  return(purrr::list_rbind(out))
}

######Results for the actual data: an HC row, BCH and IM rows for each set of clusters in hold_clus,
######and an SCPC row for each avc (held in Clusters), each with the estimated p value and confidence
######interval of the treatment, and the share of simulated p values below the estimated p value
######(sim_p) and below 0.05 (sim_05). sim_out holds the simulated p values from ols_sim().
ols_summary=function(df,eq_est,hold_clus,sim_out,avc=NULL,spcs=NULL){
  fit=fixest::feols(eq_est,data=df,weights=~wts,vcov="hetero")
  term=treatment_term(fit)
  hc=list(ci_summary_row("HC",0,fit$coeftable[term,4],coef_ci(fit),sim_out$hc))
  bch=im=list()
  for(clust in hold_clus){
    G=nlevels(clust)
    k=as.character(G)
    df1=df
    df1$clust=clust
    estimate=fixest::feols(eq_est,data=df1,weights=~wts,vcov=fixest::vcov_cluster(cluster=~clust))
    bch[[k]]=ci_summary_row("BCH",G,estimate$coeftable[term,4],coef_ci(estimate),sim_out[[paste0("bch_",G)]])
    im_coefs=im_cluster_coefs(eq_est,df,clust)
    n_missing=G-sum(is.finite(im_coefs$estimate))
    if(n_missing>0)
      warning("IM with ",G," clusters: the treatment coefficient could not be estimated in ",
              n_missing," cluster(s), which are left out of the t-test.",call.=FALSE)
    im_test=im_t_test(im_coefs)
    im[[k]]=ci_summary_row("IM",G,im_test$p.value,im_test$conf.int[1:2],sim_out[[paste0("im_",G)]])
  }
  scpc=list()
  if(length(avc)>0){
    est=scpc_stats(fit,df,avc,spcs)
    scpc=lapply(seq_along(avc),function(k){
      sim_p=scpc_sim_p(sim_out[[paste0("scpc_",avc[k])]],"SCPC",avc[k])
      ci_summary_row("SCPC",avc[k],unname(est[k,"p"]),unname(est[k,c("lower","upper")]),sim_p)
    })
  }
  return(finish_summary(c(hc,bch,im,scpc)))
}

#Moran z test for autocorr in residuals. Uses near_neigh nearest neighbours.
#Identical coordinates are fine: spdep warns about them, so those warnings are muffled.
moran=function(fm,df,near_neigh=5){
  Coords=as.matrix(df |> dplyr::select(X,Y))
  lm_1=lm(fm,df,weights=wts)
  nearest=withCallingHandlers(
    spdep::knearneigh(Coords,k=near_neigh,longlat = FALSE),
    warning=function(w){
      if(grepl("identical points",conditionMessage(w))) invokeRestart("muffleWarning")
    })
  nearest=spdep::knn2nb(nearest)   #k nearest neighbours for Moran
  nearest=spdep::nb2listw(nearest,style="W")
  moran=spdep::lm.morantest(lm_1,listw=nearest)$statistic[1,1]
  return(moran)
}

