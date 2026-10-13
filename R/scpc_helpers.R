######SCPC inference (Mueller and Watson 2022, 2023) via scpcR::scpc().
######X and Y are longitude and latitude, so distances are great-circle distances (via geodist).

######avc values must lie strictly inside (0.001, 0.99), the range scpcR accepts.
check_avc=function(avc){
  if(!is.numeric(avc)||length(avc)==0||any(!is.finite(avc))||any(avc<=0.001|avc>=0.99))
    stop("avc must be a vector of numbers between 0.001 and 0.99.",call.=FALSE)
  if(anyDuplicated(avc)>0)
    stop("avc must not contain repeated values.",call.=FALSE)
  return(sort(avc))
}

######Most of the time of scpcR::scpc() goes on the spatial weights and their eigenvectors,
######which depend only on the locations and avc. scpc_setup() computes them once for each
######avc so that they can be reused in every simulated regression, and scpc_stats() then
######repeats the rest of scpcR::scpc() with scpcR's own functions, giving the same results.
scpc_internal=function(name) utils::getFromNamespace(name,"scpcR")

scpc_setup=function(df,avc){
  setOmsWfin=scpc_internal(".setOmsWfin")
  coords=as.matrix(df[,c("X","Y")])
  lapply(avc,function(a) setOmsWfin(coords,avc0=a,latlong=TRUE))
}

######SCPC estimate, standard error, t, p value and 95% confidence interval of the
######treatment (explan_var, or the placebo sim) in a feols fit, one row per avc.
######The treatment is found by name, as feols drops collinear variables. spcs is the
######output of scpc_setup() for the same observations and avc.
scpc_stats=function(fit,df,avc,spcs=scpc_setup(df,avc)){
  treat=which(names(stats::coef(fit))%in%c("sim","explan_var"))
  if(length(treat)!=1)
    stop("The treatment variable was dropped from the regression because of collinearity.",call.=FALSE)
  if(!identical(fixest::obs(fit),seq_len(nrow(df)))){
    ##Observations were dropped, so the saved weights do not apply.
    out=t(vapply(avc,function(a){
      scpcR::scpc(fit,data=df,lon="X",lat="Y",avc=a,ncoef=treat)$scpcstats[treat,]
    },numeric(6)))
  }else{
    out=t(vapply(spcs,function(spc) scpc_treat_stats(fit,spc,treat),numeric(6)))
  }
  colnames(out)=c("Coef","Std_Err","t","p","lower","upper")
  return(out)
}

######The steps of scpcR::scpc() after .setOmsWfin(), for a fixest regression without fixed
######effects, instruments or clusters, using conditional SCPC. With the exact method, used
######below 4500 observations, only the treatment is computed; the approximate method
######draws random numbers for each coefficient in turn, so all coefficients up to the
######treatment are computed, as scpcR::scpc() does.
scpc_treat_stats=function(fit,spc,treat){
  maxrp=scpc_internal(".maxrp")
  getcv=scpc_internal(".getcv")
  S=sandwich::estfun(fit)
  n=nrow(S)
  bread_inv=sandwich::bread(fit)/n
  coef_vec=stats::coef(fit)
  model_mat=scpc_internal(".get_scpc_model_matrix")(fit)
  cond=scpc_internal(".get_conditional_projection_setup")(fit,model_mat,n=n,uncond=FALSE)
  Wfin=spc$Wfin
  q=ncol(Wfin)-1
  perm=spc$perm
  random_state=spc$random_state
  mm=if(isTRUE(spc$large_n)) cond$model_mat[perm,,drop=FALSE] else cond$model_mat

  for(j in if(isTRUE(spc$large_n)) seq_len(treat) else treat){
    coef_j=coef_vec[[j]]
    wj=(as.numeric(n*bread_inv[j,]%*%t(S))+coef_j)[perm]
    denom=sqrt(sum((t(Wfin[,-1,drop=FALSE])%*%wj)^2))
    if(!is.finite(denom)||denom<=sqrt(.Machine$double.eps))
      stop("SCPC variance projection is degenerate.",call.=FALSE)
    tau_u=as.numeric(sqrt(q)*crossprod(Wfin[,1],wj)/denom)
    SE=as.numeric(denom/(sqrt(q)*sqrt(n)))
    p_u=maxrp(spc$Omsfin,q,abs(tau_u)/sqrt(q))$max

    xj=as.numeric(n*bread_inv[j,]%*%t(cond$model_mat))[perm]
    Wx=scpc_internal(".orthogonalize_W")(Wfin,xj,sign(xj),mm,
                                         include_intercept=cond$include_intercept,fixef_id=NULL)
    if(!all(is.finite(Wx)))
      stop("SCPC conditional projection produced non-finite values.",call.=FALSE)
    omsx=scpc_internal(".get_conditional_Oms")(spc,Wx,latlong=TRUE,random_state=random_state)
    random_state=omsx$random_state
  }
  p=max(p_u,maxrp(omsx$Oms,q,abs(tau_u)/sqrt(q))$max)
  cv=max(spc$cvfin,getcv(omsx$Oms,q,0.05))
  return(c(coef_j,SE,tau_u,p,coef_j-cv*SE,coef_j+cv*SE))
}

######SCPC p values of the simulated regression j, one column per avc.
######A simulation where SCPC fails gives NA p values; these are left out by summary_scpc().
scpc_sim=function(j,Sim,eq_sim,df,avc,spcs=scpc_setup(df,avc)){
  df1=cbind.data.frame(sim=Sim[,j],df)
  p=tryCatch({
    fit=fixest::feols(eq_sim,data=df1,weights=~wts)
    scpc_stats(fit,df1,avc,spcs)[,"p"]
  },error=function(e) rep(NA_real_,length(avc)))
  return(as.data.frame(as.list(stats::setNames(p,paste0("avc_",avc)))))
}

scpc_p_values=function(Sim,eq_sim,df,avc,nSim,Parallel){
  spcs=scpc_setup(df,avc)
  scpc_out=run_sims(nSim,function(j) scpc_sim(j,Sim,eq_sim,df,avc,spcs),Parallel,fixest_single_thread=TRUE)
  return(purrr::list_rbind(scpc_out))
}

######Results table: an HC row, then one SCPC row for each avc. Confidence intervals and
######Pseudo SE, one quarter of the width of the 95% confidence interval, are in the units
######of the coefficient.
summary_scpc=function(df,eq_est,avc,scpc_out,hc_out){
  n_failed=sum(!stats::complete.cases(scpc_out))
  if(n_failed>0)
    warning("SCPC could not be computed in ",n_failed," of ",nrow(scpc_out),
            " simulations, which are left out.",call.=FALSE)
  scpc_out=scpc_out[stats::complete.cases(scpc_out),,drop=FALSE]

  fit=fixest::feols(eq_est,data=df,weights=~wts)
  est=scpc_stats(fit,df,avc)
  sim_summ=list()
  for(k in seq_along(avc)){
    ci=est[k,c("lower","upper")]
    sim_summ[[k]]=ci_summary_row("SCPC",avc[k],est[k,"p"],ci,scpc_out[,k])
    sim_summ[[k]]$pseudo_se=(ci[2]-ci[1])/4
  }

  hc=fixest::feols(eq_est,data=df,weights=~wts,vcov="hetero")
  hc_ci=confint(hc)[2,]
  hc_row=hc_summary_row(df,eq_est,hc_out)
  hc_row$pseudo_se=(hc_ci[[2]]-hc_ci[[1]])/4

  Results=purrr::list_rbind(c(list(hc_row),sim_summ)) |>
    dplyr::rename(avc=Clusters) |>
    dplyr::mutate(avc=ifelse(SE=="HC",".",as.character(avc)),
                  dplyr::across(c(est_p,sim_p,sim_05),\(x) round(x,3)),
                  pseudo_se=signif(pseudo_se,3))
  rownames(Results)=NULL
  return(list(Results=Results,Estimates=data.frame(avc=avc,est,row.names=NULL)))
}
