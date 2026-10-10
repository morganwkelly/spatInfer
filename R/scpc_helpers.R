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

######SCPC estimate, standard error, t, p value and 95% confidence interval of the
######treatment (explan_var, or the placebo sim) in a feols fit, one row per avc.
######The treatment is found by name, as feols drops collinear variables.
scpc_stats=function(fit,df,avc){
  treat=which(names(stats::coef(fit))%in%c("sim","explan_var"))
  if(length(treat)!=1)
    stop("The treatment variable was dropped from the regression because of collinearity.",call.=FALSE)
  out=t(vapply(avc,function(a){
    scpcR::scpc(fit,data=df,lon="X",lat="Y",avc=a,ncoef=treat)$scpcstats[treat,]
  },numeric(6)))
  colnames(out)=c("Coef","Std_Err","t","p","lower","upper")
  return(out)
}

######SCPC p values of the simulated regression j, one column per avc.
######A simulation where SCPC fails gives NA p values; these are left out by summary_scpc().
scpc_sim=function(j,Sim,eq_sim,df,avc){
  df1=cbind.data.frame(sim=Sim[,j],df)
  p=tryCatch({
    fit=fixest::feols(eq_sim,data=df1,weights=~wts)
    scpc_stats(fit,df1,avc)[,"p"]
  },error=function(e) rep(NA_real_,length(avc)))
  return(as.data.frame(as.list(stats::setNames(p,paste0("avc_",avc)))))
}

scpc_p_values=function(Sim,eq_sim,df,avc,nSim,Parallel){
  scpc_out=run_sims(nSim,function(j) scpc_sim(j,Sim,eq_sim,df,avc),Parallel,fixest_single_thread=TRUE)
  return(purrr::list_rbind(scpc_out))
}

######Results table: an HC row, then one SCPC row for each avc. Confidence intervals are
######normalized by the coefficient as in the other tables. Pseudo SE is one quarter of the
######width of the 95% confidence interval, in the units of the coefficient.
summary_scpc=function(df,eq_est,avc,scpc_out,hc_out){
  n_failed=sum(!stats::complete.cases(scpc_out))
  if(n_failed>0)
    warning("SCPC could not be computed in ",n_failed," of ",nrow(scpc_out),
            " simulations, which are left out.",call.=FALSE)
  scpc_out=scpc_out[stats::complete.cases(scpc_out),,drop=FALSE]

  fit=fixest::feols(eq_est,data=df,weights=~wts)
  est=scpc_stats(fit,df,avc)
  coef=abs(fit$coeftable[2,1])
  sim_summ=list()
  for(k in seq_along(avc)){
    ci=est[k,c("lower","upper")]
    sim_summ[[k]]=ci_summary_row("SCPC",avc[k],est[k,"p"],ci/coef,scpc_out[,k])
    sim_summ[[k]]$pseudo_se=(ci[2]-ci[1])/4
  }

  hc=fixest::feols(eq_est,data=df,weights=~wts,vcov="hetero")
  hc_ci=confint(hc)[2,]
  hc_row=hc_summary_row(df,eq_est,hc_out)
  hc_row$pseudo_se=(hc_ci[[2]]-hc_ci[[1]])/4

  Results=purrr::list_rbind(c(list(hc_row),sim_summ)) |>
    dplyr::rename(avc=Clusters) |>
    dplyr::mutate(avc=ifelse(SE=="HC",".",as.character(avc)),
                  dplyr::across(c(est_p,sim_p,sim_05,width_ci),\(x) round(x,3)),
                  pseudo_se=signif(pseudo_se,3))
  rownames(Results)=NULL
  return(list(Results=Results,Estimates=data.frame(avc=avc,est,row.names=NULL)))
}
