######SCPC inference (Mueller and Watson 2022, 2023), computed in scpc_core.R.
######X and Y are longitude and latitude, so distances are great-circle distances (via geodist).

######avc values must lie strictly inside (0.001, 0.99).
check_avc=function(avc){
  if(!is.numeric(avc)||length(avc)==0||any(!is.finite(avc))||any(avc<=0.001|avc>=0.99))
    stop("avc must be a vector of numbers between 0.001 and 0.99.",call.=FALSE)
  if(anyDuplicated(avc)>0)
    stop("avc must not contain repeated values.",call.=FALSE)
  return(sort(avc))
}

######The spatial weights depend only on the locations and avc, and take most of the time,
######so scpc_setup() computes them once for each avc, to be reused in every simulated regression.
scpc_setup=function(df,avc){
  coords=as.matrix(df[,c("X","Y")])
  lapply(avc,function(a) scpc_spatial(coords,a))
}

######SCPC estimate, standard error, t, p value and 95% confidence interval of the
######treatment (explan_var, or the placebo sim) in a feols fit, one row per avc.
######The treatment is found by name, as feols drops collinear variables. spcs is the
######output of scpc_setup() for the same observations and avc.
scpc_stats=function(fit,df,avc,spcs=scpc_setup(df,avc)){
  treat=match(treatment_term(fit),names(stats::coef(fit)))
  obs=fixest::obs(fit)
  if(!identical(obs,seq_len(nrow(df)))) spcs=scpc_setup(df[obs,],avc)   #observations were dropped
  out=t(vapply(spcs,function(spc) scpc_coef_stats(fit,spc,treat)[1,],numeric(6)))
  colnames(out)=c("Coef","Std_Err","t","p","lower","upper")
  return(out)
}

######Simulated SCPC p values with the draws where SCPC failed (NA) left out, and a warning.
scpc_sim_p=function(sim_p,label,avc){
  n_failed=sum(is.na(sim_p))
  if(n_failed>0)
    warning(label," could not be computed in ",n_failed," of ",length(sim_p),
            " simulations for avc = ",avc,", which are left out.",call.=FALSE)
  return(as.numeric(stats::na.omit(sim_p)))
}
