######Instrumental variable regressions with one endogenous treatment and one instrument,
######written as in fixest: outcome ~ controls | treatment ~ instrument.

######Split an IV formula into outcome, controls, treatment and instrument.
######rhs holds the controls as "+c1+c2", or "" when there are none.
parse_iv_formula=function(fm){
  form_msg="The formula must have the form outcome ~ controls | treatment ~ instrument, or outcome ~ 1 | treatment ~ instrument without controls."
  fs=stringr::str_replace_all(deparse1(fm)," ","")
  parts=stringr::str_split_1(fs,"\\|")
  if(length(parts)!=2) stop(form_msg,call.=FALSE)
  main=stringr::str_split_1(parts[1],"~")
  iv=stringr::str_split_1(parts[2],"~")
  if(length(main)!=2|length(iv)!=2|any(c(main,iv)=="")) stop(form_msg,call.=FALSE)
  if(grepl("\\+",iv[1])|grepl("\\+",iv[2]))
    stop("Only one endogenous treatment and one instrument are supported.",call.=FALSE)
  controls=main[2]
  rhs=if(controls %in% c("1","0")) "" else paste0("+",controls)
  if(rhs!=""&&any(c(iv[1],iv[2]) %in% stringr::str_split_1(controls,"\\+")))
    stop("The treatment and instrument cannot also be controls.",call.=FALSE)
  return(list(outcome=main[1],treatment=iv[1],instrument=iv[2],rhs=rhs))
}

######The OLS formula outcome ~ treatment + controls of an IV formula.
iv_ols_formula=function(fm){
  ivf=parse_iv_formula(fm)
  return(as.formula(paste0(ivf$outcome,"~",ivf$treatment,ivf$rhs),env=environment(fm)))
}

######Shared setup for iv_placebo() and iv_synth(): check inputs, set weights, rename outcome,
######treatment and instrument to dep_var, explan_var and inst_var, and add the spatial basis.
prepare_iv_data=function(fm,df,splines,pc_num,weights,max_clus){
  ivf=parse_iv_formula(fm)
  model_vars=intersect(all.vars(fm),names(df))
  missing_vars=setdiff(c(ivf$outcome,ivf$treatment,ivf$instrument),names(df))
  if(length(missing_vars)>0)
    stop("Variable(s) not in your data: ",paste(missing_vars,collapse=", "),".",call.=FALSE)
  if(anyNA(df[model_vars]))
    stop("The variables in the formula cannot have missing values.",call.=FALSE)
  prep=prepare_spatial_data(iv_ols_formula(fm),df,splines,pc_num,weights,max_clus)
  prep$df=prep$df |> dplyr::rename(inst_var=dplyr::all_of(ivf$instrument))
  return(prep)
}

######Formulas with the spatial basis added. The 2SLS equation eq_est, the reduced form eq_rf of the
######outcome on the instrument (the Anderson-Rubin regression for a zero treatment effect),
######the first stage eq_fs, and their simulated versions. In a placebo the noise `sim` replaces the
######treatment and `sim_inst` the instrument; in a synthetic outcome test `sim` replaces the outcome.
######iv_term and rf_term name the coefficients tested in the simulated equations, and rf_outcome
######the outcome of the simulated reduced form.
build_iv_formulas=function(rhs,pc,sim_replaces=c("explan_var","dep_var")){
  sim_replaces=match.arg(sim_replaces)
  env=parent.frame()
  exog=paste0(paste(names(pc),collapse="+"),rhs)
  fm=function(x) as.formula(x,env=env)
  eqs=list(eq_est=fm(paste0("dep_var~",exog,"|explan_var~inst_var")),
           eq_rf=fm(paste0("dep_var~inst_var+",exog)),
           eq_fs=fm(paste0("explan_var~inst_var+",exog)))
  if(sim_replaces=="explan_var"){
    eqs$eq_sim=fm(paste0("dep_var~",exog,"|sim~sim_inst"))
    eqs$eq_sim_rf=fm(paste0("dep_var~sim_inst+",exog))
    eqs$iv_term="fit_sim"
    eqs$rf_term="sim_inst"
    eqs$rf_outcome="dep_var"
  }else{
    eqs$eq_sim=fm(paste0("sim~",exog,"|explan_var~inst_var"))
    eqs$eq_sim_rf=fm(paste0("sim~inst_var+",exog))
    eqs$iv_term="fit_explan_var"
    eqs$rf_term="inst_var"
    eqs$rf_outcome="sim"
  }
  return(eqs)
}

######Coefficient `term` estimated separately within each cluster. Clusters where it cannot be
######estimated give NA.
iv_cluster_coefs=function(eq,df,clust,term){
  vapply(split(df,clust),function(d){
    tryCatch(unname(coef(fixest::feols(eq,data=d,weights=~wts))[term]),
             error=function(e) NA_real_)
  },numeric(1),USE.NAMES=FALSE)
}

######p values of `term` with HC standard errors, with BCH and IM for each set of clusters in
######hold_clus, and with SCPC for each avc with spatial weights spcs. A draw where SCPC fails gives
######an NA p value. Names are hc, bch_<clusters>, im_<clusters> and scpc_<avc>.
iv_test_p=function(eq,df,term,hold_clus,avc=NULL,spcs=NULL){
  fit=fixest::feols(eq,data=df,weights=~wts,vcov="hetero")
  p=c(hc=unname(fixest::pvalue(fit)[term]))
  for(clust in hold_clus){
    G=nlevels(clust)
    p[paste0("bch_",G)]=fixest::pvalue(fit,cluster=clust)[term]
    p[paste0("im_",G)]=im_t_test(list(estimate=iv_cluster_coefs(eq,df,clust,term)))$p.value
  }
  j=match(term,names(stats::coef(fit)))
  for(k in seq_along(spcs))
    p[paste0("scpc_",avc[k])]=tryCatch(scpc_coef_stats(fit,spcs[[k]],j)[,"p"],error=function(e) NA_real_)
  return(p)
}

######p values of one simulated draw: the 2SLS (Wald) test of the treatment and the Anderson-Rubin
######test of no effect (computed by ar_sim_p()), each with HC, BCH, IM and (for each avc) SCPC
######inference.
iv_sim=function(df1,eqs,hold_clus,avc=NULL,spcs=NULL){
  wald=iv_test_p(eqs$eq_sim,df1,eqs$iv_term,hold_clus,avc,spcs)
  ar=ar_sim_p(df1,eqs,hold_clus,avc,spcs)
  out=c(stats::setNames(wald,paste0("wald_",names(wald))),stats::setNames(ar,paste0("ar_",names(ar))))
  return(as.data.frame(as.list(out)))
}

######Confidence set as text, and its width: infinite for an unbounded set.
format_conf_set=function(set){
  if(length(set)==4)
    return(paste0("(-Inf, ",signif(set[2],3),"] U [",signif(set[3],3),", Inf)"))
  if(all(is.infinite(set))) return("(-Inf, Inf)")
  return(paste0("[",signif(set[1],3),", ",signif(set[2],3),"]"))
}

conf_set_row=function(SE,Clusters,est_p,set,sim_p_values){
  if(length(set)==2&&all(is.finite(set)))
    return(ci_summary_row(SE,Clusters,est_p,set,sim_p_values))
  return(data.frame(SE=SE,Clusters=Clusters,est_p,sim_p=mean(sim_p_values<est_p),
                    sim_05=mean(sim_p_values<0.05),width_ci=Inf,CI=format_conf_set(set)))
}

######Warn when the coefficient could not be estimated in some clusters.
warn_missing_clusters=function(coefs,label,G){
  n_missing=sum(!is.finite(coefs))
  if(n_missing>0)
    warning(label," with ",G," clusters: the coefficient could not be estimated in ",
            n_missing," cluster(s), which are left out of the t-test.",call.=FALSE)
}

######Add the first stage F statistic to a row of results.
with_F=function(row,F){
  row$F=F
  return(row)
}

######Results for the actual data: estimated p values and confidence intervals of the 2SLS (Wald)
######test, and confidence sets of the Anderson-Rubin test, each with HC, BCH, IM and (for each avc)
######SCPC inference, with the share of simulated p values below each
######estimated p value (sim_p) and below 0.05 (sim_05). F is the first stage F statistic with the
######same inference and number of clusters as the row, and NA for IM and SCPC. For SCPC rows,
######Clusters holds avc. sim_out holds the simulated p values from iv_sim().
iv_summary=function(df,eqs,hold_clus,sim_out,avc=NULL,spcs=NULL){
  term="fit_explan_var"
  fit=fixest::feols(eqs$eq_est,data=df,weights=~wts,vcov="hetero")
  ar_stats=ar_data_stats(df,eqs,hold_clus,avc,spcs)
  F_of=function(se,G) ar_F(Filter(function(r) r$SE==se&&r$Clusters==G,ar_stats)[[1]]$st)
  wald=list(with_F(ci_summary_row("HC",0,unname(fixest::pvalue(fit)[term]),
                                  as.numeric(confint(fit)[term,]),sim_out$wald_hc),F_of("HC",0)))
  bch=im=list()
  for(clust in hold_clus){
    G=nlevels(clust)
    k=as.character(G)
    bch[[k]]=with_F(ci_summary_row("BCH",G,unname(fixest::pvalue(fit,cluster=clust)[term]),
                                   as.numeric(confint(fit,cluster=clust)[term,]),sim_out[[paste0("wald_bch_",G)]]),
                    F_of("BCH",G))
    im_coefs=iv_cluster_coefs(eqs$eq_est,df,clust,term)
    warn_missing_clusters(im_coefs,"IM",G)
    im_test=im_t_test(list(estimate=im_coefs))
    im[[k]]=with_F(ci_summary_row("IM",G,im_test$p.value,im_test$conf.int[1:2],sim_out[[paste0("wald_im_",G)]]),
                   NA_real_)   #no first stage F for IM
  }
  scpc=lapply(seq_along(spcs),function(k){
    est=scpc_coef_stats(fit,spcs[[k]],match(term,names(stats::coef(fit))))
    sim_p=scpc_sim_p(sim_out[[paste0("wald_scpc_",avc[k])]],"2SLS SCPC",avc[k])
    with_F(ci_summary_row("SCPC",avc[k],unname(est[,"p"]),unname(est[,c("lower","upper")]),sim_p),NA_real_)   #no first stage F for SCPC
  })
  ar=lapply(ar_stats,function(r){
    col=if(r$SE=="HC") "ar_hc" else paste0("ar_",tolower(r$SE),"_",r$Clusters)
    sim_p=sim_out[[col]]
    if(r$SE=="IM"&&r$st$n_missing>0)
      warning("AR IM with ",r$Clusters," clusters: the reduced form could not be estimated in ",
              r$st$n_missing," cluster(s), which are left out of the t-test.",call.=FALSE)
    if(r$SE=="SCPC") sim_p=scpc_sim_p(sim_p,"AR SCPC",r$Clusters)
    F=if(r$SE %in% c("HC","BCH")) ar_F(r$st) else NA_real_   #no first stage F for IM or SCPC
    with_F(conf_set_row(paste("AR",r$SE),r$Clusters,r$st$p,ar_set(r$st),sim_p),F)
  })
  out=purrr::list_rbind(c(wald,bch,im,scpc,ar)) |>
    dplyr::relocate("F",.after="Clusters") |>
    dplyr::mutate(Clusters=ifelse(Clusters==0,".",as.character(Clusters))) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric)&!dplyr::any_of("width_ci"), \(x) round(x, 3)))
  rownames(out)=NULL
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

######Simulated p values for every draw in Sim (and Sim_inst for placebos), run in parallel if asked.
######spcs are the SCPC spatial weights for each avc, computed once from the locations.
iv_sim_p_values=function(df,Sim,eqs,hold_clus,nSim,Parallel,Sim_inst=NULL,avc=NULL,spcs=NULL){
  sim_fun=function(j){
    df1=cbind.data.frame(sim=Sim[,j],df)
    if(!is.null(Sim_inst)) df1$sim_inst=Sim_inst[,j]
    iv_sim(df1,eqs,hold_clus,avc,spcs)
  }
  out=run_sims(nSim,sim_fun,Parallel,fixest_single_thread=TRUE)
  return(purrr::list_rbind(out))
}


######Clusters used by the IV tests: 3 up to max_clus k-medoids clusters, as a list of factors.
iv_clusters=function(df,k_medoids,max_clus){
  hold_clus=generate_clusters(df,k_medoids,max_clus)
  return(as.list(hold_clus)[-1])
}
