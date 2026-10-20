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
######iv_term and rf_term name the coefficients tested in the simulated equations.
build_iv_formulas=function(rhs,pc,sim_replaces=c("explan_var","dep_var")){
  sim_replaces=match.arg(sim_replaces)
  env=parent.frame()
  exog=paste0(paste(names(pc),collapse="+"),rhs)
  fm=function(x) as.formula(x,env=env)
  eqs=list(eq_est=fm(paste0("dep_var~",exog,"|explan_var~inst_var")),
           eq_rf=fm(paste0("dep_var~inst_var+",exog)),
           eq_fs=fm(paste0("explan_var~inst_var+",exog)),
           eq_diff=fm(paste0("ar_diff~inst_var+",exog)))
  if(sim_replaces=="explan_var"){
    eqs$eq_sim=fm(paste0("dep_var~",exog,"|sim~sim_inst"))
    eqs$eq_sim_rf=fm(paste0("dep_var~sim_inst+",exog))
    eqs$iv_term="fit_sim"
    eqs$rf_term="sim_inst"
  }else{
    eqs$eq_sim=fm(paste0("sim~",exog,"|explan_var~inst_var"))
    eqs$eq_sim_rf=fm(paste0("sim~inst_var+",exog))
    eqs$iv_term="fit_explan_var"
    eqs$rf_term="inst_var"
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

######p values of `term` with HC standard errors, and with BCH and (if im=TRUE) IM for each set of
######clusters in hold_clus. Names are hc, bch_<clusters> and im_<clusters>.
iv_test_p=function(eq,df,term,hold_clus,im=TRUE){
  fit=fixest::feols(eq,data=df,weights=~wts,vcov="hetero")
  p=c(hc=unname(fixest::pvalue(fit)[term]))
  for(clust in hold_clus){
    G=nlevels(clust)
    p[paste0("bch_",G)]=fixest::pvalue(fit,cluster=clust)[term]
    if(im) p[paste0("im_",G)]=im_t_test(list(estimate=iv_cluster_coefs(eq,df,clust,term)))$p.value
  }
  return(p)
}

######p values of one simulated draw: the 2SLS (Wald) test of the treatment and the
######Anderson-Rubin test of the instrument in the reduced form. AR has no IM version.
iv_sim=function(df1,eq_sim,eq_sim_rf,iv_term,rf_term,hold_clus){
  wald=iv_test_p(eq_sim,df1,iv_term,hold_clus)
  ar=iv_test_p(eq_sim_rf,df1,rf_term,hold_clus,im=FALSE)
  out=c(stats::setNames(wald,paste0("wald_",names(wald))),stats::setNames(ar,paste0("ar_",names(ar))))
  return(as.data.frame(as.list(out)))
}

######Anderson-Rubin confidence set for the treatment effect b: the values where the coefficient
######of the instrument in a regression of dep_var-b*explan_var is not significant at 5%.
######a and c are the reduced form and first stage coefficients of the instrument, Vyy, Vxx and Vxy
######their variances and covariance, and crit the squared critical value of the t test. The set is
######an interval c(lower,upper), possibly c(-Inf,Inf), or the union of two rays
######c(-Inf,lower,upper,Inf) when the first stage is not significant.
ar_conf_set=function(a,c,Vyy,Vxy,Vxx,crit){
  A=c^2-crit*Vxx
  B=-2*(a*c-crit*Vxy)
  C=a^2-crit*Vyy
  disc=B^2-4*A*C
  if(A>0){
    return((-B+c(-1,1)*sqrt(max(disc,0)))/(2*A))
  }
  if(disc<=0) return(c(-Inf,Inf))
  roots=sort((-B+c(-1,1)*sqrt(disc))/(2*A))
  return(c(-Inf,roots,Inf))
}

######Reduced form, first stage, and regression of dep_var-s*explan_var on the instrument. Their
######variances under any choice of standard errors give the covariance of the reduced form and
######first stage coefficients exactly, since the variance of dep_var-b*explan_var is quadratic in b.
ar_fits=function(df,eqs){
  s=sd(df$dep_var)/sd(df$explan_var)
  df$ar_diff=df$dep_var-s*df$explan_var
  fit=function(eq) fixest::feols(eq,data=df,weights=~wts,vcov="hetero")
  return(list(y=fit(eqs$eq_rf),x=fit(eqs$eq_fs),d=fit(eqs$eq_diff),s=s))
}

######AR confidence set with HC standard errors (clust=NULL) or BCH standard errors.
ar_set_fixest=function(fits,clust=NULL){
  se2=function(fit){
    se=if(is.null(clust)) fixest::se(fit) else fixest::se(fit,cluster=clust)
    return(unname(se["inst_var"])^2)
  }
  s=fits$s
  Vyy=se2(fits$y)
  Vxx=se2(fits$x)
  Vxy=(Vyy+s^2*Vxx-se2(fits$d))/(2*s)
  t_df=if(is.null(clust)) fixest::degrees_freedom(fits$y,"t") else nlevels(clust)-1
  return(ar_conf_set(unname(coef(fits$y)["inst_var"]),unname(coef(fits$x)["inst_var"]),
                     Vyy,Vxy,Vxx,qt(0.975,t_df)^2))
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

######First stage F statistic of the instrument, with the same inference as the tests it goes with.
######With one instrument it is the squared t statistic of the instrument in the first stage: with
######HC standard errors (clust=NULL) or BCH standard errors from the clusters in clust. No F is
######given for IM.
first_stage_F=function(fit_fs,clust=NULL){
  t=if(is.null(clust)) fixest::tstat(fit_fs) else fixest::tstat(fit_fs,cluster=clust)
  return(unname(t["inst_var"])^2)
}

######Add the first stage F statistic to a row of results.
with_F=function(row,F){
  row$F=F
  return(row)
}

######Results for the actual data: estimated p values and confidence intervals of the 2SLS (Wald)
######test with HC, BCH and IM inference and the Anderson-Rubin test with HC and BCH inference (IM is
######not valid for AR), with the share of simulated p values below each estimated p value (sim_p)
######and below 0.05 (sim_05). F is the first stage F statistic with the same inference and number of
######clusters as the row, and NA for IM. sim_out holds the simulated p values from iv_sim().
iv_summary=function(df,eqs,hold_clus,sim_out){
  term="fit_explan_var"
  fit=fixest::feols(eqs$eq_est,data=df,weights=~wts,vcov="hetero")
  fits=ar_fits(df,eqs)
  F_hc=first_stage_F(fits$x)
  wald=list(with_F(ci_summary_row("HC",0,unname(fixest::pvalue(fit)[term]),
                                  as.numeric(confint(fit)[term,]),sim_out$wald_hc),F_hc))
  ar=list(with_F(conf_set_row("AR HC",0,unname(fixest::pvalue(fits$y)["inst_var"]),
                              ar_set_fixest(fits),sim_out$ar_hc),F_hc))
  bch=im=ar_bch=list()
  for(clust in hold_clus){
    G=nlevels(clust)
    k=as.character(G)
    F_bch=first_stage_F(fits$x,clust)
    bch[[k]]=with_F(ci_summary_row("BCH",G,unname(fixest::pvalue(fit,cluster=clust)[term]),
                                   as.numeric(confint(fit,cluster=clust)[term,]),sim_out[[paste0("wald_bch_",G)]]),
                    F_bch)
    im_coefs=iv_cluster_coefs(eqs$eq_est,df,clust,term)
    warn_missing_clusters(im_coefs,"IM",G)
    im_test=im_t_test(list(estimate=im_coefs))
    im[[k]]=with_F(ci_summary_row("IM",G,im_test$p.value,im_test$conf.int[1:2],sim_out[[paste0("wald_im_",G)]]),
                   NA_real_)   #no first stage F for IM
    ar_bch[[k]]=with_F(conf_set_row("AR BCH",G,unname(fixest::pvalue(fits$y,cluster=clust)["inst_var"]),
                                    ar_set_fixest(fits,clust),sim_out[[paste0("ar_bch_",G)]]),F_bch)
  }
  out=purrr::list_rbind(c(wald,bch,im,ar,ar_bch)) |>
    dplyr::relocate("F",.after="Clusters") |>
    dplyr::mutate(Clusters=ifelse(Clusters==0,".",Clusters)) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric)&!dplyr::any_of("width_ci"), \(x) round(x, 3)))
  return(out)
}

######Simulated p values for every draw in Sim (and Sim_inst for placebos), run in parallel if asked.
iv_sim_p_values=function(df,Sim,eqs,hold_clus,nSim,Parallel,Sim_inst=NULL){
  sim_fun=function(j){
    df1=cbind.data.frame(sim=Sim[,j],df)
    if(!is.null(Sim_inst)) df1$sim_inst=Sim_inst[,j]
    iv_sim(df1,eqs$eq_sim,eqs$eq_sim_rf,eqs$iv_term,eqs$rf_term,hold_clus)
  }
  out=run_sims(nSim,sim_fun,Parallel,fixest_single_thread=TRUE)
  return(purrr::list_rbind(out))
}


######Clusters used by the IV tests: 3 up to max_clus k-medoids clusters, as a list of factors.
iv_clusters=function(df,k_medoids,max_clus){
  hold_clus=generate_clusters(df,k_medoids,max_clus)
  return(as.list(hold_clus)[-1])
}
