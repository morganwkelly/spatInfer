######Anderson-Rubin (AR) inference by the reduced form approach of Chernozhukov and Hansen (2008),
######"The reduced form: A simple approach to inference with weak instruments", Economics Letters
######100(1), 68-71, with HC, BCH, IM and SCPC inference.
######
######Model. With outcome y, treatment x, instrument z and exogenous regressors W (constant,
######controls and spatial basis),
######    y = b x + W d + u,        x = c z + W g + v,
######so the reduced forms of y and x on z and W have instrument coefficients a = b c and c.
######
######The reduced form approach. For a candidate effect b0, regress y - b0 x on z and W. Its
######instrument coefficient is a - b0 c, which is zero when b0 = b, since then y - b0 x = W d + u
######does not depend on z. So any valid test of a single coefficient in this regression is a valid
######test of b = b0, and the set of b0 it does not reject is a confidence set for b. The
######regressors z and W do not involve x, so the test is valid however weak the instrument: the
######instrument's strength affects only the power of the test. At b0 = 0 it is the test of the
######instrument in the reduced form of the outcome, which is the AR test of no effect.
######
######Least squares is linear in the dependent variable. With the instrument coefficients a^ and
######c^ and residuals u^_y and u^_x of the reduced forms of y and x, the regression of y - b0 x
######has instrument coefficient a^ - b0 c^ and residuals u^_y - b0 u^_x. Each test below is of the
######form
######    reject b0  when  (a^ - b0 c^)^2 > k (S_yy - 2 b0 S_xy + b0^2 S_xx),
######where S is a 2 x 2 matrix computed from the two reduced forms (the variance of (a^, c^) under
######that inference) and the squared critical value k does not depend on b0. So the confidence set
######is where the quadratic
######    (c^2 - k S_xx) b0^2 - 2 (a^ c^ - k S_xy) b0 + (a^2 - k S_yy)
######is not positive, found exactly by ar_conf_set(): a bounded interval when c^2/S_xx > k, that is
######when the first stage is significant with the same inference, and otherwise two rays or the
######whole line. Each observation's residual enters the instrument coefficient with weight
######h_i = [(X'DX)^-1 X_i w_i]_z, where X = (z, W), D the diagonal matrix of the regression
######weights w and [.]_z the instrument's element, so that a^ - a = sum_i h_i u_yi.
######
######(a) HC. S = n/(n-K) sum_i h_i^2 (u^_yi, u^_xi)'(u^_yi, u^_xi), with K regressors, and
######    k = t(n-K)^2 at 0.975, the heteroskedasticity robust (HC1) t test of fixest.
######(b) BCH. S = G/(G-1) (n-1)/(n-K) sum_g s_g s_g', with s_g = sum_{i in g} h_i (u^_yi, u^_xi)
######    over the observations of cluster g, and k = t(G-1)^2, the clustered t test of fixest with
######    G-1 degrees of freedom (Bester, Conley and Hansen 2011).
######(c) IM. Both reduced forms are estimated separately in each cluster, giving (a^_g, c^_g). The
######    instrument coefficient of y - b0 x in cluster g is a^_g - b0 c^_g, which has mean zero in
######    every cluster when b0 = b, so the Ibragimov-Mueller t test of a zero mean of these cluster
######    estimates is valid. The estimates are the means of (a^_g, c^_g), S is their sample variance
######    divided by G, and k = t(G-1)^2. Unlike IM for 2SLS, this needs no strong instrument within
######    each cluster.
######(d) SCPC. The SCPC t statistic of a coefficient is tau = sqrt(q) W_0'm / |W_1:q' m| with
######    m_i = a^ + n h_i u^_i, where W_0 = 1/sqrt(n) and W_1:q are the spatial weights of scpc_spatial(),
######    and W_0'm = sqrt(n) a^ because sum_i h_i u^_i = 0. For y - b0 x, m(b0) = m_y - b0 m_x, so
######    tau(b0)^2 = (a^ - b0 c^)^2 / (S_yy - 2 b0 S_xy + b0^2 S_xx) with S = M'M/(nq), M = W_1:q'(m_y, m_x).
######    The critical value is the larger of the unconditional one and the conditional one of
######    Mueller and Watson (2023), which depends on the locations, regressors (z, W) and weights but
######    not on the outcome, so k is the same for every b0. The p value of no effect is that of
######    scpc_coef_stats() for the instrument in the reduced form of y.

######Reduced forms of the columns of Y (the outcome, and for confidence sets the treatment) on the
######instrument `term` and the other regressors of eq, with weights wts. Only the right hand side of
######eq is used. Returns the instrument coefficients coef, the residuals U, the weights h with which
######each residual enters the instrument coefficient, the regressors X, weights w, number of
######observations n and regressors K, and the position j of the instrument in X.
ar_reduced_forms=function(eq,df,term,Y){
  fit=fixest::feols(eq,data=df,weights=~wts)
  coef_names=names(stats::coef(fit))
  j=match(term,coef_names)
  if(is.na(j))
    stop("The instrument was dropped from the reduced form because of collinearity.",call.=FALSE)
  X=stats::model.matrix(fit)[,coef_names,drop=FALSE]
  w=stats::weights(fit)
  if(is.null(w)) w=rep(1,nrow(X))
  Y=as.matrix(Y)[fixest::obs(fit),,drop=FALSE]
  P=solve(crossprod(X,w*X),t(w*X))   #coefficients are P %*% Y
  coefs=P%*%Y
  return(list(coef=coefs[j,],U=Y-X%*%coefs,h=P[j,],X=X,w=w,n=nrow(X),K=ncol(X),j=j))
}

######AR statistics with a t test: instrument coefficients coef, their variance V, the squared
######critical value crit, and the p value of a zero coefficient of the first column (no effect).
ar_t_stats=function(coef,V,df){
  return(list(coef=coef,V=V,crit=stats::qt(0.975,df)^2,
              p=2*stats::pt(-abs(coef[[1]])/sqrt(V[1,1]),df)))
}

######(a) HC: heteroskedasticity robust variance, as fixest computes it.
ar_hc=function(rf){
  S=rf$h*rf$U
  return(ar_t_stats(rf$coef,rf$n/(rf$n-rf$K)*crossprod(S),rf$n-rf$K))
}

######(b) BCH: clustered variance with G-1 degrees of freedom, as fixest computes it.
ar_bch=function(rf,clust){
  Sg=rowsum(rf$h*rf$U,clust)
  G=nrow(Sg)
  return(ar_t_stats(rf$coef,G/(G-1)*(rf$n-1)/(rf$n-rf$K)*crossprod(Sg),G-1))
}

######(c) IM: t test of the instrument coefficients estimated in each cluster. Clusters where the
######reduced forms cannot be estimated are left out; n_missing counts them.
ar_im=function(eq,df,term,Y,clust){
  Y=as.matrix(Y)
  coefs=vapply(split(seq_len(nrow(df)),clust),function(i){
    tryCatch(ar_reduced_forms(eq,df[i,,drop=FALSE],term,Y[i,,drop=FALSE])$coef,
             error=function(e) rep(NA_real_,ncol(Y)))
  },numeric(ncol(Y)))
  coefs=matrix(coefs,ncol=ncol(Y),byrow=TRUE,dimnames=list(NULL,colnames(Y)))
  ok=rowSums(!is.finite(coefs))==0
  G=sum(ok)
  if(G<2) stop("IM needs the reduced form estimated in at least two clusters.",call.=FALSE)
  out=ar_t_stats(colMeans(coefs[ok,,drop=FALSE]),stats::cov(coefs[ok,,drop=FALSE])/G,G-1)
  out$n_missing=nrow(coefs)-G
  return(out)
}

######(d) SCPC, with the spatial weights spc of scpc_spatial() for the same observations. The
######critical value, the slowest part, is only computed when conf_set=TRUE.
ar_scpc=function(rf,spc,conf_set=TRUE){
  Wfin=spc$Wfin
  q=ncol(Wfin)-1
  if(nrow(Wfin)!=rf$n) stop("SCPC weights and regression have different observations.",call.=FALSE)
  m=sweep(rf$n*rf$h*rf$U,2,rf$coef,`+`)[spc$perm,,drop=FALSE]
  M=crossprod(Wfin[,-1,drop=FALSE],m)
  V=crossprod(M)/(q*rf$n)
  if(!is.finite(V[1,1])||V[1,1]<=.Machine$double.eps*max(1,rf$coef[[1]]^2))
    stop("SCPC variance projection is degenerate.",call.=FALSE)
  Omsx=scpc_cond_Oms(rf$X,rf$w,spc,rf$j)[[1]]
  tau=abs(rf$coef[[1]])/sqrt(V[1,1])
  p=max(scpc_maxrp(spc$Omsfin,q,tau/sqrt(q))$max,scpc_maxrp(Omsx,q,tau/sqrt(q))$max)
  crit=if(conf_set) max(spc$cvfin,scpc_cv(Omsx,q,0.05))^2 else NA_real_
  return(list(coef=rf$coef,V=V,crit=crit,p=p))
}

######AR confidence set from the statistics of one inference method, for the reduced forms of the
######outcome and the treatment.
ar_set=function(st){
  return(ar_conf_set(st$coef[[1]],st$coef[[2]],st$V[1,1],st$V[1,2],st$V[2,2],st$crit))
}

######First stage F statistic with the same inference: the squared t statistic of the instrument
######in the first stage.
ar_F=function(st) unname(st$coef[[2]]^2/st$V[2,2])

######Anderson-Rubin confidence set for the treatment effect b: the values where the coefficient
######of the instrument in a regression of dep_var-b*explan_var is not significant at 5%.
######a and c are the reduced form and first stage coefficients of the instrument, Vyy, Vxx and Vxy
######their variances and covariance, and crit the squared critical value. The set is an interval
######c(lower,upper), possibly c(-Inf,Inf), or the union of two rays c(-Inf,lower,upper,Inf) when
######the first stage is not significant.
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

######Reduced forms of the outcome and treatment in the real data.
ar_data_rf=function(df,eqs){
  ar_reduced_forms(eqs$eq_rf,df,"inst_var",cbind(y=df$dep_var,x=df$explan_var))
}

######AR statistics for the real data: HC, then BCH and IM for each set of clusters in hold_clus,
######then SCPC for each avc with spatial weights spcs. Each element has the inference SE, the
######number of clusters (or avc), the statistics and, for IM, the clusters left out.
ar_data_stats=function(df,eqs,hold_clus,avc=NULL,spcs=NULL){
  Y=cbind(y=df$dep_var,x=df$explan_var)
  rf=ar_data_rf(df,eqs)
  out=list(list(SE="HC",Clusters=0,st=ar_hc(rf)))
  bch=lapply(hold_clus,function(clust) list(SE="BCH",Clusters=nlevels(clust),st=ar_bch(rf,clust)))
  im=lapply(hold_clus,function(clust) list(SE="IM",Clusters=nlevels(clust),
                                           st=ar_im(eqs$eq_rf,df,"inst_var",Y,clust)))
  scpc=lapply(seq_along(spcs),function(k) list(SE="SCPC",Clusters=avc[k],st=ar_scpc(rf,spcs[[k]])))
  return(c(out,bch,im,scpc))
}

######p values of the AR test of no effect in one simulated draw: the instrument rf_term in the
######reduced form of the simulated outcome, with HC, BCH and IM inference for each set of clusters,
######and SCPC for each avc. A draw where SCPC fails gives an NA p value. Names are hc, bch_<G>,
######im_<G> and scpc_<avc>.
ar_sim_p=function(df1,eqs,hold_clus,avc=NULL,spcs=NULL){
  y=df1[[eqs$rf_outcome]]
  rf=ar_reduced_forms(eqs$eq_sim_rf,df1,eqs$rf_term,y)
  p=c(hc=ar_hc(rf)$p)
  for(clust in hold_clus){
    G=nlevels(clust)
    p[paste0("bch_",G)]=ar_bch(rf,clust)$p
    p[paste0("im_",G)]=ar_im(eqs$eq_sim_rf,df1,eqs$rf_term,y,clust)$p
  }
  for(k in seq_along(spcs))
    p[paste0("scpc_",avc[k])]=tryCatch(ar_scpc(rf,spcs[[k]],conf_set=FALSE)$p,error=function(e) NA_real_)
  return(p)
}
