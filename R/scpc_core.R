######Spatial correlation principal components (SCPC) inference of Mueller and Watson
######(2022, 2023), for feols regressions on longitude and latitude, without fixed effects,
######instruments or clusters, using conditional SCPC.
######
######Ported from the scpcR package (version 0.1.3, https://github.com/spatial-spur/scpcR),
######which is distributed under the following licence:
######
######  MIT License
######
######  Copyright 2025 Sascha Becker, P. David Boll, Hans-Joachim Voth
######
######  Permission is hereby granted, free of charge, to any person obtaining a copy
######  of this software and associated documentation files (the "Software"), to deal
######  in the Software without restriction, including without limitation the rights
######  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
######  copies of the Software, and to permit persons to whom the Software is
######  furnished to do so, subject to the following conditions:
######
######  The above copyright notice and this permission notice shall be included in all
######  copies or substantial portions of the Software.
######
######  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
######  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
######  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
######  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
######  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
######  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
######  SOFTWARE.

######Constants of scpcR: below 4500 observations the spatial weights are computed exactly,
######and above it from random subsets of locations, drawn with a fixed seed.
scpc_large_n=4500L
scpc_cgridfac=1.2
scpc_minavc=1e-05
scpc_seed=1

######40 point Gauss-Legendre nodes and weights on (0,1), from the eigenvalues of the
######Jacobi matrix (Golub and Welsch 1969), used to integrate rejection probabilities.
gauss_legendre_01=function(m){
  b=seq_len(m-1)/sqrt(4*seq_len(m-1)^2-1)
  J=matrix(0,m,m)
  J[cbind(1:(m-1),2:m)]=b
  J[cbind(2:m,1:(m-1))]=b
  e=eigen(J,symmetric=TRUE)
  ord=order(e$values,decreasing=TRUE)
  list(x=(e$values[ord]+1)/2,w=e$vectors[1,ord]^2)
}
scpc_gq=gauss_legendre_01(40)

######Spatial weights for one avc: the matrix W of the constant and the leading principal
######components of the spatial correlation matrix, the correlation matrices Oms of W over a
######grid of correlation decay rates, and the unconditional critical value cvfin.
######coords are longitude and latitude. large_n chooses the method for large samples.
scpc_spatial=function(coords,avc0,large_n=nrow(coords)>=scpc_large_n){
  coords=as.matrix(coords)
  n=nrow(coords)
  qmax=if(avc0>=0.05) 10 else if(avc0>=0.01) 20 else if(avc0>=0.005) 60 else 120
  distmat=NULL
  coords_use=coords
  perm=seq_len(n)
  random_t=NULL
  repeat{
    qmax=min(qmax,n-1)
    if(!large_n){
      distmat=scpc_distmat(coords)
      distv=distmat[lower.tri(distmat)]
      c0=scpc_c0_from_avc(distv,avc0)
      cmax=scpc_c0_from_avc(distv,scpc_minavc)
      W=scpc_W(distmat,c0,qmax)
      Oms=scpc_Oms(distmat,c0,cmax,W)
    }else{
      norm_s=scpc_normalize_s(coords)
      coords_use=norm_s$coords
      perm=norm_s$perm
      ln_w=scpc_ln_W(coords_use,avc0,qmax,random_t=scpc_seed)
      W=ln_w$W
      c0=ln_w$c0
      cmax=ln_w$cmax
      oms_res=scpc_ln_Oms(coords_use,c0,cmax,W,random_t=ln_w$random_t)
      Oms=oms_res$Oms
      random_t=oms_res$state
    }
    fin=scpc_final_W(Oms,W,qmax)
    if(fin$q<qmax||qmax==n-1) break
    qmax=round(qmax+qmax/2)
  }
  return(list(Wfin=fin$W,cvfin=fin$cv,Omsfin=Oms,c0=c0,cmax=cmax,coords=coords_use,
              perm=perm,distmat=distmat,large_n=large_n,random_state=random_t))
}

######Great-circle distances in units of half the earth's circumference.
scpc_distmat=function(coords){
  S=as.data.frame(coords)
  names(S)=c("lon","lat")
  return(geodist::geodist(S,measure="haversine")/(2*pi*6378137))
}

######The same distances between paired rows of S1 and S2.
scpc_distvec=function(S1,S2){
  lon1=S1[,1]*pi/180
  lat1=S1[,2]*pi/180
  lon2=S2[,1]*pi/180
  lat2=S2[,2]*pi/180
  dlon=0.5*(lon1-lon2)
  dlat=0.5*(lat1-lat2)
  return(asin(sqrt(sin(dlat)^2+cos(lat1)*cos(lat2)*sin(dlon)^2))/pi)
}

######Decay rate c0 at which the average correlation exp(-c0*distance) equals avc0.
scpc_c0_from_avc=function(dist,avc0){
  avc=function(c) mean(exp(-c*dist))
  c0=c1=10
  while(avc(c0)<avc0){
    c1=c0
    c0=0.5*c0
  }
  while(avc(c1)>avc0&c1<5000){
    c0=c1
    c1=2*c1
  }
  repeat{
    c=sqrt(c0*c1)
    if(avc(c)>avc0) c0=c else c1=c
    if(c1-c0<0.001) break
  }
  return(c)
}

scpc_demean=function(mat){
  mat=mat-rowMeans(mat)
  return(mat-matrix(colMeans(mat),nrow(mat),ncol(mat),byrow=TRUE))
}

scpc_W=function(distmat,c0,qmax){
  n=nrow(distmat)
  eig=eigen(scpc_demean(exp(-c0*distmat)),symmetric=TRUE)
  ord=order(eig$values,decreasing=TRUE)
  return(cbind(rep(1,n)/sqrt(n),eig$vectors[,ord[seq_len(qmax)],drop=FALSE]))
}

scpc_nc=function(c0,cmax) max(2,ceiling(log(cmax/c0)/log(scpc_cgridfac)))

scpc_Oms=function(distmat,c0,cmax,W){
  nc=scpc_nc(c0,cmax)
  Oms=vector("list",nc)
  c=c0
  for(i in seq_len(nc)){
    Oms[[i]]=crossprod(W,exp(-c*distmat)%*%W)
    c=c*scpc_cgridfac
  }
  return(Oms)
}

######Number of principal components q that gives the shortest expected confidence interval.
scpc_final_W=function(Oms,W,qmax){
  cvs=lengths=numeric(qmax)
  for(q in seq_len(qmax)){
    cvs[q]=scpc_cv(Oms,q,0.05)
    lengths[q]=cvs[q]*gamma(0.5*(q+1))/(sqrt(q)*gamma(0.5*q))
  }
  q_opt=which.min(lengths)
  return(list(W=W[,1:(q_opt+1),drop=FALSE],cv=cvs[q_opt],q=q_opt))
}

######Critical value of size level, the largest over the decay rates in Oms.
scpc_cv=function(Oms,q,level){
  rp=1
  i=1
  cv0=stats::qt(1-level/2,df=q)/sqrt(q)
  while(rp>level){
    cv1=cv0
    while(scpc_rp(Oms[[i]][1:(q+1),1:(q+1)],cv1)>level){
      cv0=cv1
      cv1=cv1+1/sqrt(q)
    }
    while(cv1-cv0>0.001/sqrt(q)){
      cv=0.5*(cv0+cv1)
      if(scpc_rp(Oms[[i]][1:(q+1),1:(q+1)],cv)>level) cv0=cv else cv1=cv
    }
    maxrp=scpc_maxrp(Oms,q,cv1)
    if(maxrp$i==i) break
    i=maxrp$i
    rp=maxrp$max
    cv0=cv1
  }
  return(cv1*sqrt(q))
}

######Largest rejection probability over the decay rates in Oms, and where it occurs.
scpc_maxrp=function(Oms,q,cv){
  rps=vapply(Oms,function(Om) scpc_rp(Om[1:(q+1),1:(q+1)],cv),0)
  return(list(max=max(rps),i=which.max(rps)))
}

######Rejection probability of the SCPC t test with critical value cv when the correlation
######matrix of W is Om.
scpc_rp=function(Om,cv){
  Omx=-cv^2*Om
  Omx[1,]=Om[1,]
  evals_raw=Re(eigen(Omx,only.values=TRUE)$values)
  evals=-evals_raw[evals_raw<0]
  if(!length(evals)) return(0)
  denom=max(evals_raw)
  if(!is.finite(denom)||denom<=0) return(0)
  evals=evals/denom
  tot=0
  for(j in seq_along(scpc_gq$x)){
    u=scpc_gq$x[j]
    arg=(1-u^2)*exp(sum(log1p(evals/(1-u^2))))
    if(!is.finite(arg)||arg<=0) next
    tot=tot+scpc_gq$w[j]/sqrt(arg)
  }
  return(as.numeric(tot*2/pi))
}

######Large samples: locations sorted by latitude, then longitude, and spatial weights
######computed from random subsets of locations, using scpcR's random number generator so
######that results are the same.
scpc_normalize_s=function(S){
  S[,1]=S[,1]-mean(S[,1])
  S[,1]=((S[,1]+180)%%360)-180
  perm=order(S[,2],S[,1])
  return(list(coords=S[perm,,drop=FALSE],perm=perm))
}

scpc_next_u=function(random_t){
  random_t=(64389*random_t+1)%%2^32
  return(list(value=random_t/2^32,state=random_t))
}

scpc_jumble_s=function(S,m,random_t){
  n=nrow(S)
  for(i in seq_len(m)){
    nxt=scpc_next_u(random_t)
    random_t=nxt$state
    j=floor(nxt$value*n)+1L
    tmp=S[j,]
    S[j,]=S[i,]
    S[i,]=tmp
  }
  return(list(coords=S,state=random_t))
}

scpc_ln_W=function(S,avc0,qmax,random_t,capN=20L,m=1000L){
  n=nrow(S)
  m=min(m,n)
  ms=vector("list",capN)
  block_len=m*(m-1)/2
  distvec=numeric(capN*block_len)
  r=S
  for(i in seq_len(capN)){
    jumbled=scpc_jumble_s(r,m,random_t)
    r=jumbled$coords
    random_t=jumbled$state
    ms[[i]]=list(coords=r[seq_len(m),,drop=FALSE],distmat=scpc_distmat(r[seq_len(m),,drop=FALSE]))
    distvec[((i-1L)*block_len+1L):(i*block_len)]=ms[[i]]$distmat[lower.tri(ms[[i]]$distmat)]
  }
  c0=scpc_c0_from_avc(distvec,avc0)
  cmax=scpc_c0_from_avc(distvec,scpc_minavc)
  Wall=matrix(0,n,capN*qmax)
  for(i in seq_len(capN)){
    W0=eigen(scpc_demean(exp(-c0*ms[[i]]$distmat)),symmetric=TRUE)$vectors[,seq_len(qmax),drop=FALSE]
    Wx=matrix(0,n,qmax)
    for(j in seq_len(m)){
      diff=sweep(S,2,ms[[i]]$coords[j,],`-`)
      v=exp(-c0*sqrt(rowSums(diff^2)))
      Wx=Wx+tcrossprod(v,W0[j,])
    }
    Wx=sweep(Wx,2,colMeans(Wx),`-`)
    norms=sqrt(colSums(Wx^2))
    norms[!is.finite(norms)|norms==0]=1
    Wx=sweep(Wx,2,norms,`/`)
    for(j in seq_len(qmax)) Wall[,(j-1L)*capN+i]=Wx[,j]
  }
  W=matrix(0,n,qmax)
  for(i in seq_len(qmax)){
    Wx=Wall[,seq_len(capN*i),drop=FALSE]
    evec=eigen(crossprod(Wx),symmetric=TRUE)$vectors[,1,drop=FALSE]
    W[,i]=as.numeric(Wx%*%evec)
    W[,i]=W[,i]/sqrt(sum(W[,i]^2))
    Wall=Wall-W[,i,drop=FALSE]%*%crossprod(W[,i],Wall)
  }
  return(list(W=cbind(rep(1/sqrt(n),n),W),c0=c0,cmax=cmax,random_t=random_t))
}

scpc_raninds=function(n,capM,random_t){
  v=numeric(capM+1L)
  nxt=scpc_next_u(random_t)
  random_t=nxt$state
  j=floor(n*nxt$value)
  for(i in seq_len(capM+1L)){
    v[i]=j+1L
    nxt=scpc_next_u(random_t)
    random_t=nxt$state
    j=(j+1+floor(nxt$value*(n-1)))%%n
  }
  return(list(indices=as.integer(v),state=random_t))
}

scpc_ln_Oms=function(S,c0,cmax,W,random_t,capM=1000000L){
  nc=scpc_nc(c0,cmax)
  Oms=vector("list",nc)
  n=nrow(S)
  inds_res=scpc_raninds(n,capM,random_t)
  inds=inds_res$indices
  dist=scpc_distvec(S[inds[seq_len(capM)],,drop=FALSE],S[inds[2:(capM+1L)],,drop=FALSE])
  W1=W[inds[seq_len(capM)],,drop=FALSE]
  W2=W[inds[2:(capM+1L)],,drop=FALSE]
  c=c0
  for(i in seq_len(nc)){
    cd=exp(-c*dist)
    Oms[[i]]=diag(ncol(W))+0.5*(n*(n-1)/capM)*(crossprod(W1,W2*cd)+crossprod(W2,W1*cd))
    c=c*scpc_cgridfac
  }
  return(list(Oms=Oms,state=inds_res$state))
}

######SCPC estimate, standard error, t statistic, p value and 95% confidence interval of the
######coefficients coefs (positions in coef(fit)) of a feols regression, given the spatial
######weights spc from scpc_spatial() for the same observations. The approximate method for
######large samples draws random numbers for each coefficient in turn, so all coefficients
######up to the last one wanted are computed then.
scpc_coef_stats=function(fit,spc,coefs){
  X=stats::model.matrix(fit)
  coef_vec=stats::coef(fit)
  X=X[,names(coef_vec),drop=FALSE]
  n=nrow(X)
  w=stats::weights(fit)
  if(is.null(w)) w=rep(1,n)
  ##Influence of each observation on the coefficients, and the coefficients' weights on
  ##the regressors.
  XWXinv=solve(crossprod(X,w*X))
  infl=n*XWXinv%*%t(X*(w*stats::resid(fit)))
  xw=n*XWXinv%*%t(X)

  Wfin=spc$Wfin
  q=ncol(Wfin)-1
  perm=spc$perm
  random_state=spc$random_state
  Xr=if("(Intercept)"%in%colnames(X)) X else cbind(`(Intercept)`=1,X)
  qrX=qr(Xr[perm,,drop=FALSE],tol=1e-10)
  if(qrX$rank<ncol(X)) stop("SCPC regressors are collinear.",call.=FALSE)

  out=matrix(NA_real_,length(coefs),6)
  for(j in if(spc$large_n) seq_len(max(coefs)) else coefs){
    coef_j=coef_vec[[j]]
    wj=(infl[j,]+coef_j)[perm]
    denom=sqrt(sum((t(Wfin[,-1,drop=FALSE])%*%wj)^2))
    if(!is.finite(denom)||denom<=sqrt(.Machine$double.eps))
      stop("SCPC variance projection is degenerate.",call.=FALSE)
    tau_u=as.numeric(sqrt(q)*crossprod(Wfin[,1],wj)/denom)
    SE=as.numeric(denom/(sqrt(q)*sqrt(n)))
    p_u=scpc_maxrp(spc$Omsfin,q,abs(tau_u)/sqrt(q))$max

    ##Conditional SCPC: W projected off the regressors.
    xj=xw[j,perm]
    xjs=sign(xj)
    Wx=Wfin
    Wx[,1]=Wx[,1]*xj*xjs
    if(ncol(Wx)>1L) Wx[,-1]=qr.resid(qrX,Wx[,-1,drop=FALSE]*xjs)*xj
    if(!all(is.finite(Wx)))
      stop("SCPC conditional projection produced non-finite values.",call.=FALSE)
    if(spc$large_n){
      omsx_res=scpc_ln_Oms(spc$coords,spc$c0,spc$cmax,Wx,random_t=random_state)
      Omsx=omsx_res$Oms
      random_state=omsx_res$state
    }else{
      Omsx=scpc_Oms(spc$distmat,spc$c0,spc$cmax,Wx)
    }
    if(!j%in%coefs) next
    p=max(p_u,scpc_maxrp(Omsx,q,abs(tau_u)/sqrt(q))$max)
    cv=max(spc$cvfin,scpc_cv(Omsx,q,0.05))
    out[match(j,coefs),]=c(coef_j,SE,tau_u,p,coef_j-cv*SE,coef_j+cv*SE)
  }
  dimnames(out)=list(names(coef_vec)[coefs],c("Coef","Std_Err","t","p","lower","upper"))
  return(out)
}
