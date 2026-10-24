#' Regression table comparing instrumental variable estimates.
#'
#' @description
#' `iv_regression_table()` makes a [modelsummary::modelsummary()] table of
#' instrumental variable regressions, with a column for each of the BCH, IM and
#' SCPC results supplied, optionally preceded by a column for the original
#' regression. Below the coefficients it gives, for each column, the first stage
#' F statistic and the Anderson-Rubin (AR) p value and 95% confidence set of the
#' treatment effect, all with the inference and number of clusters of that
#' column. No first stage F is given for IM or SCPC. Notes give the number of
#' clusters and `avc`.
#'
#' @param orig An optional [fixest::feols()] IV regression, typically the
#'   original regression without a spatial basis, shown in a first column
#'   labelled "Original" with the standard errors it was estimated with. Its
#'   first stage F statistic uses the same standard errors; no AR results are
#'   shown for it.
#' @param bch Output of [iv_bch_regression()], with BCH or HC standard errors. The
#'   column is labelled "BCH" or "HC".
#' @param im Output of [iv_im_regression()].
#' @param scpc Output of [iv_scpc_regression()].
#' @inheritParams regression_table
#' @param notes Further notes, added after the notes on clusters.
#'
#' @details At least one of `orig`, `bch`, `im` and `scpc` must be given. The IM
#'   estimate is the mean of the 2SLS coefficients estimated in each cluster, so
#'   it differs from the other estimates.
#'
#' @return A table from [modelsummary::modelsummary()], in the format it
#'   chooses or the one set by `output`.
#' @export
#'
#' @examples
#' \donttest{
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' fm=d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag
#' orig=fixest::feols(fm,china,weights=~timepwt48,vcov="hetero")
#' iv_bch=iv_bch_regression(fm,china,splines=6,pc_num=8,clusters=6,weights="timepwt48")
#' iv_im=iv_im_regression(fm,china,splines=6,pc_num=8,clusters=4,weights="timepwt48")
#' iv_sc=iv_scpc_regression(fm,china,splines=6,pc_num=8,avc=0.05,weights="timepwt48")
#' iv_regression_table(orig=orig,bch=iv_bch,im=iv_im,scpc=iv_sc)
#' iv_regression_plot(orig=orig,bch=iv_bch,im=iv_im,scpc=iv_sc)
#' }

iv_regression_table=function(orig=NULL,bch=NULL,im=NULL,scpc=NULL,statistic="conf.int",
                             coef_omit="Intercept|^PC[0-9]+$",gof_map="nobs",notes=NULL,...){
  res=iv_regression_models(orig,bch,im,scpc,"iv_regression_table()")
  rows=data.frame(term=c("First stage F","AR p","AR 95% set"),check.names=FALSE)
  for(name in names(res$models)){
    ar=res$ar[[name]]
    rows[[name]]=c(format_F(res$F[[name]]),
                   if(is.null(ar)) "" else format(signif(ar$AR$p.value,3)),
                   if(is.null(ar)) "" else ar$AR$CI)
  }
  return(modelsummary::modelsummary(res$models,
                                    statistic=statistic,
                                    coef_omit=coef_omit,
                                    gof_map=gof_map,
                                    add_rows=rows,
                                    notes=c(list(paste(res$notes,collapse=" ")),as.list(notes)),
                                    ...))
}

#' Plot comparing instrumental variable estimates.
#'
#' @description
#' `iv_regression_plot()` plots the coefficients and 95% confidence intervals
#' of the instrumental variable regressions supplied, one panel per
#' coefficient. For the treatment it adds the Anderson-Rubin (AR) 95%
#' confidence set of each spatial basis regression.
#' Unbounded parts of an AR set run to the edge of the panel and end in arrows.
#'
#' @inheritParams iv_regression_table
#' @param coef_omit Coefficients to leave out, as a regular expression. Defaults
#'   to the intercept and the principal components of the spatial basis.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' fm=d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag
#' iv_bch=iv_bch_regression(fm,china,splines=6,pc_num=8,clusters=6,weights="timepwt48")
#' iv_regression_plot(bch=iv_bch)

iv_regression_plot=function(orig=NULL,bch=NULL,im=NULL,scpc=NULL,coef_omit="Intercept|^PC[0-9]+$"){
  res=iv_regression_models(orig,bch,im,scpc,"iv_regression_plot()",need_modelsummary=FALSE)
  model_names=names(res$models)
  model_y=function(name) -match(name,model_names)   #first model at the top

  coefs=purrr::list_rbind(lapply(model_names,function(name){
    tidy=res$models[[name]]$tidy
    data.frame(model=name,term=tidy$term,estimate=tidy$estimate,conf.low=tidy$conf.low,conf.high=tidy$conf.high)
  }))
  coefs=coefs[!grepl(coef_omit,coefs$term),,drop=FALSE]
  coefs$y=model_y(coefs$model)+0.12
  coefs$interval="2SLS 95% CI"

  segs=purrr::list_rbind(lapply(names(res$ar),function(name){
    ar=res$ar[[name]]
    ar_segments(ar$AR_Set,ar$AR$term,model_y(name)-0.12)
  }))
  if(is.null(segs)) segs=data.frame(term=character(0),x=numeric(0),xend=numeric(0),y=numeric(0))
  segs=segs[segs$term %in% coefs$term,,drop=FALSE]
  segs$arrow=character(nrow(segs))
  if(nrow(segs)>0){
    #unbounded ends run to the edge of the treatment panel, beyond everything else shown in it
    treat=segs$term[1]
    finite=c(unlist(coefs[coefs$term==treat,c("conf.low","conf.high")]),segs$x[is.finite(segs$x)],
             segs$xend[is.finite(segs$xend)])
    span=diff(range(finite))
    if(span==0) span=max(1,abs(finite[1]))
    lo=min(finite)-0.15*span
    hi=max(finite)+0.15*span
    segs$arrow=ifelse(is.infinite(segs$x),ifelse(is.infinite(segs$xend),"both","first"),
                      ifelse(is.infinite(segs$xend),"last","none"))
    segs$x=pmax(segs$x,lo)
    segs$xend=pmin(segs$xend,hi)
    segs$interval="AR 95% set"
  }

  term_levels=unique(coefs$term)   #panels in the order of the coefficients
  coefs$term=factor(coefs$term,levels=term_levels)
  segs$term=factor(segs$term,levels=term_levels)

  colours=c("2SLS 95% CI"="#0072B2","AR 95% set"="#D55E00")
  linetypes=c("2SLS 95% CI"="solid","AR 95% set"="dashed")
  p=ggplot2::ggplot()+
    ggplot2::geom_vline(xintercept=0,linetype="dashed",colour="grey50",linewidth=0.3)+
    ggplot2::geom_linerange(data=coefs,ggplot2::aes(xmin=conf.low,xmax=conf.high,y=y,
                                                     colour=interval,linetype=interval),linewidth=0.7)+
    ggplot2::geom_point(data=coefs,ggplot2::aes(x=estimate,y=y,colour=interval),size=2)
  for(ends in c("none","first","last","both")){
    d=segs[segs$arrow==ends,,drop=FALSE]
    if(nrow(d)==0) next
    arrow=if(ends=="none") NULL else ggplot2::arrow(length=ggplot2::unit(0.08,"inches"),ends=ends,type="closed")
    p=p+ggplot2::geom_segment(data=d,ggplot2::aes(x=x,xend=xend,y=y,yend=y,colour=interval,linetype=interval),
                              linewidth=0.7,arrow=arrow)
  }
  p=p+ggplot2::facet_wrap(~term,scales="free_x")+
    ggplot2::scale_y_continuous(breaks=-seq_along(model_names),labels=model_names,
                                limits=c(-length(model_names)-0.5,-0.5))+
    ggplot2::scale_colour_manual(values=colours,breaks=names(colours),name=NULL)+
    ggplot2::scale_linetype_manual(values=linetypes,breaks=names(linetypes),name=NULL)+
    ggplot2::labs(x="Coefficient",y=NULL,caption=paste(strwrap(paste(res$notes,collapse=" "),width=100),collapse="\n"))+   #long notes wrapped to fit
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="bottom",panel.grid.minor=ggplot2::element_blank(),
                   panel.grid.major.y=ggplot2::element_blank())
  return(p)
}

######Line segments of an AR confidence set at height y: an interval, two rays or the whole line.
######Infinite ends are kept as -Inf and Inf, for iv_regression_plot() to clip and mark with arrows.
ar_segments=function(set,term,y){
  if(length(set)==4)
    return(data.frame(term=term,x=c(-Inf,set[3]),xend=c(set[2],Inf),y=y))
  return(data.frame(term=term,x=set[1],xend=set[2],y=y))
}

format_F=function(F) if(is.na(F)) "" else format(signif(F,3))

######Check the regressions supplied and put them in a named list modelsummary accepts, with their
######first stage F statistics, AR results and notes on clusters.
iv_regression_models=function(orig,bch,im,scpc,caller,need_modelsummary=TRUE){
  if(need_modelsummary&&!requireNamespace("modelsummary",quietly=TRUE))
    stop(caller," needs the modelsummary package. Install it with install.packages(\"modelsummary\").",
         call.=FALSE)
  if(is.null(orig)&&is.null(bch)&&is.null(im)&&is.null(scpc))
    stop("Give at least one of orig, bch, im and scpc.",call.=FALSE)

  models=list()
  F=list()
  ar=list()
  model_notes=character(0)
  if(!is.null(orig)){
    if(!inherits(orig,"fixest"))
      stop("orig must be a fixest::feols() regression.",call.=FALSE)
    models$Original=as_modelsummary_list(broom::tidy(orig,conf.int=TRUE),stats::nobs(orig))
    F$Original=if(isTRUE(orig$is_iv)) fixest::fitstat(orig,"ivwald1",simplify=TRUE)$stat else NA
  }
  if(!is.null(bch)){
    if(!inherits(bch,"iv_regression")||inherits(bch,c("iv_im_regression","iv_scpc_regression")))
      stop("bch must be the output of iv_bch_regression().",call.=FALSE)
    name=bch$Inference
    models[[name]]=as_modelsummary_list(bch$Coefs,bch$N)
    F[[name]]=bch$F
    ar[[name]]=bch
    model_notes=c(model_notes,if(name=="HC") "HC: heteroskedasticity robust standard errors."
                  else paste0("BCH: ",bch$Clusters," clusters."))
  }
  if(!is.null(im)){
    if(!inherits(im,"iv_im_regression"))
      stop("im must be the output of iv_im_regression().",call.=FALSE)
    models$IM=as_modelsummary_list(im$Coefs,im$N)
    F$IM=im$F   #NA: no first stage F for IM
    ar$IM=im
    model_notes=c(model_notes,paste0("IM: ",im$Clusters," clusters."))
  }
  if(!is.null(scpc)){
    if(!inherits(scpc,"iv_scpc_regression"))
      stop("scpc must be the output of iv_scpc_regression().",call.=FALSE)
    models$SCPC=as_modelsummary_list(scpc$Coefs,scpc$N)
    F$SCPC=scpc$F   #NA: no first stage F for SCPC
    ar$SCPC=scpc
    model_notes=c(model_notes,paste0("SCPC: avc = ",scpc$avc,"."))
  }
  model_notes=c(model_notes,paste0("First stage F and Anderson-Rubin (AR) results use the inference of each column",
                                   if(!is.null(im)||!is.null(scpc))
                                     paste0("; no F is given for ",paste(c(if(!is.null(im)) "IM",if(!is.null(scpc)) "SCPC"),collapse=" or "))
                                   else "","."))
  return(list(models=models,F=F,ar=ar,notes=model_notes))
}
