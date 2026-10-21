#' Table combining chosen IV placebo and synthetic outcome results.
#'
#'@description
#' `iv_combined_table()` is the instrumental variable version of [combined_table()]. It makes
#' one summary table from an [iv_placebo()] test and the matching [iv_synth()] test, with a block
#' of 2SLS results followed by a block of Anderson-Rubin (AR) results. Each block starts with the
#' HC row and then, for each inference method chosen, the row with the chosen number of clusters
#' (BCH and IM) or SCPC bound `avc`, with the synthetic outcome p value from the same row of the
#' synth.
#'
#' @param plbo An object created by [iv_placebo()].
#' @param syn An object created by [iv_synth()] for the same regression and spatial basis.
#' @param select A named numeric vector giving, for each inference method to show, the number of
#'   clusters (`BCH` and `IM`) or the value of `avc` (`SCPC`) of the 2SLS row, for example
#'   `c(BCH=7, IM=5, SCPC=0.04)`. Rows appear in this order. SCPC needs `avc` in [iv_placebo()]
#'   and [iv_synth()].
#' @param ar_select The same for the AR rows. Defaults to `select`.
#' @param caption An optional caption for the table.
#'
#' @return Generates a tinytable object with the 2SLS and AR blocks. Each row gives the first stage
#'   F statistic with its inference (blank for IM and SCPC), the estimated, placebo and synthetic
#'   outcome p values, the proportion of placebo regressions significant at 5%, and the confidence
#'   interval (an AR confidence set for AR rows) and its width. The notes give the z-score of a
#'   Moran test of the reduced form residuals, and the R2, effective range and structure used to
#'   simulate noise for the treatment, the instrument and the outcome.
#' @export
#'
#' @examples
#' \donttest{
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' fm=d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag
#' # Use 100 simulations and turn off parallel processing.
#' plbo=iv_placebo(fm,china,splines=6,pc_num=13,nSim=100,weights="timepwt48",Parallel=FALSE,avc=0.04)
#' syn=iv_synth(fm,china,splines=6,pc_num=13,nSim=100,weights="timepwt48",Parallel=FALSE,avc=0.04)
#' iv_combined_table(plbo,syn,select=c(BCH=6,IM=5,SCPC=0.04))
#' }

iv_combined_table=function(plbo, syn, select, ar_select=select, caption=""){
check_iv_combined_input(plbo,"Plac p","plbo","iv_placebo()")
check_iv_combined_input(syn,"Synth p","syn","iv_synth()")
check_same_regression(list(plbo,syn))
check_iv_select(select,"select")
check_iv_select(ar_select,"ar_select")

#####Block of rows of one test: HC, then the chosen row of each inference method, with the synth
#####p value of the same row. prefix is "" for 2SLS and "AR " for Anderson-Rubin.
block=function(sel,prefix){
  rows=lapply(c("HC",names(sel)),function(se){
    value=if(se=="HC") NA else sel[[se]]
    row=pick_row(plbo$Results,paste0(prefix,se),value)
    row$`Synth p`=pick_row(syn$Results,paste0(prefix,se),value)$`Synth p`
    row$SE=se
    names(row)[2]="Level"
    row[,c("SE","Level","F","Est p","Plac p","Plac 5%","Synth p","CI Width","CI")]
  })
  do.call(rbind,rows)
}
wald=block(select,"")
ar=block(ar_select,"AR ")
rr=rbind(wald,ar)
rownames(rr)=NULL
F=rep("",nrow(rr))   #IM and SCPC rows have no first stage F
F[!is.na(rr$F)]=trimws(format(rr$F[!is.na(rr$F)],digits=3))
rr$F=F
rr=rr |>
  dplyr::rename(Adj=SE) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric)&!dplyr::any_of("CI Width"),~{round(.x,3)}))   #CI widths keep their significant digits
methods=union(names(select),names(ar_select))
names(rr)[2]=if(all(methods=="SCPC")) "avc" else if("SCPC" %in% methods) "Clusters/avc" else "Clusters"

sp=plbo$Spatial_Params|>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric),~{round(.x,2)}))
ss=syn$Spatial_Params|>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric),~{round(.x,2)}))

comb_tab=tinytable::tt(rr,theme="striped",align="r",digits=3,
              notes=list(
                paste0("Moran=",sp$Moran,
                       ", N=",sp$N,
                       ", Splines=",sp$Splines,
                       ", PCs=",sp$PCs,"."
                ),
                paste0("Treatment: Structure=",sp$Structure,
                       ", Effective Range=",sp$Effective_Range,
                       ", R2=",sp$R2,". ",
                       "Instrument: Structure=",sp$Structure_Inst,
                       ", Effective Range=",sp$Effective_Range_Inst,
                       ", R2=",sp$R2_Inst,"."
                ),
                paste0("Outcome: Structure=",ss$Structure,
                       ", Effective Range=",ss$Effective_Range,
                       ", R2=",ss$R2,"."
                ),
                "2SLS estimates and Anderson-Rubin tests of a zero treatment effect, with estimated and
                placebo p values, proportion of placebo regressions significant at 5%, synthetic outcome
                p values, and confidence intervals or Anderson-Rubin confidence sets, for the chosen number
                of clusters or SCPC average correlation bound avc. F is the first stage F statistic, not
                given for IM or SCPC. Treatment R2 gives the explanatory power of a regression of the
                treatment on the principal components, instrument R2 of the instrument on the treatment,
                and outcome R2 of the outcome on a quadratic in longitude and latitude."
              ),
              caption=caption
)

comb_tab= comb_tab |>
  tinytable::style_tt(i = 0, bold=TRUE ) |>
  tinytable::style_tt(j=1, bold=TRUE ) |>
  tinytable::style_tt(j = c("Plac 5%"), color = "orange") |>
  tinytable::group_tt(i=stats::setNames(list(1,nrow(wald)+1),c("2SLS","Anderson-Rubin")))

return(comb_tab)
}

#####iv_placebo() or iv_synth() results, with the expected simulated p value column.
check_iv_combined_input=function(x,p_col,arg,fun){
  ok=is.list(x) && is.data.frame(x$Results) && p_col %in% names(x$Results) &&
    is_iv_results(x$Results) && !is.null(x$Spatial_Params)
  if(!ok) stop("`",arg,"` is not the output of ",fun,".",call.=FALSE)
}

#####Named numbers of clusters or avc values, one for each of BCH, IM and SCPC at most.
check_iv_select=function(sel,arg){
  ok=is.numeric(sel) && length(sel)>0 && !anyNA(sel) && !is.null(names(sel)) &&
    all(names(sel) %in% c("BCH","IM","SCPC")) && !anyDuplicated(names(sel))
  if(!ok)
    stop("`",arg,"` must be a named vector of numbers of clusters or avc values, such as c(BCH=7, IM=5, SCPC=0.04).",
         call.=FALSE)
}
