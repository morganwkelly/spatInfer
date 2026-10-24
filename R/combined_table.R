#' Table combining chosen placebo and synthetic outcome results.
#'
#'@description
#' `combined_table()` makes one summary table from a [placebo()] test and the matching [synth()]
#' test, or from the output of [placebo_synth()]. It starts with the HC row and then, for each inference method chosen, gives the row with
#' the chosen number of clusters (BCH and IM) or SCPC bound `avc`, with the synthetic outcome p
#' value from the same row of the synth. [iv_combined_table()] does the same for instrumental
#' variable regressions.
#'
#' @param plbo An object created by [placebo()], or by [placebo_synth()], in which case `syn` is
#'   left out.
#' @param syn An object created by [synth()] for the same regression and spatial basis. Leave it
#'   out when `plbo` comes from [placebo_synth()].
#' @param select A named numeric vector giving, for each inference method to show, the number of
#'   clusters (`BCH` and `IM`) or the value of `avc` (`SCPC`) of its row, for example
#'   `c(BCH=7, IM=5, SCPC=0.04)`. Rows appear in this order. SCPC needs `avc` in [placebo()] and
#'   [synth()].
#' @param caption An optional caption for the table.
#'
#' @return Generates a tinytable object. Its first row gives the HC regression, placebo and
#'   synthetic outcome p values, and each following row those of a chosen inference method, with the
#'   regression confidence interval and its width. The notes give the z-score of a Moran test of the
#'   regression residuals, using 5 nearest neighbours, and, separately for the treatment and the
#'   outcome, the R2, effective range and structure used to simulate noise: the R2 of a regression of
#'   the treatment on the principal components of the tensor, and of the outcome on a quadratic in
#'   longitude and latitude.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' # Use 100 observations and 100 simulations to speed things up.
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=100)
#' fm=mobility~single_mothers+short_commute+gini+dropout_rate+social_cap+dropout_na
#' plbo=placebo(fm, opportunity, splines=6, pc_num=15, nSim=100, Parallel=FALSE, max_clus=7, avc=0.04)
#' synt=synth(fm, opportunity, splines=6, pc_num=15, nSim=100, Parallel=FALSE, max_clus=7, avc=0.04)
#' combined_table(plbo, synt, select=c(BCH=5, IM=4, SCPC=0.04),
#'   caption="Placebo and synthetic outcome significance levels for single mothers variable.")
#' # Or both tests from one command.
#' ps=placebo_synth(fm, opportunity, splines=6, pc_num=15, nSim=100, Parallel=FALSE, max_clus=7,
#'   avc=0.04, progress=FALSE)
#' combined_table(ps, select=c(BCH=5, IM=4, SCPC=0.04))


combined_table=function(plbo, syn=NULL, select, caption=""){
if(missing(select)) select=NULL
inp=unpack_placebo_synth(plbo,syn,select,"placebo_synth()")
plbo=inp$plbo; syn=inp$syn; select=inp$select
check_combined_input(plbo,"Plac p","plbo","placebo()")
check_combined_input(syn,"Synth p","syn","synth()")
check_same_regression(list(plbo,syn))
check_select(select,"select")

#####HC row first, then the chosen row of each inference method with the synth p value of the same row.
rr=combined_rows(plbo,syn,select,"")
rr=rr |>
  dplyr::rename(Adj=SE) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric)&!dplyr::any_of("CI Width"),~{round(.x,3)}))   #CI widths keep their significant digits
names(rr)[2]=level_name(names(select))

sp=plbo$Spatial_Params|>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric),~{round(.x,2)}))
ss=syn$Spatial_Params|>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric),~{round(.x,2)}))

comb_tab=tinytable::tt(rr,theme="striped",align="r",digits=3,
              notes=list(
                paste0("Moran=",sp$Moran,
                       ", N=",ss$N,
                       ", Splines=",sp$Splines,
                       ", PCs=",sp$PCs,"."
                ),
                paste0("Treatment: Structure=",sp$Structure,
                       ", Effective Range=",sp$Effective_Range,
                       ", R2=",sp$R2,"."
                ),
                paste0("Outcome: Structure=",ss$Structure,
                       ", Effective Range=",ss$Effective_Range,
                       ", R2=",ss$R2,"."
                ),
                "Estimated and placebo p values and proportion of placebo regressions significant at 5%,
                synthetic outcome p values, and confidence intervals, for the chosen number of clusters
                or SCPC average correlation bound avc. Treatment R2 gives the explanatory power of a
                regression of the treatment on the principal components, and outcome R2 of a regression
                of the outcome on a quadratic in longitude and latitude."
              ),
              caption=caption
)

comb_tab= comb_tab |>
  tinytable::style_tt(i = 0, bold=TRUE ) |>
  tinytable::style_tt(j=1, bold=TRUE ) |>
  tinytable::style_tt(j = c("Plac 5%"), color = "orange")

return(comb_tab)
}

#####combined_table() and iv_combined_table() take a placebo and a synth, or the output of
#####placebo_synth() alone, in which case select may also be given second, in the place of syn.
unpack_placebo_synth=function(plbo,syn,select,fun){
  is_ps=is.list(plbo) && is.null(plbo$Results) && all(c("Placebo","Synth") %in% names(plbo))
  if(!is_ps) return(list(plbo=plbo,syn=syn,select=select))
  if(is.null(select) && is.numeric(syn)){
    select=syn
    syn=NULL
  }
  if(!is.null(syn))
    stop("`syn` must be left out when `plbo` is the output of ",fun,".",call.=FALSE)
  list(plbo=plbo$Placebo,syn=plbo$Synth,select=select)
}

#####Placebo and synth results with the expected simulated p value column. iv is TRUE for
#####iv_combined_table(), which takes only IV results; combined_table() takes none.
check_combined_input=function(x,p_col,arg,fun,iv=FALSE){
  ok=is.list(x) && is.data.frame(x$Results) && p_col %in% names(x$Results) &&
    !is.null(x$Spatial_Params)
  if(!ok || (iv && !is_iv_results(x$Results)))
    stop("`",arg,"` is not the output of ",fun,".",call.=FALSE)
  if(!iv && is_iv_results(x$Results))
    stop("combined_table() does not take iv_placebo() or iv_synth() results; use iv_combined_table().",
         call.=FALSE)
}

#####Named numbers of clusters or avc values, one for each of BCH, IM and SCPC at most.
check_select=function(sel,arg){
  ok=is.numeric(sel) && length(sel)>0 && !anyNA(sel) && !is.null(names(sel)) &&
    all(names(sel) %in% c("BCH","IM","SCPC")) && !anyDuplicated(names(sel))
  if(!ok)
    stop("`",arg,"` must be a named vector of numbers of clusters or avc values, such as c(BCH=7, IM=5, SCPC=0.04).",
         call.=FALSE)
}

#####Header of the column of numbers of clusters or avc values, given the inference methods shown.
level_name=function(methods){
  if(all(methods=="SCPC")) "avc" else if("SCPC" %in% methods) "Clusters/avc" else "Clusters"
}

#####All tests must come from the same regression and spatial basis, so that they share
#####one HC row and one set of notes.
check_same_regression=function(xs){
  hc_p=vapply(xs,function(x) x$Results$`Est p`[x$Results$SE=="HC"][1],0)
  basis=vapply(xs,function(x) paste(x$Spatial_Params$Splines,x$Spatial_Params$PCs),"")
  if(any(abs(hc_p-hc_p[1])>1e-12) || any(basis!=basis[1]))
    stop("All placebos and synths must come from the same regression, with the same splines and PCs.",
         call.=FALSE)
}

#####Row of results whose SE label is `se` and whose clusters or avc (the second column, stored as
#####character with "." for HC) equal `value`. value is NA for HC rows, which have no level.
pick_row=function(results,se,value){
  level=results[[2]]
  rows=results$SE==se
  if(!is.na(value)) rows=rows & (abs(suppressWarnings(as.numeric(level))-value)<1e-8) %in% TRUE
  rows=which(rows %in% TRUE)
  if(length(rows)==0){
    available=paste(level[results$SE==se],collapse=", ")
    stop("No ",se," row with ",names(results)[2]," = ",value,". Available values: ",available,".",
         call.=FALSE)
  }
  results[rows[1],,drop=FALSE]
}

#####HC row, then the chosen row of each inference method in sel, with the synth p value of the same
#####row. prefix is "" for least squares and 2SLS, and "AR " for Anderson-Rubin, whose rows are
#####labelled by their inference alone. F is kept for IV results.
combined_rows=function(plbo,syn,sel,prefix){
  rows=lapply(c("HC",names(sel)),function(se){
    value=if(se=="HC") NA else sel[[se]]
    row=pick_row(plbo$Results,paste0(prefix,se),value)
    row$`Synth p`=pick_row(syn$Results,paste0(prefix,se),value)$`Synth p`
    row$SE=se
    names(row)[2]="Level"
    row[,intersect(c("SE","Level","F","Est p","Plac p","Plac 5%","Synth p","CI Width","CI"),names(row))]
  })
  rr=do.call(rbind,rows)
  rownames(rr)=NULL
  return(rr)
}
