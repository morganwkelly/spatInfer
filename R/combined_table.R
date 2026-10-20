#' Table combining chosen placebo and synthetic outcome results.
#'
#'@description
#' `combined_table()` makes one summary table from several placebo tests and their matching
#' synthetic outcome tests. From each placebo it takes the row with a chosen number of clusters
#' (for [placebo()] and [placebo_im()]) or SCPC bound `avc` (for [placebo_scpc()]), and adds the
#' synthetic outcome p value from the same row of the matching synth. The first row gives the
#' HC results.
#'
#' @param placebos A list of objects created by [placebo()], [placebo_im()] or [placebo_scpc()],
#'   or a single such object.
#' @param synths A list of objects created by [synth()], [synth_im()] or [synth_scpc()], in the
#'   same order as `placebos` and using the same estimator as each, or a single such object.
#' @param select A numeric vector, one entry for each placebo, giving the number of clusters (BCH
#'   and IM) or the value of `avc` (SCPC) of the row to show.
#' @param caption An optional caption for the table.
#'
#' @return Generates a tinytable object. Its first row gives the HC regression, placebo and
#'   synthetic outcome p values, and each following row those of a chosen placebo and synth, with the
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
#' plbo_bch=placebo(fm, opportunity, splines=6, pc_num=15, nSim=100, Parallel=FALSE, max_clus=7)
#' synt_bch=synth(fm, opportunity, splines=6, pc_num=15, nSim=100, Parallel=FALSE, max_clus=7)
#' plbo_im=placebo_im(fm, opportunity, splines=6, pc_num=15, nSim=100, Parallel=FALSE, max_clus=7)
#' synt_im=synth_im(fm, opportunity, splines=6, pc_num=15, nSim=100, Parallel=FALSE, max_clus=7)
#' combined_table(list(plbo_bch, plbo_im), list(synt_bch, synt_im), select=c(5, 4),
#'   caption="Placebo and synthetic outcome significance levels for single mothers variable.")


combined_table=function(placebos, synths, select, caption=""){
if(!is.null(placebos$Results)) placebos=list(placebos)   #a single placebo
if(!is.null(synths$Results)) synths=list(synths)
if(length(placebos)!=length(synths) || length(placebos)!=length(select))
  stop("`placebos`, `synths` and `select` must have the same length.", call.=FALSE)
if(!is.numeric(select) || anyNA(select))
  stop("`select` must give a number of clusters or an `avc` value for each placebo.", call.=FALSE)
for(i in seq_along(placebos)){
  check_combined_input(placebos[[i]],"Plac p","placebos",i)
  check_combined_input(synths[[i]],"Synth p","synths",i)
  if(estimator(placebos[[i]])!=estimator(synths[[i]]))
    stop("Placebo ",i," uses ",estimator(placebos[[i]])," but synth ",i," uses ",
         estimator(synths[[i]]),".", call.=FALSE)
}
check_same_regression(c(placebos,synths))

#####HC row first, then the chosen row of each placebo with the synth p value of the same row.
hc=combined_row(placebos[[1]],synths[[1]],"HC",NA)
rows=lapply(seq_along(placebos),function(i)
  combined_row(placebos[[i]],synths[[i]],estimator(placebos[[i]]),select[i]))
rr=do.call(rbind,c(list(hc),rows))
rownames(rr)=NULL
rr=rr |>
  dplyr::rename(Adj=SE) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric)&!dplyr::any_of("CI Width"),~{round(.x,3)}))   #CI widths keep their significant digits
is_scpc=vapply(placebos,estimator,"")=="SCPC"
names(rr)[2]=if(all(is_scpc)) "avc" else if(any(is_scpc)) "Clusters/avc" else "Clusters"

sp=placebos[[1]]$Spatial_Params|>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric),~{round(.x,2)}))
ss=synths[[1]]$Spatial_Params|>
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

#####Estimator of placebo or synth results: the SE label of their non-HC rows (BCH, IM or SCPC).
estimator=function(x) setdiff(x$Results$SE,"HC")[1]

#####Placebo and synth results (not IV) with the expected simulated p value column.
check_combined_input=function(x,p_col,arg,i){
  ok=is.list(x) && is.data.frame(x$Results) && p_col %in% names(x$Results) &&
    !is.null(x$Spatial_Params)
  if(!ok)
    stop("Element ",i," of `",arg,"` is not the output of ",
         if(arg=="placebos") "placebo(), placebo_im() or placebo_scpc()"
         else "synth(), synth_im() or synth_scpc()",".", call.=FALSE)
  if(is_iv_results(x$Results))
    stop("combined_table() does not take iv_placebo() or iv_synth() results; use placebo_table() and synth_table().",
         call.=FALSE)
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

#####Row of a placebo whose SE label is `se` and whose clusters or avc equal `value` (any value for
#####HC), with the synth p value of the matching synth row.
combined_row=function(plbo,syn,se,value){
  pick=function(results){
    level=results[[2]]   #Clusters or avc, stored as character with "." for HC
    rows=results$SE==se
    if(se!="HC") rows=rows & (abs(suppressWarnings(as.numeric(level))-value)<1e-8) %in% TRUE
    rows=which(rows %in% TRUE)
    if(length(rows)==0){
      available=paste(level[results$SE==se],collapse=", ")
      stop("No ",se," row with ",names(results)[2]," = ",value,". Available values: ",available,".",
           call.=FALSE)
    }
    results[rows[1],,drop=FALSE]
  }
  row=pick(plbo$Results)
  row$`Synth p`=pick(syn$Results)$`Synth p`
  names(row)[2]="Level"
  row[,c("SE","Level","Est p","Plac p","Plac 5%","Synth p","CI Width","CI")]
}
