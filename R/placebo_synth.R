#' Placebo and synthetic outcome tests in one command.
#'
#'@description
#' `placebo_synth()` runs the placebo test of [placebo()] and the synthetic
#' outcome test of [synth()] on the same regression, with a progress bar
#' showing the step being run and the number of simulated regressions done.
#' The data, clusters and SCPC spatial weights are prepared once and shared by
#' both tests, and the results are the same as those of separate calls to
#' [placebo()] and [synth()].
#'
#' @inheritParams placebo
#' @inheritParams synth
#' @param nSim The number of placebos, and of synthetic outcomes, to generate.
#'   Defaults to 1000 but lower values should be used first to get an idea of
#'   how the regression is behaving.
#' @param progress Show a progress bar. Defaults to `TRUE`.
#' @details See [placebo()] for the inference methods.
#' @return A list with `Placebo`, the output of [placebo()], and `Synth`, the
#'   output of [synth()]. Show them in one table with [placebo_synth_table()].
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' # Use 100 observations and 100 simulations to speed things up. Turn off parallel processing.
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=100)
#' ps=placebo_synth(mobility~single_mothers+short_commute+gini+dropout_rate+social_cap,
#' opportunity,
#' splines=6,
#' pc_num=15,
#' nSim=100,
#' Parallel=FALSE,
#' max_clus = 7,
#' avc=c(0.02,0.04,0.06)
#' )
#' placebo_synth_table(ps)
#' placebo_synth_table(ps, select=c(BCH=6, IM=6, SCPC=0.04))

placebo_synth=function(fm,df,splines,pc_num,
                       nSim=1000,weights=NULL,max_clus=6,
                       Parallel=TRUE,exact_cholesky=TRUE,k_medoids=TRUE,xy_interaction=FALSE,avc=NULL,
                       scpc_weighting="scpcR",progress=TRUE){
#
  withr::local_preserve_seed()  #fixed simulation seeds do not change the user's random numbers
  inp=sim_inputs(fm,df,splines,pc_num,weights,max_clus,k_medoids,avc,scpc_weighting)
  pb=if(isTRUE(progress)) sim_progress(2*nSim) else NULL
  obj=list(
    Placebo=placebo_run(inp,nSim,Parallel,exact_cholesky,pb),
    Synth=synth_run(inp,nSim,Parallel,exact_cholesky,xy_interaction,pb)
  )
  return(obj)
}

#' Placebo and synthetic outcome tests for instrumental variable regressions in one command.
#'
#'@description
#' `iv_placebo_synth()` runs the placebo test of [iv_placebo()] and the
#' synthetic outcome test of [iv_synth()] on the same instrumental variable
#' regression, with a progress bar showing the step being run and the number of
#' simulated regressions done. The data, clusters and SCPC spatial weights are
#' prepared once and shared by both tests, and the results are the same as
#' those of separate calls to [iv_placebo()] and [iv_synth()].
#'
#' @inheritParams iv_placebo
#' @inheritParams iv_synth
#' @inheritParams placebo_synth
#' @details See [iv_placebo()] for the inference methods and Anderson-Rubin
#'   confidence sets.
#' @return A list with `Placebo`, the output of [iv_placebo()], and `Synth`,
#'   the output of [iv_synth()]. Show them in one table with
#'   [placebo_synth_table()].
#' @export
#'
#' @examples
#' \donttest{
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' # Use 100 simulations and turn off parallel processing.
#' ps=iv_placebo_synth(d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag,
#' china,
#' splines=6,
#' pc_num=8,
#' nSim=100,
#' weights="timepwt48",
#' Parallel=FALSE
#' )
#' placebo_synth_table(ps)
#' }

iv_placebo_synth=function(fm,df,splines,pc_num,
                          nSim=1000,weights=NULL,max_clus=6,
                          Parallel=TRUE,exact_cholesky=TRUE,k_medoids=TRUE,xy_interaction=FALSE,avc=NULL,
                          scpc_weighting="scpcR",progress=TRUE){
#
  withr::local_preserve_seed()  #fixed simulation seeds do not change the user's random numbers
  inp=iv_sim_inputs(fm,df,splines,pc_num,weights,max_clus,k_medoids,avc,scpc_weighting)
  pb=if(isTRUE(progress)) sim_progress(2*nSim) else NULL
  obj=list(
    Placebo=iv_placebo_run(inp,nSim,Parallel,exact_cholesky,pb),
    Synth=iv_synth_run(inp,nSim,Parallel,exact_cholesky,xy_interaction,pb)
  )
  return(obj)
}

#' Table of placebo and synthetic outcome results.
#'
#'@description
#' `placebo_synth_table()` makes one summary table from the output of [placebo_synth()] or
#' [iv_placebo_synth()]. By default it shows every inference method, number of clusters and value
#' of `avc`, as [placebo_table()] does, with the synthetic outcome p value of each row added. With
#' `select` it shows only the HC row and the chosen rows, as [combined_table()] and
#' [iv_combined_table()] do.
#'
#' @param ps An object created by [placebo_synth()] or [iv_placebo_synth()].
#' @param select `NULL` (the default) for every row, or a named numeric vector giving, for each
#'   inference method to show, the number of clusters (`BCH` and `IM`) or the value of `avc`
#'   (`SCPC`) of its row, for example `c(BCH=7, IM=5, SCPC=0.04)`.
#' @param ar_select For [iv_placebo_synth()] results with `select`, the rows chosen for the
#'   Anderson-Rubin tests, as in [iv_combined_table()]. Defaults to `select`.
#' @param caption An optional caption for the table.
#'
#' @return Generates a tinytable object giving, for each row, the estimated and placebo p values,
#'   the proportion of placebo regressions significant at 5%, the synthetic outcome p value and the
#'   regression confidence interval and its width. For IV results the table has a block of 2SLS rows
#'   followed by a block of Anderson-Rubin rows, with first stage F statistics. The notes give the
#'   z-score of a Moran test, and the structure, effective range and R2 used to simulate noise for
#'   the treatment and the outcome (and for IV placebos the instrument).
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' # Use 100 observations and 100 simulations to speed things up. Turn off parallel processing.
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=100)
#' ps=placebo_synth(mobility~single_mothers+short_commute+gini+dropout_rate+social_cap,
#' opportunity, splines=6, pc_num=15, nSim=100, Parallel=FALSE, max_clus=7, avc=c(0.02,0.04))
#' placebo_synth_table(ps, caption="Placebo and synthetic outcome significance levels.")
#' placebo_synth_table(ps, select=c(BCH=6, IM=6, SCPC=0.04))

placebo_synth_table=function(ps, select=NULL, ar_select=select, caption=""){
if(!is.list(ps) || !all(c("Placebo","Synth") %in% names(ps)))
  stop("`ps` is not the output of placebo_synth() or iv_placebo_synth().",call.=FALSE)
is_iv=is.data.frame(ps$Placebo$Results) && is_iv_results(ps$Placebo$Results)
fun=if(is_iv) "iv_placebo_synth()" else "placebo_synth()"
check_combined_input(ps$Placebo,"Plac p","ps$Placebo",fun,iv=is_iv)
check_combined_input(ps$Synth,"Synth p","ps$Synth",fun,iv=is_iv)
if(!is.null(select)){
  if(is_iv) return(iv_combined_table(ps$Placebo,ps$Synth,select,ar_select,caption))
  return(combined_table(ps$Placebo,ps$Synth,select,caption))
}
return(sim_table(ps$Placebo,"placebo",caption,syn=ps$Synth))
}
