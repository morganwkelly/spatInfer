#' Synthetic outcome tests for instrumental variable regressions.
#'@description
#'
#'  `iv_synth()` tests the null hypothesis that the outcome of an instrumental
#'  variable regression is spatial noise, and therefore unrelated to the
#'  treatment. The outcome is replaced by synthetic outcomes: a quadratic trend
#'  in longitude and latitude, plus noise with the same spatial correlation as
#'  the residuals from that trend. It gives the estimated and synthetic outcome
#'  significance levels of the treatment for 2SLS and for the Anderson-Rubin
#'  (AR) test, each with HC, BCH, IM and (optionally) SCPC inference, for each
#'  number of large clusters or SCPC bound `avc`.
#'
#'@inheritParams iv_placebo
#'@inheritParams synth
#'@param nSim The number of synthetic outcomes to generate. Defaults to 1000
#'  but lower values should be used first to get an idea of how the regression
#'  is behaving.
#'@param max_clus The maximum number of large clusters to examine. Defaults to
#'  6.
#'@details See [iv_placebo()] for the inference methods and AR confidence sets.
#'
#'@return  A list containing `Results`, which gives the estimated and
#'  synthetic outcome p values and the first stage F statistic for each kind of
#'  inference and number of clusters (except IM and SCPC), and `Spatial_Params`, giving the Moran
#'  statistic of the reduced form residuals, the R2 of the outcome on a
#'  quadratic in longitude and latitude, and the effective range and structure
#'  of the synthetic noise. Show the results with [synth_table()], which gives
#'  the 2SLS and AR results in one table.
#'@export
#'
#' @examples
#' \donttest{
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' # Use 100 simulations and turn off parallel processing.
#' syn=iv_synth(d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag,
#' china,
#' splines=6,
#' pc_num=8,
#' nSim=100,
#' weights="timepwt48",
#' Parallel=FALSE
#' )
#' synth_table(syn)
#' }

iv_synth=function(fm,df,splines,pc_num,
                  nSim=1000,weights=NULL,max_clus=6,
                  Parallel=TRUE,exact_cholesky=TRUE,k_medoids=TRUE,xy_interaction=FALSE,avc=NULL){
#
  withr::local_preserve_seed()  #fixed simulation seeds do not change the user's random numbers
  if(!is.null(avc)) avc=check_avc(avc)
  prep=prepare_iv_data(fm,df,splines,pc_num,weights,max_clus)
  df=prep$df
  rhs=prep$rhs
  pc=prep$pc
  rm(prep)

# 2SLS, reduced form and first stage equations, and synthetic outcome equations where `sim`
# replaces the outcome. All have the spatial basis added.
eqs=build_iv_formulas(rhs,pc,sim_replaces="dep_var")
rm(pc)

# Regress outcome on quadratic in longitude and latitude to generate synthetic outcomes.
lm_res=synth_trend(df,xy_interaction)

# Synthetic outcomes: trend plus noise with the spatial correlation of the residuals.
noise_sim=Noise_Sim(df,lm_res,nSim,exact_cholesky,Parallel)
Sim=noise_sim$Sim

Spatial_Params=data.frame(Moran=moran(eqs$eq_rf,df),   #Moran statistic of reduced form residuals
                          R2=summary(lm_res)$r.squared,   #explanatory power of quadratic trend for outcome
                          Effective_Range=noise_sim$matern_params$Effective_Range,
                          Structure=noise_sim$matern_params$Structure)
Spatial_Params=round(Spatial_Params,3)
Spatial_Params=cbind.data.frame(Spatial_Params,N=nrow(df),Splines=splines,PCs=pc_num)
rm(noise_sim)

# k-medoids clusters from 3 up to max_clus
hold_clus=sim_clusters(df,k_medoids,max_clus)

# SCPC spatial weights for each avc, which depend only on the locations.
spcs=if(is.null(avc)) NULL else scpc_setup(df,avc)

# p values of the treatment (2SLS) and of the instrument in the reduced form (AR) with synthetic outcomes.
sim_out=iv_sim_p_values(df,Sim,eqs,hold_clus,nSim,Parallel,avc=avc,spcs=spcs)

Results=iv_summary(df,eqs,hold_clus,sim_out,avc,spcs) |>
  dplyr::select(-sim_05) |>
  dplyr::rename(`Synth p`= sim_p,
                `Est p`=est_p,`CI Width`=width_ci)

obj=list(
  Results=Results,
  Spatial_Params=Spatial_Params
)
return(obj)

}
