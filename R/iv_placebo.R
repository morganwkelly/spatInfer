#' Placebo significance levels for instrumental variable regressions.
#'@description
#'
#'  `iv_placebo()` runs a spatial noise placebo test on an instrumental
#'  variable regression with one endogenous treatment and one instrument. It
#'  gives the estimated and placebo significance levels of the treatment for
#'  2SLS with HC, BCH and IM inference, and for the Anderson-Rubin (AR) test
#'  with HC and BCH inference, for each number of large clusters.
#'
#'  The placebo treatment matches the spatial trends of the treatment by
#'  regressing it on the principal components chosen by [optimal_basis()] and
#'  adding noise with the same correlation structure as the residuals from that
#'  regression. The placebo instrument is generated from the placebo treatment
#'  using the estimated regression of the instrument on the treatment, plus
#'  noise with the same correlation structure as the residuals from that
#'  regression. The outcome is unchanged, so placebo regressions estimate the
#'  effect of a treatment that has none, through an instrument related to it as
#'  strongly as the real one.
#'
#'@param fm  The equation formula, written as in [fixest::feols()]:
#'  `outcome ~ controls | treatment ~ instrument`, or
#'  `outcome ~ 1 | treatment ~ instrument` without controls. There must be one
#'  treatment and one instrument.
#'@param df  The dataset used. Longitude and latitude must be named `X` and `Y`,
#'  and the variables in the formula cannot have missing values.
#'@param splines The dimension of the linear tensor used, from [optimal_basis()].
#'@param pc_num The number of principal components to include, again from
#'  [optimal_basis()].
#'@param nSim The number of placebos to generate. Defaults to 1000 but lower
#'  values should be used first to get an idea of how the regression is
#'  behaving.
#'@param weights The name of the weighting variable in the dataset, as a
#'  string, for example `weights="pop"`. Leave it out for an unweighted
#'  regression.
#'@param max_clus The maximum number of large clusters to examine with the
#'  placebo test. Defaults to 6.
#'@param Parallel Run the simulations in parallel. Set to `FALSE` if there are
#'  memory problems, which can happen especially with large datasets.
#'@param exact_cholesky For very large datasets, setting this to `FALSE` will use
#'  the `BRISC` Cholesky approximation. This requires the `BRISC` package.
#'@param k_medoids For large datasets, set to `FALSE` to use a faster
#'  approximation to generate medoids.
#'@details BCH uses clustered standard errors, and IM a t-test of the
#'  coefficients estimated separately in each cluster. The AR test of a zero
#'  treatment effect is the test of the instrument in the reduced form
#'  regression of the outcome on the instrument, controls and spatial basis.
#'  Its confidence set contains the treatment effects `b` where the instrument
#'  is not significant at 5% in a regression of `outcome - b*treatment`, and it
#'  is valid however weak the instrument. When the first stage is not
#'  significant the set is unbounded: either two rays, shown as
#'  `(-Inf, a] U [b, Inf)`, or the whole line, and its width is `Inf`. IM is
#'  not valid for the AR test, so no IM version of it is given.
#'
#'  `F` in each row of the results is the first stage F statistic of the
#'  instrument with HC or BCH standard errors. It is not given (`NA`) for IM.
#'
#'@return  A list containing `Results`, which summarizes the placebo values
#'  and gives the first stage F statistic for each kind of inference and number
#'  of clusters (except IM), and `Spatial_Params`, giving the Moran statistic of the reduced
#'  form residuals, the R2 of the treatment on the principal components and the
#'  effective range and structure of the treatment noise, and the R2 of the
#'  instrument on the treatment and the effective range and structure of the
#'  instrument noise (`_Inst`). For each kind of inference choose the number of
#'  clusters where the proportion of placebo regressions significant at 5% is
#'  in the region of 0.05 to 0.07. Show the results with [placebo_table()],
#'  using `type` to choose the BCH, IM or AR table.
#'@export
#'
#' @examples
#' \donttest{
#' library(spatInfer)
#' data(autor_china)
#' china=subset(autor_china,t2==1)
#' # Use 100 simulations and turn off parallel processing.
#' plbo=iv_placebo(d_sh_empl_mfg~l_shind_manuf_cbp+l_sh_popedu_c|d_tradeusch_pw~d_tradeotch_pw_lag,
#' china,
#' splines=6,
#' pc_num=8,
#' nSim=100,
#' weights="timepwt48",
#' Parallel=FALSE
#' )
#' placebo_table(plbo)
#' }

iv_placebo=function(fm,df,splines,pc_num,
                    nSim=1000,weights=NULL,max_clus=6,
                    Parallel=TRUE,exact_cholesky=TRUE,k_medoids=TRUE){
#
  withr::local_preserve_seed()  #fixed simulation seeds do not change the user's random numbers
  prep=prepare_iv_data(fm,df,splines,pc_num,weights,max_clus)
  df=prep$df
  rhs=prep$rhs
  pc=prep$pc
  rm(prep)

# 2SLS, reduced form and first stage equations, and placebo equations where `sim` replaces the
# treatment and `sim_inst` the instrument. All have the spatial basis added.
eqs=build_iv_formulas(rhs,pc,sim_replaces="explan_var")

# Regress treatment on spatial basis to get trend of placebo and residuals.
lm_res=lm(as.formula(paste("explan_var",paste(names(pc),collapse="+"),sep="~")),
          data=df)
rm(pc)

# Placebo treatments: trend plus noise with the spatial correlation of the residuals.
noise_sim=Noise_Sim(df,lm_res,nSim,exact_cholesky,Parallel)
Sim=noise_sim$Sim

# Placebo instruments: the estimated relationship between instrument and treatment applied to the
# placebo treatments, plus noise with the spatial correlation of its residuals. The noise has its
# own seeds so that it is independent of the placebo treatment noise.
lm_inst=lm(inst_var~explan_var,data=df)
noise_inst=Noise_Sim(df,lm_inst,nSim,exact_cholesky,Parallel,seeds=c(5678,567),trend=FALSE)
Sim_inst=coef(lm_inst)[1]+coef(lm_inst)[2]*Sim+noise_inst$Sim

Spatial_Params=data.frame(Moran=moran(eqs$eq_rf,df),   #Moran statistic of reduced form residuals
                          R2=summary(lm_res)$r.squared,   #explanatory power of principal components for treatment
                          Effective_Range=noise_sim$matern_params$Effective_Range,
                          Structure=noise_sim$matern_params$Structure,
                          R2_Inst=summary(lm_inst)$r.squared,   #explanatory power of treatment for instrument
                          Effective_Range_Inst=noise_inst$matern_params$Effective_Range,
                          Structure_Inst=noise_inst$matern_params$Structure)
Spatial_Params=round(Spatial_Params,3)
Spatial_Params=cbind.data.frame(Spatial_Params,N=nrow(df),Splines=splines,PCs=pc_num)
rm(noise_sim,noise_inst)

# k-medoids clusters from 3 up to max_clus
hold_clus=iv_clusters(df,k_medoids,max_clus)

# p values of placebo treatments (2SLS) and placebo instruments (AR) in simulated regressions.
sim_out=iv_sim_p_values(df,Sim,eqs,hold_clus,nSim,Parallel,Sim_inst=Sim_inst)

Results=iv_summary(df,eqs,hold_clus,sim_out) |>
  dplyr::rename(`Plac p`= sim_p, `Plac 5%`=sim_05,
                `Est p`=est_p,`CI Width`=width_ci)

obj=list(
  Results=Results,
  Spatial_Params=Spatial_Params
)
return(obj)

}
