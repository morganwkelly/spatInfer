#' Placebo significance levels using SCPC inference.
#'@description
#'
#'  `placebo_scpc()` runs a spatial noise placebo test to choose the average
#'  correlation bound `avc` for SCPC inference (Mueller and Watson) and obtain
#'  the placebo significance level of the treatment variable. The placebo
#'  matches the spatial trends of the explanatory variable by regressing it on
#'  the principal components chosen by [optimal_basis()] and then simulates
#'  noise with the same correlation structure as the residuals from that
#'  regression. For each value of `avc` it gives the regression and placebo
#'  significance levels.
#'
#' @param fm  The equation formula. The variable of interest is the first one on
#'   the right hand side.
#' @param df  The name of the dataset use. Longitude and latitude must be named
#'   X and Y and have no missing values.
#' @param splines The dimension of the linear tensor used, from optimal_basis
#'   function
#' @param pc_num The number of principal component to include, again from
#'   optimal_basis function.
#' @param nSim The number of placebos to generate. Defaults to 1000 but lower
#'   values should be used first to get an idea of how the regression is
#'   behaving.
#' @param weights Set weights=TRUE if the regression is weighted. The weighting
#'   variable in the dataset must be named weights.
#' @param avc The values of the SCPC upper bound on the average pairwise
#'   correlation to examine with the placebo test. Each must lie between 0.001
#'   and 0.99. Defaults to 0.02, 0.04, 0.06, 0.08 and 0.1.
#' @param Parallel Speed things up by running simulations in parallel. Set to
#'   `FALSE` if this creates problems.
#' @param exact_cholesky Use an exact Cholesky decomposition to generate
#'   synthetic noise. For very large datasets, setting this to `FALSE` will use
#'   the BRISC Cholesky approximation. This requires the `BRISC` package.
#' @param jitter_coords If some sites have identical coordinates, jitter them for the
#'   Moran test by adding Gaussian noise with standard deviation 0.01 degrees
#'   (about 1 km). Set to `FALSE` to use the coordinates as they are, in which
#'   case spdep warns about identical points. Only the Moran test is affected.
#'
#' @details SCPC inference uses [scpcR::scpc()] with X and Y as longitude and
#'   latitude, so distances are great-circle distances. Larger values of `avc`
#'   allow for stronger spatial correlation and give wider confidence intervals.
#'   `Pseudo SE` is one quarter of the width of the 95% confidence interval, in
#'   the units of the treatment coefficient, so that the interval is roughly the
#'   estimate plus or minus two pseudo standard errors. A simulation where SCPC
#'   cannot be computed is left out, with a warning.
#'
#' @return A list containing Results which summarizes the placebo values,
#'   Spatial_Params giving the Moran test value and the range and structure used
#'   to generate the placebos, and Estimates giving the SCPC estimate, standard
#'   error, t statistic, p value and 95% confidence interval of the treatment
#'   for each value of `avc`. Choose the smallest `avc` where the proportion of
#'   placebo regressions significant at 5% is close to 0.05.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' # Use 100 observations and 100 simulations to speed things up.
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=100)
#' # Use the number of splines and PCs indicated by `optimal_basis()`. Turn off parallel processing.
#' plbo_scpc=placebo_scpc(mobility~single_mothers+short_commute+gini+dropout_rate+social_cap,
#' opportunity,
#' splines=6,
#' pc_num=15,
#' nSim=100,
#' Parallel=FALSE
#' )
#' placebo_table(plbo_scpc, caption="SCPC placebo values for single mothers variable.")
#'


placebo_scpc=function(fm,df,splines,pc_num,
                      nSim=1000,weights=FALSE,avc=seq(0.02,0.1,by=0.02),
                      Parallel=TRUE,exact_cholesky=TRUE,jitter_coords=TRUE){
#
  withr::local_preserve_seed()  #fixed simulation seeds do not change the user's random numbers
  avc=check_avc(avc)
  prep=prepare_spatial_data(fm,df,splines,pc_num,weights)
  df=prep$df
  rhs=prep$rhs
  pc=prep$pc
  rm(prep)


# Define two equations: the estimated equation using the true explanatory variable and
# the simulated equation that uses a placebo instead.
# Both equations have a spatial basis of principal components added.
eqs=build_formulas(rhs,pc,sim_replaces="explan_var")
eq_est=eqs$eq_est
eq_sim=eqs$eq_sim
rm(eqs)

# Regress explanatory variable on spatial basis to get trend of placebo and residuals.
lm_res=lm(as.formula(paste("explan_var",paste(names(pc),collapse="+"),sep="~")),
          data=df)

rm(pc)

# Get spatial correlation pattern of residuals and generate placebo noise with the same structure.
noise_sim=Noise_Sim(df,lm_res,nSim,exact_cholesky,Parallel)
#Collect Output

#Moran test using 5 nearest neighbours
Moran=moran(eq_est,df,jitter_coords=jitter_coords)

Spatial_Params=data.frame(Moran,R2=summary(lm_res)$r.squared,   #explanatory power of principal components for x
                          Effective_Range=noise_sim$matern_params$Effective_Range,   #fraction of 95th perc distance bw coords
                          Structure=noise_sim$matern_params$Structure)
Spatial_Params=round(Spatial_Params,3)
Spatial_Params=cbind.data.frame(Spatial_Params,Splines=splines,PCs=pc_num)
Sim=noise_sim$Sim
rm(noise_sim)

# Get p values of placebo variables in simulated regressions both for SCPC and HC standard errors.
scpc_out=scpc_p_values(Sim,eq_sim,df,avc,nSim,Parallel)
hc_out=hc_p_values(Sim,eq_sim,df,nSim,Parallel)

# Summary Placebo results
summ=summary_scpc(df,eq_est,avc,scpc_out,hc_out)
Results=summ$Results |>
                 dplyr::rename(`Plac p`= sim_p, `Plac 5%`=sim_05,
                               `Est p`=est_p,`CI Width`=width_ci,`Pseudo SE`=pseudo_se)

obj=list(
  Results=Results,
  Spatial_Params=Spatial_Params,
  Estimates=summ$Estimates
)
return(obj)

}
