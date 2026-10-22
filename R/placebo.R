#' Placebo significance levels with HC, BCH, IM and SCPC inference.
#'@description
#'
#'  `placebo()` runs a spatial noise placebo test on a regression to choose the
#'  number of large clusters for BCH and IM inference, and the average
#'  correlation bound `avc` for SCPC inference, and to obtain the placebo
#'  significance level of the treatment variable. The placebo matches
#'  the spatial trends of the explanatory variable by regressing it on the
#'  principal components chosen by [optimal_basis()] and then simulates noise
#'  with the same correlation structure as the residuals from that regression.
#'  For HC inference, for each number of clusters with BCH and IM, and for each
#'  value of `avc` with SCPC, it gives the regression and placebo significance
#'  levels.
#'
#'@import stats
#'
#'@param fm  The equation formula. The variable of interest should be the first
#'  one on the right hand side.
#'@param df  The name of the dataset use. Longitude and latitude must be named
#'  `X` and `Y` and have no missing values.
#'@param splines The dimension of the linear tensor used, from [optimal_basis()]
#'  function.
#'@param pc_num The number of principal component to include, again from
#'  optimal_basis function.
#'@param nSim The number of placebos to generate. Defaults to 1000 but lower
#'  values should be used first to get an idea of how the regression is
#'  behaving.
#'@param weights The name of the weighting variable in the dataset, as a
#'  string, for example `weights="pop"`. Leave it out for an unweighted
#'  regression.
#'@param max_clus The maximum number of large clusters to examine with the
#'  placebo test, for BCH and IM. Defaults to 6.
#'@param Parallel Run the simulations in parallel. Set to `FALSE` if there are memory problems,
#'  which can happen especially with large datasets.
#'@param exact_cholesky For very large datasets, setting this to `FALSE` will use
#'  the `BRISC` Cholesky approximation. This requires the
#'  `BRISC` package.
#'@param k_medoids For large datasets, set to `FALSE` to use a faster approximation to generate
#'  medoids.
#'@param avc The values of the SCPC upper bound on the average pairwise
#'  correlation to examine, each between 0.001 and 0.99, for example
#'  `avc=seq(0.02,0.1,by=0.02)`. The default, `NULL`, gives no SCPC results,
#'  which are the slowest to compute.
#'@details BCH uses clustered standard errors. In IM inference the treatment
#'  coefficient is estimated separately in each cluster and the cluster
#'  coefficients are tested with a t-test. Clusters where it cannot be
#'  estimated, for example because the treatment does not vary within the
#'  cluster, are left out of the t-test, with a warning when this happens in
#'  the regression on the actual data.
#'
#'  The SCPC code is based on the scpcR package of Sascha Becker, P. David Boll
#'  and Hans-Joachim Voth (<https://github.com/spatial-spur/scpcR>, MIT licence),
#'  ported into spatInfer so that scpcR is not needed. It gives the same results
#'  as scpcR for unweighted regressions. With `weights`, conditional SCPC is
#'  computed for the regression rescaled by the square root of the weights, the
#'  least squares regression equivalent to weighted least squares, whereas scpcR
#'  leaves the regressors unweighted.
#'
#'  With 4,500 or more observations, SCPC uses scpcR's large sample approximation, which
#'  computes the spatial weights from random subsets of 1,000 locations and extends them
#'  to the rest. As in scpcR, this extension measures the distance between locations in
#'  degrees of longitude and latitude, whereas the rate at which correlation decays with
#'  distance is calibrated on great-circle distances in fractions of the earth's
#'  circumference. The weights therefore decay faster with distance than intended, and
#'  results can differ from the exact calculation used for smaller samples. The
#'  calculation is left unchanged so that results match scpcR.
#'  SCPC inference uses X and Y as longitude and
#'  latitude, so distances are great-circle distances. Larger values of `avc`
#'  allow for stronger spatial correlation and give wider confidence intervals.
#'  A simulation where SCPC cannot be computed is left out, with a warning.
#'
#'@return  A list containing `Results` which summarizes the placebo values and
#'  `Spatial_Params` giving the Moran test value and the range and structure
#'  used to generate the placebos. `Results` has an HC row, then BCH and IM rows
#'  for 3 up to `max_clus` clusters, then (when `avc` is given) an SCPC row for
#'  each value of `avc`, which `Clusters` holds. For BCH and IM choose the number of clusters where the
#'  proportion of placebo regressions significant at 5% is in the region of
#'  0.05 to 0.07; given the conservatism of IM, the highest proportion is
#'  sometimes only 0.03--0.04. For SCPC choose the smallest `avc` where it is
#'  close to 0.05. Show the results with [placebo_table()].
#'@export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' # Use 100 observations and 100 simulations to speed things up.
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=100)
#' # Use the number of splines and PCs indicated by [optimal_basis()] and turn off parallel processing.
#' plbo=placebo(mobility~single_mothers+short_commute+gini+dropout_rate+social_cap,opportunity,
#' splines=6,
#' pc_num=15,
#' nSim=100,
#' Parallel=FALSE,
#' max_clus = 7,
#' avc=c(0.02,0.04,0.06)
#' )
#' placebo_table(plbo)


placebo=function(fm,df,splines,pc_num,
                     nSim=1000,weights=NULL,max_clus=6,
                     Parallel=TRUE,exact_cholesky=TRUE,
                 k_medoids=TRUE,avc=NULL){
#
  withr::local_preserve_seed()  #fixed simulation seeds do not change the user's random numbers
  if(!is.null(avc)) avc=check_avc(avc)
  prep=prepare_spatial_data(fm,df,splines,pc_num,weights,max_clus)
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
Moran=moran(eq_est,df)

Spatial_Params=data.frame(Moran,R2=summary(lm_res)$r.squared,   #explanatory power of principal components for x
                          Effective_Range=noise_sim$matern_params$Effective_Range,   #fraction of 95th perc distance bw coords
                          Structure=noise_sim$matern_params$Structure)
Spatial_Params=round(Spatial_Params,3)
Spatial_Params=cbind.data.frame(Spatial_Params,N=nrow(df),Splines=splines,PCs=pc_num)
Sim=noise_sim$Sim
rm(noise_sim)

# k-medoids clusters from 3 up to max_clus
hold_clus=sim_clusters(df,k_medoids,max_clus)

# SCPC spatial weights for each avc, which depend only on the locations.
spcs=if(is.null(avc)) NULL else scpc_setup(df,avc)

# p values of placebo variables in simulated regressions with HC, BCH, IM and SCPC inference.
sim_out=ols_sim_p_values(df,Sim,eq_sim,hold_clus,nSim,Parallel,avc,spcs)

# Summary Placebo results
Results=ols_summary(df,eq_est,hold_clus,sim_out,avc,spcs)
Results=Results  |>
                 dplyr::rename(`Plac p`= sim_p, `Plac 5%`=sim_05,
                               `Est p`=est_p,`CI Width`=width_ci)

obj=list(
  Results=Results,
  Spatial_Params=Spatial_Params

)
return(obj)

}
