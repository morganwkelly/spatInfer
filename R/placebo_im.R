#' IM placebo tests
#'
#' Run placebo tests on regressions to choose the optimal number of large
#' clusters for IM inference and obtain the placebo significance level of the
#' treatment variable.
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
#' @param weights Set weights=T if the regression is weighted. The weighting
#'   variable in the dataset must be named weights.
#' @param max_clus The maximum number of clusters to examine with the placebo
#'   test. Defaults to 6.
#' @param Parallel Speed things up by running simulations in parallel. Set to F
#'   if this creates problems.
#' @param exact_cholesky Use an exact Cholesky decomposition to generate
#'   synthetic noise. For very large datasets, setting this to F will use the
#'   BRISC Cholesky approximation.
#' @param k_medoids Use k-medoids clustering (PAM). For large datasets, set to F
#'   to use Clara to generate medoids.
#' @param jitter_coords If some sites have identical coordinates, jitter by
#'   adding Gaussian noise with standard deviation of 0.01 (10 km) to allow
#'   Moran test to be calculated.
#'
#' @return A list containing Results which summarizes the placebo values and
#'   Spatial_Params giving the Moran test value and the range and structure used
#'   to generate the placebos. Choose the number of clusters where the
#'   proportion of placebo regressions significant at 5\% is in the region of
#'   0.05 to 0.07. Sometimes, however, given the conservatism of IM, the highest
#'   proportion will only be 0.03--0.04.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' # Use 100 observations and 100 simulations to speed things up
#'
#'set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=100)
#' # Use the number of splines and PCs indicated by optimal_basis()
#' plbo_im=placebo_im(mobility~racial_seg+single_mom,  opportunity,
#' splines=4, pc_num=3,nSim=100,Parallel=FALSE)
#'
#' placebo_table(plbo_im,adjust="IM")
#'


placebo_im=function(fm,df,splines,pc_num,
                    nSim=1000,weights=FALSE,max_clus=6,
                    Parallel=FALSE,exact_cholesky=TRUE,
                    k_medoids=TRUE,jitter_coords=TRUE){
#
  if(is.null(df$X)|is.null(df$Y))
    stop("You must have longitude and latitude variables named X and Y")
  if(sum(is.na(df$X))>0|sum(is.na(df$Y))>0)
    stop("You cannot have missing values in longitude and latitude.")
  if(max_clus<3)
    stop("Your maximum number of clusters max_clus must be greater than 2.")

  if(!weights){
    df$wts=1
  }else{
    if(is.null(df$weights)){
      stop("There is no variable called weights in your data.")
    }else{
    df$wts=df$weights}
  }

#rename dependent and explanatory variables as dep_var and explan_var and list all other variables in a string called rhs
  new_names=set_names(fm)
  df=df |> dplyr::rename(dep_var=new_names$gg[1],
                   explan_var=stringr::str_split_1(new_names$gg[2],"\\+")[1])
  rhs=new_names$rhs

#get the principal components that minimise BIC and add them to the dataset.
pc=prin_comp(df,splines,pc_num)
df=cbind.data.frame(df,pc)


# Define two equations: the estimated equation using the true explanatory variable and
# the simulated equation that uses a placebo instead.
# Both equations have a spatial basis of principal components added.
eq_est=paste0("dep_var~explan_var",rhs)
eq_sim=paste0("dep_var~sim",rhs)
eq_sim=as.formula(paste(eq_sim,paste(names(pc),collapse="+"),sep="+"))
eq_est=as.formula(paste(eq_est,paste(names(pc),collapse="+"),sep="+"))

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
Spatial_Params=cbind.data.frame(Spatial_Params,Splines=splines,PCs=pc_num)
Sim=noise_sim$Sim
rm(noise_sim)

#generate k-medoids clusters of different sizes up to max_clus
hold_clus=generate_clusters(df,k_medoids,max_clus)

# Get p values of placebo variables in simulated regressions both for BCH and HC standard errors.
im_out=im_p_values(Sim,eq_sim,df,hold_clus,nSim,max_clus,Parallel)
hc_out=hc_p_values(Sim,eq_sim,df,nSim,Parallel)

# Summary Placebo results
Results=summary_im(df,eq_est,hold_clus,max_clus,im_out,hc_out)
Results=Results  |>
                 dplyr::rename(`Plac p`= sim_p, `Plac 5%`=sim_05,
                               `Est p`=est_p,`CI Width`=width_ci)

# #Collect Output
# Spatial_Params=data.frame(R2=summary(lm_res)$r.squared,   #explanatory power of principal components for x
#                           Effective_Range=noise_sim$matern_params$Effective_Range,   #fraction of 95th perc distance bw coords
#                           Structure=noise_sim$matern_params$Structure)
# Spatial_Params=round(Spatial_Params,3)

obj=list(
  Results=Results,
  Spatial_Params=Spatial_Params

)
#class(obj)=c("placebo_im")
return(obj)

}


#' #'@export
#' print.placebo_im=function(obj){
#'
#'   cat("Regression and placebo significance levels for BCH and HC standard errors.\nPlacebo_05 gives the percentage of placebo simulations significant at 5%.\nConf Int gives the width of the 95% confidence interval.\n")
#'   results=obj$Results #|>
#'    # dplyr::rename(`Placebo p`= sim_p, `Placebo 5%`=sim_05,`Est p`=est_p,`CI Width`=width_ci,`Var Ratio`=var_range)
#'   print(results)
#'   cat("\n\nSpatial parameters of treatment variable.\nSpatial R2 is the R2 of the regression of the treatment on the spatial basis.\n\n")
#'   print(obj$Spatial_Params)
#'   }

