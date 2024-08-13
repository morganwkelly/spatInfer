#' Placebo significance level.
#'@description
#'
#'  `plot_basis()` gives a 3d plot of the fitted tensor spline surface of the outcome. This is just a 
#'  wrapper around the mgcv command vis.gam()
#'  
#' @param fm  Regression formula
#' @param df Name of dataframe
#' @param splines Number of tensor splines given by optimal_basis()
#' @param grd_num Number of grid lines in plot.
#' @param phi   Figure viewing height in degrees
#' @param theta Figure rotation in degrees.
#' @param Title String giving title of plot
#'
#' @return A 3d plot of the fitted outcome surface.#' 
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' # Use 100 observations and 100 simulations to speed things up.
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=100)

#' # Use the number of splines indicated by [optimal_basis()].
#' plot_basis(mobility~racial_seg+single_mom,  opportunity,
#' splines=4)
#' 


plot_basis=function(fm,df,splines,grd_num=20,phi=50,theta=30, Title=""){

  
  new_names=set_names(fm)
  df=df |> dplyr::rename(dep_var=new_names$gg[1],
                         explan_var=stringr::str_split_1(new_names$gg[2],"\\+")[1])
  rhs=new_names$rhs
  gm_2=mgcv::bam(dep_var~
                   te(X,Y,bs=c("bs"),
                      k=splines,
                      m=1),
                 data=df,
                 discrete=T)
  
  mgcv::vis.gam(gm_2,n.grid=grd_num,theta=30,phi=50,ticktype="detailed", #   color="gray",   #phi does view height, theta rotation
          xlab="Long",ylab="Lat",zlab="Outcome",main=Title)
}
  
#plot_basis(fm,df,splines=5)
