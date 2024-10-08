% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synth_table.R
\name{synth_table}
\alias{synth_table}
\title{Makes a summary table of estimated and synthetic outcome p values from output generated by BHC or IM simulations.}
\usage{
synth_table(syn)
}
\arguments{
\item{syn}{An object containing simulation results created by synth or synth_im.
#'}
}
\value{
Generates a tinytable object. This gives the regression and synthetic outcome p values as the number of clusters increases.
The regression confidence interval and its width are also reported. The next part of the table gives
the z-score of a Moran test, using 5 nearest neighbours, the R2 of a regression of the outcome on a
quadratic in longitude and latitude, and the effective range and structure of this regression's residuals.
}
\description{
Makes a summary table of estimated and synthetic outcome p values from output generated by BHC or IM simulations.
}
\examples{
library(spatInfer)
data(opportunity)
# Choose a sample of 250 observations and run only 100 placebo simulations to speed things up.
set.seed(123)
opportunity=opportunity |> dplyr::slice_sample(n=250)
# Use the number of splines and PCs indicated by optimal_basis()
syn=synth(mobility~racial_seg+single_mom,  opportunity,
splines=4, pc_num=3,nSim=100)

synth_table(syn,adjust="BCH")
}
