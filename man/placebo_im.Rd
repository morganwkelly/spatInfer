% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/placebo_im.R
\name{placebo_im}
\alias{placebo_im}
\title{IM placebo tests}
\usage{
placebo_im(
  fm,
  df,
  splines,
  pc_num,
  nSim = 1000,
  weights = FALSE,
  max_clus = 6,
  Parallel = TRUE,
  exact_cholesky = TRUE,
  k_medoids = TRUE,
  jitter_coords = TRUE
)
}
\arguments{
\item{fm}{The equation formula. The variable of interest is the first one on
the right hand side.}

\item{df}{The name of the dataset use. Longitude and latitude must be named
X and Y and have no missing values.}

\item{splines}{The dimension of the linear tensor used, from optimal_basis
function}

\item{pc_num}{The number of principal component to include, again from
optimal_basis function.}

\item{nSim}{The number of placebos to generate. Defaults to 1000 but lower
values should be used first to get an idea of how the regression is
behaving.}

\item{weights}{Set weights=T if the regression is weighted. The weighting
variable in the dataset must be named weights.}

\item{max_clus}{The maximum number of clusters to examine with the placebo
test. Defaults to 6.}

\item{Parallel}{Speed things up by running simulations in parallel. Set to F
if this creates problems.}

\item{exact_cholesky}{Use an exact Cholesky decomposition to generate
synthetic noise. For very large datasets, setting this to F will use the
BRISC Cholesky approximation.}

\item{k_medoids}{Use k-medoids clustering (PAM). For large datasets, set to F
to use Clara to generate medoids.}

\item{jitter_coords}{If some sites have identical coordinates, jitter by
adding Gaussian noise with standard deviation of 0.01 (10 km) to allow
Moran test to be calculated.}
}
\value{
A list containing Results which summarizes the placebo values and
Spatial_Params giving the Moran test value and the range and structure used
to generate the placebos. Choose the number of clusters where the
proportion of placebo regressions significant at 5\\% is in the region of
0.05 to 0.07. Sometimes, however, given the conservatism of IM, the highest
proportion will only be 0.03--0.04.
}
\description{
Run placebo tests on regressions to choose the optimal number of large
clusters for IM inference and obtain the placebo significance level of the
treatment variable.
}
\examples{
library(spatInfer)
data(opportunity)
# Use 100 observations and 100 simulations to speed things up

set.seed(123)
opportunity=opportunity |> dplyr::slice_sample(n=100)
# Use the number of splines and PCs indicated by optimal_basis()
plbo_im=placebo_im(mobility~racial_seg+single_mom,  opportunity,
splines=4, pc_num=3,nSim=100,Parallel=FALSE)

placebo_table(plbo_im,adjust="IM")

}
