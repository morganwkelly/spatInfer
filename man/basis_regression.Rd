% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_regression.R
\name{basis_regression}
\alias{basis_regression}
\title{Carry out spatial basis regressions with BCH standard errors.}
\usage{
basis_regression(fm, df, splines, pc_num, clusters, weights = F, cov = "BCH")
}
\arguments{
\item{fm}{The equation formula. The variable of interest is the first one on the right hand side.}

\item{df}{The name of the dataset use. Longitude and latitude must be named X and Y and have no missing values.}

\item{splines}{The dimension of the linear tensor used, from optimal_basis function.}

\item{pc_num}{The number of principal component to include, again from optimal_basis.}

\item{clusters}{Number of k-medoids clusters to use based on the placebo test.}

\item{weights}{Set weights=T if the regression is weighted. The weighting variable in the dataset must be named weights.}

\item{cov}{Defaults to BCH. It gives heteroskedasticity consistent standard errors otherwise.}
}
\value{
feols object that can be printed and exported using the modelsummary package.
Given that t-statistics from large cluster regressions are not immediately interpretable due to their low degrees
of freedom the regression table gives a confidence interval and p-value, although this is trivial to change in
the modelsummary output.
}
\description{
Carry out spatial basis regressions with BCH standard errors.
}
\examples{
library(spatInfer)
data(opportunity)
# Choose a sample of 250 observations and run only 100 placebo simulations
# to speed things up.
set.seed(123)
opportunity=opportunity |> dplyr::slice_sample(n=250)
# Use the number of splines and PCs indicated by optimal_basis and
# the number of clusters from the placebo test.
ck=basis_regression(mobility~racial_seg+single_mom,  opportunity,
splines=4,pc_num=3,
clusters=5,cov="BCH")

#Now generate a regression output table.
modelsummary::modelsummary(list(CK=ck),
statistic = c("conf.int","p = {p.value}"),
coef_omit = c("Intercept|PC*"), #omit basis and intercept
gof_map = c("nobs", "r.squared"),fmt=2)
}
