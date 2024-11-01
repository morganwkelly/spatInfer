% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_basis.R
\name{plot_basis}
\alias{plot_basis}
\title{Placebo significance level.}
\usage{
plot_basis(fm, df, splines, grd_num = 20, phi = 50, theta = 30, Title = "")
}
\arguments{
\item{fm}{Regression formula}

\item{df}{Name of dataframe}

\item{splines}{Number of tensor splines given by optimal_basis()}

\item{grd_num}{Number of grid lines in plot.}

\item{phi}{Figure viewing height in degrees}

\item{theta}{Figure rotation in degrees.}

\item{Title}{String giving title of plot}
}
\value{
A 3d plot of the fitted outcome surface.
}
\description{
\code{plot_basis()} gives a 3d plot of the fitted tensor spline surface of the outcome. This is just a
wrapper around the mgcv command vis.gam()
}
\examples{
library(spatInfer)
data(opportunity)
# Use the number of splines indicated by `optimal_basis()`.
plot_basis(mobility~single_mothers+short_commute+
gini+dropout_rate+social_cap+dropout_na, opportunity,
splines=6, 
Title="6x6 Tensor Surface of Mobility")

}
