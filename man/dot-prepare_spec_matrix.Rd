% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SandwichLayerVariance.R
\name{.prepare_spec_matrix}
\alias{.prepare_spec_matrix}
\title{(Internal) Helper function for design-based meat matrix
calculation}
\usage{
.prepare_spec_matrix(x, ...)
}
\arguments{
\item{x}{a fitted \code{teeMod} model}
}
\value{
a \eqn{m \times (p+2)} matrix of cluster sums of design-based
estimating equations scaled by \eqn{\sqrt{m_{b0}m_{b1}}/m_{b}}. Here
\eqn{m} is the number of clusters, \eqn{p} is the number of covariates used
in the prior covariance adjustment (excluding intercept)
}
\description{
(Internal) Helper function for design-based meat matrix
calculation
}
\keyword{internal}
