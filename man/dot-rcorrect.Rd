% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teeMod.R
\name{.rcorrect}
\alias{.rcorrect}
\title{Bias correct residuals contributing to standard errors of a \code{teeMod}}
\usage{
.rcorrect(resids, x, model, type, ...)
}
\arguments{
\item{resids}{numeric vector of residuals to correct}

\item{x}{teeMod object}

\item{model}{string indicating which model the residuals are from. \code{"itt"}
indicates correction to the residuals of \code{x}, and \code{"cov_adj"}
indicates correction to the residuals of the covariance adjustment model.
This informs whether corrections should use information from \code{x} or the
\code{fitted_covariance_model} slot of the \code{SandwichLayer} object in the
\code{offset} for corrections}

\item{type}{string indicating the desired bias correction. Can be one of
\code{"(HC/CR/MB)0"}, \code{"(HC/CR/MB)1"}, or \code{"(HC/CR/MB)2"}}

\item{...}{additional arguments passed from up the call stack; in particular,
the \code{cluster_cols} argument, which informs whether to cluster and provide
CR2 corrections instead of HC2 corrections, as well as the correction for the
number of clusters in the CR1 correction. This may also include a \code{by}
argument.}
}
\description{
Bias correct residuals contributing to standard errors of a \code{teeMod}
}
\keyword{internal}
