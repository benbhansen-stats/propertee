% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StudySpecificationConverters.R
\name{as_rct_spec}
\alias{as_rct_spec}
\alias{as_obs_spec}
\alias{as_rd_spec}
\title{Convert \code{StudySpecification} between types}
\usage{
as_rct_spec(StudySpecification, ..., loseforcing = FALSE)

as_obs_spec(StudySpecification, ..., loseforcing = FALSE)

as_rd_spec(StudySpecification, data, ..., forcing)
}
\arguments{
\item{StudySpecification}{a \code{StudySpecification} to convert}

\item{...}{Ignored.}

\item{loseforcing}{converting from RD to another \code{StudySpecification}
type will error to avoid losing the forcing variable. Setting
\code{loseforcing = TRUE} allows the conversion to automatically drop the
forcing variable. Default \code{FALSE}.}

\item{data}{converting to an RD requires adding a \code{forcing} variable,
which requires access to the original data.}

\item{forcing}{converting to an RD requires adding a \code{forcing}
variable. This should be entered as a formula which would be passed to
\code{update}, e.g. \code{forcing = . ~ . + forcing(forcevar)}.}
}
\value{
\code{StudySpecification} of the updated type
}
\description{
Convert a \code{StudySpecification} between a observational
study, a randomized control trial, and a regression discontinuity (created
from \code{obs_spec}, \code{rct_spec} and \code{rd_spec} respectively).
}
\examples{
spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
spec
as_obs_spec(spec)
as_rd_spec(spec, simdata, forcing = ~ . + forcing(force))
spec2 <- rd_spec(o ~ uoa(uoa1, uoa2) + forcing(force), data = simdata)
spec2
# as_rct_spec(spec2) # this will produce an error
as_rct_spec(spec2, loseforcing = TRUE)
}
