% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StudySpecificationAccessors.R
\name{.specification_accessors_newdata_validate}
\alias{.specification_accessors_newdata_validate}
\title{(Internal) Checks newdata/by argument for specification accessors}
\usage{
.specification_accessors_newdata_validate(newdata, by)
}
\arguments{
\item{newdata}{newdata argument from e.g. \code{treatment()},
\code{blocks()}, etc}

\item{by}{from e.g. \code{treatment()}, \code{blocks()}, etc. See
\code{.check_by()}}
}
\value{
Invisibly \code{TRUE}. Warns or errors as appropriate.
}
\description{
(Internal) Checks newdata/by argument for specification accessors
}
\keyword{internal}
