% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StudySpecificationConverters.R
\name{.remove_forcing}
\alias{.remove_forcing}
\title{(Internal) Removes the forcing column entirely from a
\code{StudySpecification}}
\usage{
.remove_forcing(spec)
}
\arguments{
\item{spec}{A \code{StudySpecification}}
}
\value{
The \code{StudySpecification} without any forcing variable
}
\description{
In preparation for converting an RD \code{StudySpecification} to another
\code{StudySpecification}, this will strip the forcing variable entirely. It
is removed from the data (both \code{@structure} and \code{@column_index}),
as well as from the formula stored in \code{@call}.
}
\details{
Note that the output \code{StudySpecification} will fail a validity check
(with \code{validObject()}) due to an RD \code{StudySpecification} requiring
a forcing variable, so change the \code{@type} immediately.
}
\keyword{internal}
