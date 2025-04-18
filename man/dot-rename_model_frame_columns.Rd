% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StudySpecification.R
\name{.rename_model_frame_columns}
\alias{.rename_model_frame_columns}
\title{(Internal) Rename columns to strip function calls}
\usage{
.rename_model_frame_columns(modframe)
}
\arguments{
\item{modframe}{A \code{data.frame}.}
}
\value{
The \code{data.frame} with function calls removed
}
\description{
After calling \code{model.frame()} on the formula input to
\code{.new_StudySpecification()}, the names of the columns will include
function names, e.g. "block(blockvar)". This function strips all these.
}
\keyword{internal}
