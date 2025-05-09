% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StudySpecification.R
\name{var_table}
\alias{var_table}
\alias{var_names}
\title{Extract Variable Names from \code{StudySpecification}}
\usage{
var_table(specification, compress = TRUE, report_all = FALSE)

var_names(specification, type, implicitBlocks = FALSE)
}
\arguments{
\item{specification}{a \code{StudySpecification} object}

\item{compress}{should multiple variables be compressed into a
comma-separated string? Default \code{TRUE}. If \code{FALSE}, multiple
columns can be created instead.}

\item{report_all}{should we report all possible structures even if they
don't exist in the \code{StudySpecification}? Default \code{FALSE}.}

\item{type}{one of "t", "u", "b", "f"; for "treatment",
"unit_of_assignment", "block", and "forcing" respectively}

\item{implicitBlocks}{If the \code{StudySpecification} is created without
blocks, setting this to \code{TRUE} will return "\code{.blocks_internal}"
as the variable name corresponding to the blocks.}
}
\value{
\code{var_table} returns the requested table. \code{var_names}
returns a vector of variable names.
}
\description{
Methods to extract the variable names to the elements of the
structure of the \code{StudySpecification} (e.g. treatment, unit of
analysis, etc)
}
\details{
When \code{compress} is \code{TRUE}, the result will always have
two columns. When \code{FALSE}, the result will have number of columns
equal to the largest number of variables in a particular role, plus one.
E.g., a call such as \code{rct_spec(z ~ unitid(a, b, c, d) ...} will have
4+1=5 columns in the output matrix with \code{compress = FALSE}.

When \code{report_all} is \code{TRUE}, the matrix is guaranteed to have 3
rows (when the \code{specification} is an RCT or Obs) or 4 rows (when the
\code{specification} is a RD), with empty variable entries as appropriate.
When \code{FALSE}, the matrix will have minimum 2 rows (treatment and unit
of assignment/unitid/cluster), with additional rows for blocks and forcing
if included in the \code{StudySpecification}.
}
\examples{
spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
var_table(spec)
var_table(spec, compress = FALSE)
var_names(spec, "t")
var_names(spec, "u")
var_names(spec, "b")
}
