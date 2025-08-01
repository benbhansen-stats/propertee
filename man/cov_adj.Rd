% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_adj.R
\name{cov_adj}
\alias{cov_adj}
\title{Covariance adjustment of \code{teeMod} model estimates}
\usage{
cov_adj(model, newdata = NULL, specification = NULL, by = NULL)
}
\arguments{
\item{model}{any model that inherits from a \code{glm}, \code{lm}, or
\code{robustbase::lmrob} object}

\item{newdata}{a dataframe of new data. Default is NULL, in which case a
dataframe is sought from higher up the call stack.}

\item{specification}{a \code{StudySpecification} object. Default is NULL, in
which case a \code{StudySpecification} object is sought from higher up the
call stack.}

\item{by}{optional; a string or named vector of unique identifier columns in
the data used to create \code{specification} and the data used to fit the
covariance adjustment model. Default is NULL, in which case unit of
assignment columns are used for identification (even if they do not
uniquely identify units of observation). If a named vector is provided,
names should represent variables in the data used to create
\code{specification}, while values should represent variables in the
covariance adjustment data.}
}
\value{
A \code{SandwichLayer} if \code{specification} is not NULL or a
\code{StudySpecification} object is found in the call stack, otherwise a
\code{PreSandwichLayer} object
}
\description{
\code{cov_adj()} takes a fitted covariance model and returns
the information necessary for adjusting direct adjustment model estimates
and associated standard errors for covariates. Standard errors will
reflect adjustments made to the outcomes as well as contributions to
sampling variability arising from the estimates of the covariance
adjustment model coefficients.
}
\details{
Prior to generating adjustments, \code{cov_adj()} identifies the
treatment variable specified in the \code{StudySpecification} object
passed to \code{specification} and replaces all values with a reference
level. If the treatment has logical type, this reference level is
\code{FALSE}, and if it has numeric type, this is the smallest
non-negative value (which means 0 for 0/1 binary). Factor treatments are
not currently supported for \code{StudySpecification} objects.\cr\cr The
values of the output vector represent adjustments for the outcomes in
\code{newdata} if \code{newdata} is provided; adjustments for the outcomes
in the data used to fit a \code{teeMod} model if \code{cov_adj()} is
called within the \code{offset} argument of the model fit; or they are the
fitted values from \code{model} if no relevant dataframe can be extracted
from the call stack. The length of the output of \code{cov_adj()} will
match the number of rows of the dataframe used.
}
\examples{
data("STARplus")

##' A prognostic model fitted to experimental + non-experimental controls
y0hat_read <- lm(read_yr1 ~ gender*dob +dobNA + race,
                 data = STARplus,
                 subset = cond_at_entry!="small")

STARspec <- rct_spec(cond_at_entry ~ unit_of_assignment(stdntid) +
                         block(grade_at_entry, school_at_entry),
                     subset=!is.na(grade_at_entry),# excludes non-experimentals
                     data = STARplus)
ett_wts    <- ett(STARspec, data = STARplus,
                  dichotomy= cond_at_entry =="small" ~.)

ett_read <- lm(read_yr1 ~ assigned(dichotomy= cond_at_entry =="small" ~.),
               offset = cov_adj(y0hat_read),
               data = STARplus,
               weights = ett_wts)
coef(ett_read)
ett_read |> as.lmitt() # brings in control-group means of outcome, predictions

ate_read <- lmitt(read_yr1 ~ 1, STARspec, STARplus,
                  dichotomy= cond_at_entry =="small" ~.,
                  offset = cov_adj(y0hat_read),
                  weights = "ate")
show(ate_read)
vcov(ate_read, type = "HC0", cov_adj_rcorrect = "HC0") |> unname()

ate_read_loc <-
    lmitt(read_yr1 ~ race, STARspec, STARplus,
          dichotomy= cond_at_entry =="small" ~.,
          offset = cov_adj(y0hat_read, newdata = STARplus),
          weights = "ate")
show(ate_read_loc)
}
