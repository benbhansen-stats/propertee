##' @title Produce confidence intervals for linear models
##' @details \code{.confint_lm()} is a copy of \code{stats::confint.lm} but
##'   passes arguments in \code{...} to the \code{vcov()} call. When called on a
##'   \code{teeMod} model, this produces confidence intervals where standard
##'   errors are computed based on the desired formulation of the
##'   \code{vcov_tee()} call.
##'
##' @param object a fitted \code{teeMod} model
##' @param parm a specification of which parameters are to be given confidence
##'   intervals, either a vector of numbers or a vector of names. If missing,
##'   all parameters are considered.
##' @param level the confidence level required.
##' @param ... additional arguments to pass to \code{vcov.teeMod()}
##' @return A matrix (or vector) with columns giving lower and upper confidence
##'   limits for each parameter. These will be labelled as (1-level)/2 and 1 -
##'   (1-level)/2 in % (by default 2.5% and 97.5%)
##' @keywords internal
.confint_lm <- function(object, parm, level = 0.95, ...) {
  cf <- stats::coef(object)
  ses <- sqrt(diag(vcov(object, ...)))
  pnames <- names(ses)
  if (is.matrix(cf))
    cf <- stats::setNames(as.vector(cf), pnames)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- stats::qt(a, object$df.residual)
  pct <- paste(round(a, 3) * 100, "%")
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm,
                                                                   pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}
