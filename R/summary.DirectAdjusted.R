##' If a \code{DirectAdjusted} object is fit with a \code{SandwichLayer} offset,
##' then the usual \code{stats::summary.lm()} output is enhanced by the use of
##' covariance-adjusted sandwich standard errors, with t-test values
##' recalculated to reflect the new standard errors.
##'
##' @title Create summary of \code{DirectAdjusted} object
##' @param object DirectAdjusted object
##' @param vcov.type A string indicating the desired variance estimator.
##'   Currently accepts "CR0", "MB0", or "HC0"
##' @param ... Additional arguments to \code{vcovDA()}, such as the desired
##'   finite sample heteroskedasticity-robust standard error adjustment.
##' @return object of class \code{summary.DirectAdjusted}
##' @export
##' @method summary DirectAdjusted
summary.DirectAdjusted <- function(object,
                                   vcov.type = c("CR0", "MB0", "HC0"),
                                   ...) {
  out <- summary(as(object, "lm"))
  covmat <- vcovDA(object, type = vcov.type, ...)
  out$coefficients[, 2L] <- sqrt(diag(covmat))
  out$coefficients[, 3L] <- out$coefficients[, 1L] / out$coefficients[, 2L]
  out$coefficients[, 4L] <- 2*stats::pt(abs(out$coefficients[, 3L]),
                                        object$df.residual,
                                        lower.tail = FALSE)
  class(out) <- "summary.DirectAdjusted"
  out$vcov.type <- attr(covmat, "type")
  out$call <- object@lmitt_call
  return(out)
}


##' @title Print summary of \code{DirectAdjusted} object
##' @param x \code{summary.DirectAdjusted} object
##' @param digits the number of significant digits to use when printing.
##' @param signif.stars logical. If ‘TRUE’, ‘significance stars’ are printed for
##'   each coefficient.
##' @param ... Other args
##' @return object, invisibly
##' @importFrom stats pt printCoefmat
##' @export
print.summary.DirectAdjusted <- function(x, digits =
                                              max(3L, getOption("digits") - 3L),
                                         signif.stars =
                                           getOption("show.signif.stars"),
                                         ...) {

  to_report <- x$coefficients
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Treatment Effects:\n")
  stats::printCoefmat(to_report, digits = digits)
  cat(paste0("Std. Error calculated via type \"", x$vcov.type, "\"\n\n"))
  invisible(x)
}
