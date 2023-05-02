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

  if (object$rank > 0) {
    covmat <- vcovDA(object, type = vcov.type, ...)
    out$coefficients[, 2L] <- sqrt(diag(covmat))
    out$coefficients[, 3L] <- out$coefficients[, 1L] / out$coefficients[, 2L]
    out$coefficients[, 4L] <- 2*stats::pt(abs(out$coefficients[, 3L]),
                                          object$df.residual,
                                          lower.tail = FALSE)
    out$vcov.type <- attr(covmat, "type")
  }
  class(out) <- "summary.DirectAdjusted"
  out$lmitt <- object
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
print.summary.DirectAdjusted <- function(x,
                                         digits =
                                           max(3L, getOption("digits") - 3L),
                                         signif.stars =
                                           getOption("show.signif.stars"),
                                         ...) {

  cat("\nCall:\n", paste(deparse(x$lmitt@lmitt_call), sep = "\n",
                         collapse = "\n"), "\n", sep = "")

  df <- x$df

  if (x$lmitt@lmitt_fitted) {
    coefmatname <- "Treatment Effects"
  } else {
    coefmatname <- "Coefficients"
  }


  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    if (nsingular <- df[3L] - df[1L])
      cat("\n", coefmatname, ": (", nsingular,
          " not defined because of singularities)\n",
          sep = "")
    else cat("\n", coefmatname, ":\n")
    coefs <- x$coefficients
    if (any(aliased <- x$aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4,
                      dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    stats::printCoefmat(coefs,
                        digits = digits,
                        signif.stars = signif.stars,
                        na.print = "NA", ...)
  }

  if (nrow(x$coefficients) > 0) {
    # Only print if we estimate at least one coefficient
    cat(paste0("Std. Error calculated via type \"", x$vcov.type, "\"\n\n"))
  } else {
    cat("\n")
  }
  invisible(x)
}
