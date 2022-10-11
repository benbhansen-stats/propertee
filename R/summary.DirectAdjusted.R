##' If a \code{DirectAdjusted} object is fit with a \code{SandwichLayer} offset,
##' then the usual \code{stats::summary.lm()} output is enhanced by the use of
##' covariance-adjusted sandwich standard errors, with t-test values
##' recalculated to reflect the new standard errors.
##'
##' @title Create summary of \code{DirectAdjusted} object
##' @param object DirectAdjusted object
##' @param ... Additional arguments to \code{vcovDA()}, such as the desired
##' finite sample heteroskedasticity-robust standard error adjustment.
##' @return object of class \code{summary.DirectAdjusted}
##' @export
##' @method summary DirectAdjusted
summary.DirectAdjusted <- function(object, ...) {
  out <- summary(as(object, "lm"))
  if (inherits(object$model$`(offset)`, "SandwichLayer")) {
    out$coefficients[, 2L] <- sqrt(diag(vcovDA(object, ...)))
    out$coefficients[, 3L] <- out$coefficients[, 1L] / out$coefficients[, 2L]
    out$coefficients[, 4L] <- 2*stats::pt(abs(out$coefficients[, 3L]),
                                          object$df.residual,
                                          lower.tail = FALSE)
  }
  class(out) <- "summary.DirectAdjusted"
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
print.summary.DirectAdjusted <-
  function(x,
           digits = max(3L, getOption("digits") - 3L),
           signif.stars = getOption("show.signif.stars"),
           ...) {

  to_report <- x$coefficients
  to_report <- to_report[!grepl("^\\.absorbed\\(", rownames(to_report)), ]
  stats::printCoefmat(to_report, digits = digits)
  invisible(x)
}
