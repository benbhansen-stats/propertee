##' @title Create summary of \code{DirectAdjusted} object
##' @param object DirectAdjusted object
##' @param ... Other args
##' @return object of class \code{summary.DirectAdjusted}
##' @export
##' @method summary DirectAdjusted
summary.DirectAdjusted <- function(object, ...) {
  out <- object
  class(out) <- "summary.DirectAdjusted"
  return(out)
}


##' @title Print summary of \code{DirectAdjusted} object
##' @param x \code{summary.DirectAdjusted} object
##' @param ... Other args
##' @param max_unit_print Maximum number of treatment levels to print in
##'   treatment table
##' @return object, invisibly
##' @importFrom stats pt printCoefmat
##' @export
print.summary.DirectAdjusted <- function(x, ..., max_unit_print = 3) {
  # temporary code - produces an annoying warning. Eventually to
  # be replaced with manual calculations
  cf <- x$coefficients[!grepl("^\\.absorbed\\(", names(x$coefficients))]
  se <- sqrt(diag(vcov(x)))
  se <- se[!grepl("^\\.absorbed\\(", names(se))]
  tv <- cf/se
  pv <- 2 * stats::pt(abs(tv), x$df.residual, lower.tail = FALSE)

  to_report <- cbind(cf, se, tv, pv)
  stats::printCoefmat(to_report, digits = 3)
  invisible(x)
}
