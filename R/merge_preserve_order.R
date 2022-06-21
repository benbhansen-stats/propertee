##' @title (Internal) Merge \code{data.frame}s ensuring order of first
##'   \code{data.frame} is maintained
##' @param x \code{data.frame} whose ordering is to be maintained
##' @param ... Additional arguments to \code{merge()}, particularly a second
##'   \code{data.frame} and a \code{by=} argument.
##' @return Merged \code{data.frame} with the same ordering as \code{x}.
##' @keywords internal
.merge_preserve_order <- function(x, ...) {
  x$..orig_ordering.. <- seq_len(nrow(x))
  x <- merge(x, ...)
  x <- x[order(x$..orig_ordering..), ]
  x$..orig_ordering.. <- NULL
  return(x)
}
