##' @param treatment_binary Should the treatment be dichotomized if
##'   \code{object} contains a \code{dichotomy}? Ignored if \code{object} does
##'   not contain a \code{dichotomy}.
##' @export
##' @method summary Design
##' @rdname Design_summary
summary.Design <- function(object, ..., treatment_binary = TRUE) {
  out <- list()
  out$Design <- object
  out$treatment_table <- dtable(object, "treatment", ...)
  class(out) <- "summary.Design"
  return(out)
}

##' @title Summarizing \code{Design} objects
##' @description [summary()] method for class \code{Design}.
##' @param object \code{Design} object, usually a result of a call to
##'   [rct_design()], [obs_design()], or [rd_design()].
##' @param x \code{summary.Design} object, usually as a result of a call to
##'   [summary.Design()]
##' @param ... Ignored
##' @param max_unit_print Maximum number of treatment levels to print in
##'   treatment table
##' @return The \code{Design} or \code{summary.Design}object, invisibly
##' @export
##' @rdname Design_summary
print.summary.Design <- function(x, ..., max_unit_print = 3) {
  show(x$Design)

  cat("Number of units per Treatment group: \n")
  tt <- x$treatment_table
  # The value below defines the max number of treatment groups to print
  max_print_table <- min(max_unit_print, length(tt))
  tdf <- as.data.frame(tt[seq_len(max_print_table)])
  colnames(tdf) <- c("Txt Grp", "Num Units")
  # knitr::kable(tt, align = "cc", format = "simple") looks real nice if we want
  # to add that dependency
  if (length(tt) > max_print_table) {
    tdf[, 1] <- as.character(tdf[, 1])
    tdf <- rbind(tdf, c("...", ""))
  }
  print(tdf, row.names = FALSE)
  if (length(tt) > max_print_table) {
    if (length(tt) - max_print_table == 1) {
      group <- "group"
    } else {
      group <- "groups"
    }
    cat(paste0(length(tt) - max_print_table, " smaller treatment ",
               group, " excluded.\n"))
    cat("Use `dtable` function to view full results.")
  }
  cat("\n")
  invisible(x)
}
