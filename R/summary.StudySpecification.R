##' @param treatment_binary Should the treatment be dichotomized if
##'   \code{object} contains a \code{dichotomy}? Ignored if \code{object} does
##'   not contain a \code{dichotomy}.
##' @export
##' @method summary StudySpecification
##' @rdname StudySpecification_summary
summary.StudySpecification <- function(object, ..., treatment_binary = TRUE) {
  out <- list()
  out$StudySpecification <- object
  out$treatment_table <- stable(object, "treatment", ...)
  class(out) <- "summary.StudySpecification"
  return(out)
}

##' @title Summarizing \code{StudySpecification} objects
##' @description [summary()] method for class \code{StudySpecification}.
##' @param object \code{StudySpecification} object, usually a result of a call
##'   to [rct_spec()], [obs_spec()], or [rd_spec()].
##' @param x \code{summary.StudySpecification} object, usually as a result of a
##'   call to [summary.StudySpecification()]
##' @param ... Ignored
##' @param max_unit_print Maximum number of treatment levels to print in
##'   treatment table
##' @return The \code{StudySpecification} or
##'   \code{summary.StudySpecification}object, invisibly
##' @export
##' @rdname StudySpecification_summary
print.summary.StudySpecification <- function(x, ..., max_unit_print = 3) {
  show(x$StudySpecification)

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
    cat("Use `specification_table` function to view full results.")
  }
  cat("\n")
  invisible(x)
}
