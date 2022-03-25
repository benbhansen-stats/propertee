##' @title Summary of Design object
##' @param object Design object
##' @param ... Other args
##' @return object of class "summary.Design"
##' @export
summary.Design <- function(object, ...) {
  out <- list()
  out$Design <- object
  out$treatment_table <- treatment_table(object)
  class(out) <- "summary.Design"
  out
}

##' @title Print summary object
##' @param x Design object
##' @param ... Other args
##' @param max_unit_print Maximum number of treatment levels to print in treatment
##'   table
##' @return object, invisibly
##' @export
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
    cat("Use `treatment_table` function to view full results.")
  }
  cat("\n")
  invisible(x)
}
