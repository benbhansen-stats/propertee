##' @title Summary of Design object
##' @param object Design object
##' @param ... Other args
##' @return object of class "summary.Design"
##' @export
summary.Design <- function(object, ...) {
  out <- list()
  out$Design <- object
  out$treatmentTable <- treatmentTable(object)
  class(out) <- "summary.Design"
  out
}

##' @title Print summary object
##' @param x Design object
##' @param ... Other args
##' @param maxUnitPrint Maximum number of treatment levels to print in treatment
##'   table
##' @return object, invisibly
##' @export
print.summary.Design <- function(x, ..., maxUnitPrint = 3) {
  show(x$Design)

  cat("Number of units per Treatment group: \n")
  tt <- x$treatmentTable
  # The value below defines the max number of treatment groups to print
  maxPrintTable <- min(maxUnitPrint, length(tt))
  tdf <- as.data.frame(tt[seq_len(maxPrintTable)])
  colnames(tdf) <- c("Txt Grp", "Num Units")
  # knitr::kable(tt, align = "cc", format = "simple") looks real nice if we want
  # to add that dependency
  if (length(tt) > maxPrintTable) {
    tdf[,1] <- as.character(tdf[,1])
    tdf <- rbind(tdf, c("...", ""))
  }
  print(tdf, row.names = FALSE)
  if (length(tt) > maxPrintTable) {
    if (length(tt) - maxPrintTable == 1) {
      group <- "group"
    } else {
      group <- "groups"
    }
    cat(paste0(length(tt) - maxPrintTable, " smaller treatment ",
               group, " excluded.\n"))
    cat("Use `treatmentTable` function to view full results.")
  }
  cat("\n")
  invisible(x)
}
