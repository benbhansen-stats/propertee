##' @title Summary of Design object
##' @param object Design object
##' @param ... Other args
##' @return object, invisibly
##' @export
summary.Design <- function(object, ...) {
  destype <- switch(object@type,
                    "RCT" = "Randomized Control Trial",
                    "RD" = "Regression Discontinuity Design",
                    "Obs" = "Observational Study")

  cat("\n")
  cat(destype)
  cat("\n\n")

  cat(paste("Contains ", nrow(object@structure) -
                           sum(duplicated(object@structure[object@columnIndex == "u"])),
            paste0(" ", object@unitOfAssignmentType, "s ") ,"(`",
            paste(varNames(object, "u"), collapse = "`, `"),
            "`)",
            sep = ""))

  if (length(varNames(object, "b")) > 0) {
    cat(paste(" nested within ",
              nrow(object@structure) -
                sum(duplicated(object@structure[object@columnIndex == "b"])),
              " blocks (`",
              paste(varNames(object, "b"), collapse = "`, `"),
              "`)",
              sep = ""))
  }

  # Without a binary treatment variable, what do we want to print here? #11
  #cat("\n")
  #ntreat <- sum(object@structure[object@columnIndex == "t"])
  #nctrl <- sum(1 - object@structure[object@columnIndex == "t"])

  #cat(paste(ntreat, " units  assigned to treatment, ",
  #          nctrl, " units assigned to control (`",
  #          varNames(object, "t"), "`)",
  #          sep = ""))
  cat("\n")

  if (length(varNames(object, "f")) > 0) {
    cat(paste("Forcing variable(s): `",
              paste(varNames(object, "f"), collapse = "`, `"),
              "`",
              sep = ""))
    cat("\n")
  }

  cat("\n")
  invisible(object)
}
