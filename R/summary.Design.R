##' @title Summary of Design object
##' @param object Design object
##' @param ... Other args
##' @return object of class "summary.Design"
##' @export
summary.Design <- function(object, ...) {
  out <- list()
  out$Design <- object
  class(out) <- "summary.Design"
  out
}

##' @title Print summary object
##' @param object Design object
##' @param ... Other args
##' @return object, invisibly
##' @export
print.summary.Design <- function(object, ...) {
  des <- object$Design

  destype <- switch(des@type,
                    "RCT" = "Randomized Control Trial",
                    "RD" = "Regression Discontinuity Design",
                    "Obs" = "Observational Study")

  cat("\n")
  cat(destype)
  cat("\n\n")

  cat(paste("Contains ", nrow(des@structure) -
                           sum(duplicated(des@structure[des@columnIndex == "u"])),
            paste0(" ", des@unitOfAssignmentType, "s ") ,"(`",
            paste(varNames(des, "u"), collapse = "`, `"),
            "`)",
            sep = ""))

  if (length(varNames(des, "b")) > 0) {
    cat(paste(" nested within ",
              nrow(des@structure) -
                sum(duplicated(des@structure[des@columnIndex == "b"])),
              " blocks (`",
              paste(varNames(des, "b"), collapse = "`, `"),
              "`)",
              sep = ""))
  }

  # Without a binary treatment variable, what do we want to print here? #11
  #cat("\n")
  #ntreat <- sum(des@structure[des@columnIndex == "t"])
  #nctrl <- sum(1 - des@structure[des@columnIndex == "t"])

  #cat(paste(ntreat, " units  assigned to treatment, ",
  #          nctrl, " units assigned to control (`",
  #          varNames(des, "t"), "`)",
  #          sep = ""))
  cat("\n")

  if (length(varNames(des, "f")) > 0) {
    cat(paste("Forcing variable(s): `",
              paste(varNames(des, "f"), collapse = "`, `"),
              "`",
              sep = ""))
    cat("\n")
  }

  cat("\n")
  invisible(des)
}
