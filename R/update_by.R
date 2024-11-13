##' Helper function used to update the variable names of a
##' \code{StudySpecification} when user passes a \code{by=} argument to align
##' variable names between data sets.
##' @title (Internal) Use \code{by} to update \code{StudySpecification} with new
##'   variable names
##' @param specification A \code{StudySpecification}
##' @param data \code{Data set}
##' @param by named vector or list connecting names of unit of
##'   assignment/unitid/cluster variables in \code{specification} to unit of
##'   assignment/unitid/cluster variables in \code{data}. Names represent
##'   variables in the StudySpecification; values represent variables in the
##'   data.
##' @return A \code{StudySpecification} with updated variable names
##' @keywords internal
.update_by <- function(specification, data, by) {
  .check_by(by)

  # Ensure all names and replacements are valid
  missingnames <- !(names(by) %in% colnames(specification@structure))
  if (any(missingnames)) {
    warning(paste("by labels not found in StudySpecification. unknown elements:",
                  paste(names(by)[missingnames], collapse = ", ")))
  }
  missingdata <-  !(by %in% colnames(data))
  if (any(missingdata)) {
    warning(paste("by replacement values not found in data. unknown elements:",
                  paste(by[missingnames], collapse = ", ")))
  }

  # if we have any names or replacements missing in the specification or data,
  # there's a warning, and then don't try to replace that element
  by <- by[!missingnames & !missingdata]

  newnames <- vapply(colnames(specification@structure), function(x) {
    pos <- names(by) == x
    if (any(pos)) {
      return(by[[which(pos)]])
    }
    return(x)
  }, "character")

  colnames(specification@structure) <- newnames
  validObject(specification)
  return(specification)
}

##' Thie ensures that the \code{by=} argument is of the proper type, is named,
##' and consists of only unique entries.
##'
##' @title (Internal) A few checks to ensure \code{by=} is valid
##' @param by named vector or list connecting names of unit of
##'   assignment/unitid/cluster variables in \code{specification} to unit of
##'   assignment/unitid/cluster variables in \code{data}. Names represent
##'   variables in the StudySpecification; values represent variables in the
##'   data.
##' @return \code{NULL} if no errors are found
##' @keywords internal
.check_by <- function(by) {
  if (!(is.vector(by) || is.list(by)) ||
        is.null(names(by)) ||
        any(names(by) == "")) {
    stop("'by' must be named vector or named list")
  }
  if (any(duplicated(names(by))) || any(duplicated(by))) {
    stop("all entries in 'by' must be unique")
  }
  invisible(NULL)
}
