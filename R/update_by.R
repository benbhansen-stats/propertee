##' Helper function used to update the variable names of a \code{Design} when
##' user passes a \code{by=} argument to align variable names between data sets.
##' @title (Internal) Use \code{by} to update \code{Design} with new variable
##'   names
##' @param design A \code{Design}
##' @param data \code{Data set}
##' @param by named vector or list connecting names of cluster/unit of
##'   assignment variables in \code{design} to cluster/unit of assignment
##'   variables in \code{data}. Names represent variables in the Design; values
##'   represent variables in the data.
##' @return A \code{Design} with updated variable names
##' @keywords internal
.update_by <- function(design, data, by) {
  .check_by(by)

  # Ensure all names and replacements are valid
  missingnames <- !(names(by) %in% colnames(design@structure))
  if (any(missingnames)) {
    warning(paste("by labels not found in Design. unknown elements:",
                  paste(names(by)[missingnames], collapse = ", ")))
  }
  missingdata <-  !(by %in% colnames(data))
  if (any(missingdata)) {
    warning(paste("by replacement values not found in data. unknown elements:",
                  paste(by[missingnames], collapse = ", ")))
  }

  # if we have any names or replacements missing in the design or data,
  # there's a warning, and then don't try to replace that element
  by <- by[!missingnames & !missingdata]

  newnames <- vapply(colnames(design@structure), function(x) {
    pos <- names(by) == x
    if (any(pos)) {
      return(by[[which(pos)]])
    }
    return(x)
  }, "character")

  colnames(design@structure) <- newnames
  validObject(design)
  return(design)
}

# Internal function to throw errors is `by` is inappropriate.
.check_by <- function(by) {
  if (!(is.vector(by) || is.list(by)) ||
        is.null(names(by)) ||
        any(names(by) == "")) {
    stop("'by' must be named vector or named list")
  }
  if (any(duplicated(names(by))) || any(duplicated(by))) {
    stop("all entries in 'by' must be unique")
  }
  NULL
}
