# Internal function to expand uoa-level weights to the level of the data
.join_design_weights <- function(weights, design, target, data = NULL) {
  # Merge uoa data with weights at uoa level
  uoadata <- design@structure[, var_names(design, "u"), drop = FALSE]
  uoadata$design_weights <- weights


  # Merge with data to expand weights to unit of analysis level
  weights <- .merge_preserve_order(data, uoadata,
                                  by = var_names(design, "u"))$design_weights

  WeightedDesign(weights,
                 Design = design,
                 target = target)
}

# Internal function to merge data.frames ensuring order of first data.frame is
# maintained
.merge_preserve_order <- function(x, ...) {
  x$..orig_ordering.. <- seq_len(nrow(x))
  x <- merge(x, ...)
  x <- x[order(x$..orig_ordering..), ]
  x$..orig_ordering.. <- NULL
  return(x)
}

# Internal function to use by to update the design with new variable names
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
