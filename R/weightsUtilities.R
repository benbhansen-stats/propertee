# Internal function to expand uoa-level weights to the level of the data
.join_design_weights <- function(weights, design, target, data = NULL) {

  if (nrow(data) != nrow(design@structure)) {
    # Merge uoa data with weights at uoa level
    uoadata <- design@structure[, design@columnIndex == "u", drop = FALSE]
    uoadata$Design_weights <- weights

    # Merge with data to expand weights to unit of analysis level
    merged <- merge(data, uoadata, by = colnames(uoadata)[-ncol(uoadata)])

    # Extract weights from merged data
    weights <- merged$Design_weights
  }

  WeightedDesign(weights, Design = design, target = target)
}

# Internal function to use by to update the design with new variable names
.update_by <- function(design, data, by) {
  if (!is.list(by) ||
        is.null(names(by)) ||
        any(names(by) == "")) {
    stop("by must be named list")
  }
  if (any(duplicated(names(by))) || any(duplicated(by))) {
    stop("by must be unique")
  }

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
