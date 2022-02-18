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

# Internal function to use varLinks to update the design with new variable names
.update_varLinks <- function(design, data, varLinks) {
  if (!is.list(varLinks) ||
        is.null(names(varLinks)) ||
        any(names(varLinks) == "")) {
    stop("varLinks must be named list")
  }
  if (any(duplicated(names(varLinks))) || any(duplicated(varLinks))) {
    stop("varLinks must be unique")
  }

  # Ensure all names and replacements are valid
  missingnames <- !(names(varLinks) %in% colnames(design@structure))
  if (any(missingnames)) {
    warning(paste("varLinks labels not found in Design. unknown elements:",
                  paste(names(varLinks)[missingnames], collapse = ", ")))
  }
  missingdata <-  !(varLinks %in% colnames(data))
  if (any(missingdata)) {
    warning(paste("varLinks replacement values not found in data. unknown elements:",
                  paste(varLinks[missingnames], collapse = ", ")))
  }

  # if we have any names or replacements missing in the design or data,
  # there's a warning, and then don't try to replace that element
  varLinks <- varLinks[!missingnames & !missingdata]

  newnames <- vapply(colnames(design@structure), function(x) {
    pos <- names(varLinks) == x
    if (any(pos)) {
      return(varLinks[[which(pos)]])
    }
    return(x)
  }, "character")

  colnames(design@structure) <- newnames
  return(design)
}
