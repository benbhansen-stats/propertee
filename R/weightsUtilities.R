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

# Internal function to use unitOfAssignmentIds to update the design with new
# variable names
.update_unitOfAssignmentIds <- function(design, data, unitOfAssignmentIds) {
  if (!is.list(unitOfAssignmentIds) ||
        is.null(names(unitOfAssignmentIds)) ||
        any(names(unitOfAssignmentIds) == "")) {
    stop("unitOfAssignmentIds must be named list")
  }
  if (any(duplicated(names(unitOfAssignmentIds))) || any(duplicated(unitOfAssignmentIds))) {
    stop("unitOfAssignmentIds must be unique")
  }

  # Ensure all names and replacements are valid
  missingnames <- !(names(unitOfAssignmentIds) %in% colnames(design@structure))
  if (any(missingnames)) {
    warning(paste("unitOfAssignmentIds labels not found in Design. unknown elements:",
                  paste(names(unitOfAssignmentIds)[missingnames], collapse = ", ")))
  }
  missingdata <-  !(unitOfAssignmentIds %in% colnames(data))
  if (any(missingdata)) {
    warning(paste("unitOfAssignmentIds replacement values not found in data. unknown elements:",
                  paste(unitOfAssignmentIds[missingnames], collapse = ", ")))
  }

  # if we have any names or replacements missing in the design or data,
  # there's a warning, and then don't try to replace that element
  unitOfAssignmentIds <- unitOfAssignmentIds[!missingnames & !missingdata]

  newnames <- vapply(colnames(design@structure), function(x) {
    pos <- names(unitOfAssignmentIds) == x
    if (any(pos)) {
      return(unitOfAssignmentIds[[which(pos)]])
    }
    return(x)
  }, "character")

  colnames(design@structure) <- newnames
  return(design)
}
