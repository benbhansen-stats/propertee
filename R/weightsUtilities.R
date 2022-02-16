# Internal function to ensure agreement in treatment levels
.treatment_concordance <- function(design, data) {
  treatvar <- colnames(design@structure[which(design@columnIndex == "t")])
  ldes <- levels(design@structure[, treatvar, drop = TRUE])
  if (!treatvar %in% names(data)) {
    stop(paste0("Treatment variable '", treatvar, "' not found in data."))
  }
  datatreatment <- .convert_treatment_to_factor(data[, treatvar, drop = TRUE])
  ldat <- levels(datatreatment)

  if (!identical(ldes, ldat)) {
    if (!all(ldes %in% ldat)) {
      warning("Some levels of treatment in Design not found in data")
    }
    if (!all(ldat %in% ldes)) {
      stop("Some levels of treatment in data not found in Design")
    }
  }
}

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
