WeightedDesign <- setClass("WeightedDesign",
                           contains = "numeric",
                           slots = c(Design = "Design",
                                     target = "character"))

setValidity("WeightedDesign", function(object) {
  if (any(object < 0)) {
    return("Weights must be non-negative")
  }
  if (all(object == 0)) {
    return("At least one weight must be positive")
  }
  if (!object@target %in% c("ate", "ett")) {
    return(paste("@target must be one of [ate,ett]. unknown @target:", object@target))
  }
  TRUE
})



##' @title Generate Direct Adjusted Weights
##' @param design a Design object created by one of `RCT_Design`, `RD_Design`,
##'   or `Obs_Design`.
##' @param data optionally the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the `weights` argument of a
##'   model.
##' @param unitOfAssignmentIds optional; list connecting names of unit of
##'   assignment variables in `design` to unit of assignment variables in `data`
##' @return a WeightedDesign object
##' @export
##' @rdname WeightCreators
ett <- function(design, data = NULL, unitOfAssignmentIds = NULL) {
  if (is.null(data)) {
    data <- .get_data_from_model(design@call$formula, unitOfAssignmentIds)
  }

  if (!is.null(unitOfAssignmentIds)) {
    design <- .update_unitOfAssignmentIds(design, data, unitOfAssignmentIds)
  }

  .treatment_concordance(design, data)

  #### generate weights
  weights <- rev(seq_len(nrow(design@structure)))

  .join_design_weights(weights, design, target = "ett", data = data)
}

##' @export
##' @rdname WeightCreators
ate <- function(design, data = NULL, unitOfAssignmentIds = NULL) {
  if (is.null(data)) {
    data <- .get_data_from_model(design@call$formula, unitOfAssignmentIds)
  }

  if (!is.null(unitOfAssignmentIds)) {
    design <- .update_unitOfAssignmentIds(design, data, unitOfAssignmentIds)
  }

  .treatment_concordance(design, data)

  #### generate weights
  weights <- seq_len(nrow(design@structure))

  .join_design_weights(weights, design, target = "ate", data = data)
}

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

# Internal function to use unitOfAssignmentIds to update the design with new variable
# names
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


##' @title Show a WeightedDesign
##' @param object WeightedDesignDesign object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "WeightedDesign", function(object) {
  print(object@.Data)
  invisible(object)
})


##' WeightedDesigns do not support addition or subtraction, but do support all
##' other reasonable operations.
##'
##' @title WeightedDesign Ops
##' @param e1 WeightedDesign or numeric
##' @param e2 numeric or WeightedDesign
##' @rdname WeightedDesignOps
##' @export
setMethod("+", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) addsubtracterror() )

##' @rdname WeightedDesignOps
##' @export
setMethod("+", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) addsubtracterror() )

##' @rdname WeightedDesignOps
##' @export
setMethod("-", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) addsubtracterror() )

##' @rdname WeightedDesignOps
##' @export
setMethod("-", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) addsubtracterror() )

##' @rdname WeightedDesignOps
##' @export
setMethod("*", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) {
            e1@.Data <- e1@.Data*e2
            validObject(e1)
            e1
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("*", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) {
            e2@.Data <- e1*e2@.Data
            validObject(e2)
            e2
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("/", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) {
            e1@.Data <- e1@.Data/e2
            validObject(e1)
            e1
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("/", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) {
            e2@.Data <- e1/e2@.Data
            validObject(e2)
            e2
          })

addsubtracterror <- function() {
  stop("Cannot perform addition or subtraction on WeightedDesigns")
}


setGeneric("weights")

##' @title Extract Weights from WeightedDesign
##' @param object WeightedDesign object
##' @param ... Ignored
##' @return Weights
##' @export
setMethod("weights", "WeightedDesign", function(object, ...) {
  as.numeric(object)
})
