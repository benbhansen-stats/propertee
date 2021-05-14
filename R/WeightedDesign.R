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
##' @param clusterIds optional; list connecting names of cluster variables in
##'   `design` to cluster variables in `data`
##' @return a WeightedDesign object
##' @export
##' @rdname WeightCreators
ett <- function(design, data = NULL, clusterIds = NULL) {
  if (is.null(data)) {
    data <- get_data_from_model()
  }

  if (!is.null(clusterIds)) {
    design <- update_clusterIds(design, data, clusterIds)
  }


  #### generate weights
  weights <- rev(seq_len(nrow(design@structure)))

  joinDesignWeights(weights, design, target = "ett", data = data)
}

##' @export
##' @rdname WeightCreators
ate <- function(design, data = NULL, clusterIds = NULL) {
  if (is.null(data)) {
    data <- get_data_from_model()
  }

  if (!is.null(clusterIds)) {
    design <- update_clusterIds(design, data, clusterIds)
  }

  #### generate weights
  weights <- seq_len(nrow(design@structure))

  joinDesignWeights(weights, design, target = "ate", data = data)
}

# Internal function to use clusterIds to update the design with new variable
# names
update_clusterIds <- function(design, data, clusterIds) {
  if (!is.list(clusterIds) ||
        is.null(names(clusterIds)) ||
        any(names(clusterIds) == "")) {
    stop("clusterIds must be named list")
  }
  if (any(duplicated(names(clusterIds))) || any(duplicated(clusterIds))) {
    stop("clusterIds must be unique")
  }

  # Ensure all names and replacements are valid
  missingnames <- !(names(clusterIds) %in% colnames(design@structure))
  if (any(missingnames)) {
    warning(paste("clusterIds labels not found in Design. unknown elements:",
                  paste(names(clusterIds)[missingnames], collapse = ", ")))
  }
  missingdata <-  !(clusterIds %in% colnames(data))
  if (any(missingdata)) {
    warning(paste("clusterIds replacement values not found in data. unknown elements:",
                  paste(clusterIds[missingnames], collapse = ", ")))
  }

  # if we have any names or replacements missing in the design or data,
  # there's a warning, and then don't try to replace that element
  clusterIds <- clusterIds[!missingnames & !missingdata]


  newnames <- vapply(colnames(design@structure), function(x) {
    pos <- names(clusterIds) == x
    if (any(pos)) {
      return(clusterIds[[which(pos)]])
    }
    return(x)
  }, "character")

  colnames(design@structure) <- newnames
  return(design)
}

# Internal functoin to try and retrieve the data from the model when `ate` or
# `ett` are called without a data argument
get_data_from_model <- function() {
  data <- get("data", envir = sys.frame(-4))
  return(data)
}

# Internal function to expand cluster-level weights to the level of the data
joinDesignWeights <- function(weights, design, target, data = NULL) {

  if (nrow(data) != nrow(design@structure)) {
    # Merge cluster data with weights at cluster level
    clusterdata <- design@structure[, design@columnIndex == "c", drop = FALSE]
    clusterdata$Design_weights <- weights

    # Merge with data to expand weights to unit of analysis level
    merged <- merge(data, clusterdata, by = colnames(clusterdata)[-ncol(clusterdata)])

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
