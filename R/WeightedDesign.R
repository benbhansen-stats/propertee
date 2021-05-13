WeightedDesign <- setClass("WeightedDesign",
                           contains = "numeric",
                           slots = c(Design = "Design",
                                     target = "character"))



##' @title Generate Direct Adjusted Weights
##' @param design a Design object created by one of `RCT_Design`, `RD_Design`,
##'   or `Obs_Design`.
##' @param data optionally the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the `weights` argument of a
##'   model.
##' @return a WeightedDesign object
##' @export
##' @rdname WeightCreators
ett <- function(design, data = NULL) {
  #### generate weights
  weights <- nrow(design@structure):1

  joinDesignWeights(weights, design, target = "ett", data = data)
}

##' @export
##' @rdname WeightCreators
ate <- function(design, data = NULL) {
  #### generate weights
  weights <- 1:nrow(design@structure)

  joinDesignWeights(weights, design, target = "ate", data = data)
}


joinDesignWeights <- function(weights, design, target, data = NULL) {

  if (is.null(data)) {
    data <- get("data", envir = sys.frame(-4))
  }

  if (nrow(data) != nrow(design@structure)) {
    # Merge cluster data with weights at cluster level
    clusterdata <- design@structure[, design@columnIndex == "c", drop = FALSE]
    clusterdata$Design_weights = weights

    # Merge with data to expand weights to unit of analysis level
    merged <- merge(data, clusterdata, by = colnames(clusterdata)[-ncol(clusterdata)])

    # Extract weights from merged data
    weights <- merged$Design_weights
  }

  WeightedDesign(weights, Design = design, target = target)
}
