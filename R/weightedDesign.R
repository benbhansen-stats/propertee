WeightedDesign <- setClass("WeightedDesign",
                           contains = "numeric",
                           slots = c(Design = "Design",
                                     estimand = "character"))



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
  weights <- rep(1, nrow(design@structure))

  joinDesignWeights(weights, design, estimand = "ett", data = data)
}

##' @export
##' @rdname WeightCreators
ate <- function(design, data = NULL) {
  #### generate weights
  weights <- rep(1, nrow(design@structure))

  joinDesignWeights(weights, design, estimand = "ate", data = data)
}


joinDesignWeights <- function(weights, design, estimand, data = NULL) {

  # If data is NULL, extract from environment

  # Join weights & design@structure, then join with data.

  # extract weights separately.

  WeightedDesign(weights, Design = design, estimand = estimand)
}
