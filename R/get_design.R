# (Internal) Find the inclusion of ate/ett weights and extract the Design
.get_design <- function() {

  design <- NULL

  .find.design <- function(type) {
    stopifnot(type %in% c("weights", "offset"))

    design <- NULL

    # Identify all frames with the appropriate argument
    keyframes <- lapply(sys.calls(), `[[`, type)

    if (length(keyframes) == 0) {
      return(NULL)
    }
    # Loop over each frame which has a weight argument.
    # Its most likely the first frame, but perhaps not.
    for (i in which(!vapply(keyframes, is.null, TRUE))) {
      possible_design_holder <- get(type, sys.frame(i))
      if (is(possible_design_holder, "WeightedDesign") ||
          is(possible_design_holder, "CovAdjPrediction")) {
        # If we have a WeightedDesign, save it and break
        design <- possible_design_holder@Design
        break()
      }
    }
    if (!is(design, "Design")) {
      return(NULL)
    }
    return(design)
  }


  weight_design <- NULL
  covadj_design <- NULL
  if (sys.call(-1)[[1]] != ".weights_calc") {
    weight_design <- .find.design("weights")
  }
  if (sys.call(-1)[[1]] != "cov_adj") {
    covadj_design <- .find.design("offset")
  }

  if (is.null(weight_design) && is.null(covadj_design)) {
    stop(paste("Unable to locate Design in call stack, please use the",
               " `design` argument to pass a Design object."))
  }
  if (is(weight_design, "Design") && is(covadj_design, "Design")) {
    if (!identical(weight_design, covadj_design)) {
      stop("`Design`s foundn in both `cov_adj` and weights but differ")
    }
    return(weight_design)
  }

  if (is.null(weight_design)) {
    return(covadj_design)
  } else {
    return(weight_design)
  }

}
