##' @title Covariance Adjustment for Treatment Estimation
##' @param model Any model which supports a `predict` function
##' @param newdata New data
##' @param design Optional `Design`.
##' @return Covariate adjusted outcomes
##' @export
##' @example inst/examples/cov_adj.R
cov_adj <- function(model, newdata = NULL, design =  NULL) {
  if (is.null(design)) {
    # Identify all frames with a weights argument
    weights_args <- lapply(sys.calls(), `[[`, "weights")
    # Loop over each frame which has a weight argument.
    # Its most likely the first frame, but perhaps not.
    for (i in which(!vapply(weights_args, is.null, logical(1)))) {
      possible_design <- get("weights", sys.frame(i))
      if (is(possible_design, "WeightedDesign")) {
        # If we have a WeightedDesign, save it and break
        design <- possible_design
        break()
      }
    }
  }
  # If we weren't able to find the Design, send an error.
  if (is.null(design)) {
    stop("Unable to locate Design in call stack, please use the `design` argument to pass a Design object.")
  }

  # TODO: support predict(..., type = "response"/"link"/other?)
  ca <- tryCatch(stats::predict(model, type = "response", newdata = newdata),
                 error = function(e) {
                   stop(paste("covariate adjustment model",
                              "must support predict function"))
                 })
  ca
}
