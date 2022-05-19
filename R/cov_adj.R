##' @title Covariance Adjustment for Treatment Estimation
##' @param model Any model which supports a \code{predict} function
##' @param newdata New data
##' @param design Optional \code{Design}.
##' @return Covariate adjusted outcomes
##' @export
##' @example inst/examples/cov_adj.R
cov_adj <- function(model, newdata = NULL, design =  NULL) {
  if (is.null(design)) {
    design <- .get_design()
  }

  # TODO: support predict(..., type = "response"/"link"/other?)
  ca <- tryCatch(stats::predict(model, type = "response", newdata = newdata),
                 error = function(e) {
                   stop(paste("covariate adjustment model",
                              "must support predict function"))
                 })
  return(new("CovAdjPrediction",
             ca,
             Design = design))
}
