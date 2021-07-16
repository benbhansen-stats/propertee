##' @title Covariance Adjustment for Treatment Estimation
##' @param model Any model which supports a `predict` function
##' @return Covariate adjusted outcomes
##' @export
##' @example inst/examples/cov_adj.R
cov_adj <- function(model) {
  # TODO: support predict(..., type = "response"/"link"/other?)
  covAdj <- tryCatch(stats::predict(model, type = "response"),
                     error = function(e) {
                       stop(paste("covariate adjustment model",
                                  "must support predict function"))
                     })
  covAdj
}
