##' @title Covariance Adjustment for Treatment Estimation
##' @param model Any model which supports a `predict` function
##' @return Covariate adjusted outcomes
##' @export
cov_adj <- function(model) {
  # TODO: support predict(..., type = "response"/"link"/other?)
  covAdj <- tryCatch(stats::predict(model, type = "response"),
                     error = function(e) {
                       stop("covariate adjustment model must support predict function")
                     })
  covAdj
}
