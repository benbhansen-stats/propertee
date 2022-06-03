#' @include SandwichLayer.R
NULL

##' @title Covariance Adjustment for Treatment Estimation
##' @param model Any model of class \code{glm} or \code{lm} (excluding those from
##' the \CRANpkg{gam} package) that supports \code{predict} and \code{model.matrix}
##' methods
##' @param newdata Optional; a data.frame of new data
##' @param design Optional \code{Design}.
##' @return Covariate adjusted outcomes
##' @export
##' @example inst/examples/cov_adj.R
cov_adj <- function(model, newdata = NULL, design =  NULL) {
  if (is.null(design)) {
    design <- .get_design()
  }

  if (is.null(newdata)) {
    form <- model$call$formula
    newdata <- tryCatch(
      .get_data_from_model("weights", form),
       error = function(e) {
         stop(paste("cov_adj must be called with a newdata argument if not called",
                    "as an offset argument"))
       })
  }
  ca_and_grad <- .get_ca_and_prediction_gradient(model, newdata)
  psl <- new("PreSandwichLayer",
             ca_and_grad$ca,
             fitted_covariance_model = model,
             prediction_gradient = ca_and_grad$prediction_gradient)
  
  return(psl)
}
