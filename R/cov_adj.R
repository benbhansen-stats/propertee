#' @include SandwichLayer.R
NULL

##' @title Covariance Adjustment for Treatment Estimation
##' @param model Any model that inherits from a \code{glm}, \code{lm}, or \code{
##' robustbase::lmrob} object
##' @param newdata Optional; a data.frame of new data
##' @param design Optional \code{Design}.
##' @return Covariate adjusted outcomes
##' @export
##' @example inst/examples/cov_adj.R
cov_adj <- function(model, newdata = NULL, design =  NULL) {
  if (is.null(newdata)) {
    form <- model$call$formula
    newdata <- tryCatch(
      .get_data_from_model("weights", form),
       error = function(e) {
         warning(paste("Could not find quasiexperimental data in the call stack,",
                       "using the covariance model data to generate the covariance",
                       "adjustments"))
         stats::model.frame(model)
       })
  }


  ca_and_grad <- .get_ca_and_prediction_gradient(model, newdata)
  psl <- new("PreSandwichLayer",
             ca_and_grad$ca,
             fitted_covariance_model = model,
             prediction_gradient = ca_and_grad$prediction_gradient)

  if (is.null(design)) {
    design <- .get_design(NULL_on_error = TRUE)
  }

  if (is.null(design)) {
    return(psl)
  } else {
    return(as.SandwichLayer(psl, design))
  }
}
