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
    # TODO: replace this with the new .get_design() call
    design <- tryCatch(
      error = function(e) {
        message(paste("Unable to locate Design in call stack, use the `design`",
                      "argument in `lmitt` or `as.DirectAdjusted` to pass a",
                      "Design object."))
      },
      message = function(m) {
        NULL
      },
      .get_design())
  }

  # if (is.null(design)) {
  #   return(psl)
  # } else {
  #   return(as.SandwichLayer(psl, design))
  # }
  return(psl)
}
