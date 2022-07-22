#' @include SandwichLayer.R
NULL

##' Covariance Adjustment for Treatment Estimation
##'
##' Prior to obtaining predicted values, \code{cov_adj()} tries to identify the
##' treatment variable (as specified in the \code{design}) and replace it with
##' the reference level. If the treatment is binary, this is \code{FALSE}. If
##' treatment is numeric, it is the smallest non-negative value (note that this
##' means for 0/1 binary, it uses a 0). Factor treatments are not currently
##' supported, but if we add them, it will use the first \code{level()} of the
##' factor, you may change this by using \code{relevel()} to adjust.
##' @param model Any model that inherits from a \code{glm}, \code{lm}, or
##'   \code{robustbase::lmrob} object
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
         warning(paste("Could not find quasiexperimental data in the call",
                       "stack, using the covariance model data to generate",
                       "the covariance adjustments"))
         stats::model.frame(model)
       })
  }

  if (is.null(design)) {
    design <- .get_design(NULL_on_error = TRUE)
  }

  if (!is.null(design)) {
    trt_name <- var_names(design,'t')
    if (trt_name %in% names(newdata))
      if (is.numeric(treatment(design)[, 1])) {
        newdata[[trt_name]] <- min(abs(treatment(design)[, 1]))
      } else if (is.logical(treatment(design)[, 1])) {
        newdata[[trt_name]] <- FALSE
      } else if (is.factor(treatment(design)[, 1])) {
        newdata[[trt_name]] <- levels(treatment(design)[, 1])[1]
      } else {
        warning(paste("The treatment variable is in the covariance adjustment",
                      "model, and is neither logical or numeric; for now,",
                      "partial residuals only implemented for logical or",
                      "numeric treatments"))
      }
  }

  ca_and_grad <- .get_ca_and_prediction_gradient(model, newdata)
  psl <- new("PreSandwichLayer",
             ca_and_grad$ca,
             fitted_covariance_model = model,
             prediction_gradient = ca_and_grad$prediction_gradient)


  if (is.null(design)) {
    return(psl)
  } else {
    return(as.SandwichLayer(psl, design))
  }
}
