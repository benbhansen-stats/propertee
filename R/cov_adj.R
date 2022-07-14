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

  if (is.null(design)) { ### I moved this up here
    design <- .get_design(NULL_on_error = TRUE)
  }

  if(!is.null(design)) ### added this
    ##if(design@type=='RD') ## I think we want to set Z=0 for all design types, right?
  {
    trt_name <- var_names(design,'t')
    if(trt_name %in% names(newdata))
      if(is.numeric(newdata[[trt_name]])){
        newdata[[trt_name]] <- 0 ## are the scenarios where 0 is not the control value?
      } else if(is.logical(newdata[[trt_name]])){
        newdata[[trt_name]] <- FALSE
      } else warning(paste("The treatment variable is in the covariance adjustment model,",
                           "and is neither logical or numeric; for now, partial residuals",
                           "only implemented for logical or numeric treatments"))
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
