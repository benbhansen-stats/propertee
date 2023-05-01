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
##' @param design Optional \code{Design}. If not provided, the function will
##' search through the call stack to find one.
##' @param by optional; vector or list connecting names of cluster/unit of
##' assignment variables in \code{design} to cluster/unit of assignment
##' variables in the covariance adjustment data. If a named vector, names should
##' represent variables in the \code{Design} object and values should represent
##' variables in the data. Only needed if: 1) columns not related to the design should
##' be used for merging the covariance adjustment and quasiexperimental samples, or
##' 2) the column names differ between the datasets used to fit both models.
##' @return \code{SandwichLayer} or \code{PreSandwichLayer} object; the former if
##' `design` is provided or a `design` can be found in the call stack, otherwise
##' the latter. The values represent the covariance adjustments for the
##' observations in `newdata`, if `newdata` is provided or found as an argument to
##' \code{lmitt.formula}, or the fitted values from `model`. The length of the
##' output of \code{cov_adj()} varies depending on this logic.
##' @export
##' @example inst/examples/cov_adj.R
cov_adj <- function(model, newdata = NULL, design =  NULL, by = NULL) {
  if (is.null(design)) {
    design <- .get_design(NULL_on_error = TRUE)
  }
  
  if (is.null(newdata)) {
    form <- .update_ca_model_formula(model, by, design)
    newdata <- tryCatch(
      .get_data_from_model("cov_adj", form),
      error = function(e) {
        warning(paste("Could not find quasiexperimental data in the call stack,",
                      "or it did not contain the columns specified in `by`.",
                      "Using the covariance adjustment data to generate",
                      "the covariance adjustments"), call. = FALSE)
        tryCatch({
          data_call <- model$call$data
          if (is.null(data_call)) {
            stop("`model` must be fit using a `data` argument")
          }
          data <- eval(data_call, envir = environment(formula(model)))
          stats::model.frame(form, data, na.action = na.pass)
        }, error = function(e) {
          stop(paste(e$message, "in covariance adjustment data"), call. = FALSE)
        })
      })
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
    return(as.SandwichLayer(psl, design, by = by, Q_data = newdata))
  }
}

##' (Internal) Add columns for merging covariance adjustment and quasiexperimental samples
##' to model formula
##' @details This function is typically used prior to \code{.get_data_from_model}
##' and incorporates information provided in a `by` vector to ensure the necessary
##' columns for merging the two samples are included in any \code{model.frame} calls.
##' @inheritParams cov_adj
##' @return formula
##' @keywords internal
.update_ca_model_formula <- function(model, by = NULL, design = NULL) {
  form <- deparse(formula(model))
  if (!is.null(design)) {
    form <- paste(c(form, paste(var_names(design, "u"), collapse = " + ")), collapse = " + ")
  }
  if (!is.null(by) && is.null(names(by))) {
    form <- paste(c(form, paste(by, collapse = " + ")), collapse = " + ")
  } else if (!is.null(names(by))) {
    names(by)[names(by) == ""] <- by[names(by) == ""]
    form <- paste(c(form, paste(names(by), collapse = " + ")), collapse = " + ")
  }

  form <- as.formula(form)
  return(form)
}
