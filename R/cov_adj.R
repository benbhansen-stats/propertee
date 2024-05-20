#' @include SandwichLayer.R
NULL

##' @title Covariance adjustment of \code{teeMod} model estimates
##' @description
##'  \code{cov_adj()} takes a fitted covariance model and returns the information
##'  necessary for adjusting direct adjustment model estimates and associated
##'  standard errors for covariates. Standard errors will reflect adjustments made
##'  to the outcomes as well as contributions to sampling variability arising from
##'  the estimates of the covariance adjustment model coefficients.
##' @details Prior to generating adjustments, \code{cov_adj()} identifies the
##'   treatment variable specified in the \code{Design} object passed to
##'   \code{design} and replaces all values with a reference level. If the treatment has
##'   logical type, this reference level is \code{FALSE}, and if it has numeric type,
##'   this is the smallest non-negative value (which means 0 for 0/1 binary).
##'   Factor treatments are not currently supported for \code{Design} objects.\cr\cr
##'   The values of the output vector represent adjustments for the outcomes in
##'   \code{newdata} if \code{newdata} is provided; adjustments for the outcomes
##'   in the data used to fit a \code{teeMod} model if \code{cov_adj()}
##'   is called within the \code{offset} argument of the model fit; or they are
##'   the fitted values from \code{model} if no relevant dataframe can be extracted
##'   from the call stack. The length of the output of \code{cov_adj()} will match
##'   the number of rows of the dataframe used.
##' @param model any model that inherits from a \code{glm}, \code{lm}, or
##'   \code{robustbase::lmrob} object
##' @param newdata a dataframe of new data. Default is NULL, in which case
##'   a dataframe is sought from higher up the call stack.
##' @param design a \code{Design} object. Default is NULL, in which case a
##'   \code{Design} object is sought from higher up the call stack.
##' @inheritParams as.SandwichLayer
##' @return A \code{SandwichLayer} if \code{design} is not NULL or a \code{Design}
##' object is found in the call stack, otherwise a \code{PreSandwichLayer} object
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
        warning(paste("Could not find direct adjustment data in the call stack,",
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

##' @title (Internal) Add columns for merging covariance adjustment and direct
##' adjustment samples to model formula
##' @details This function is typically used prior to \code{.get_data_from_model()}
##' and incorporates information provided in a \code{by} argument to ensure the necessary
##' columns for merging the two samples are included in any \code{model.frame()} calls.
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
