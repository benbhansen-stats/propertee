#' @include Design.R DesignAccessors.R
NULL

setClass("PreSandwichLayer",
         contains = "numeric",
         slots = c(fitted_covariance_model = "ANY",
                   prediction_gradient = "matrix"))

setValidity("PreSandwichLayer", function(object) {
  if (!("terms" %in% attr(object@fitted_covariance_model, "names"))) {
    return("Fitted covariance model must have a 'terms' attribute")
  }
  
  tryCatch({
    sandwich::bread(object@fitted_covariance_model)
    sandwich::estfun(object@fitted_covariance_model)
    }, error = function(e) {
      stop(paste("Functions for extracting vcov elements not applicable to",
                 "fitted covariance model"))
  })
  
  if (!is.numeric(object@prediction_gradient)) {
    return("Prediction gradient must be a numeric matrix")
  }
  if (dim(object@prediction_gradient)[1] != length(object)) {
    msg <- paste0("Prediction gradient matrix (",
                  paste(dim(object@prediction_gradient), collapse = ", "),
                  ") does not have the same dimension along axis 1 as the ",
                  "covariance adjustment vector (",
                  length(object))
    return(msg)
  }
  if (dim(object@prediction_gradient)[2] != length(object@fitted_covariance_model$coefficients)) {
    return(paste0("Prediction gradient does not have the same number of columns as ",
                  "predictors in the covariance model"))
  }
  TRUE
})

setClass("SandwichLayer",
         contains = "numeric",
         slots = c(fitted_covariance_model = "ANY",
                   prediction_gradient = "matrix",
                   keys = "data.frame",
                   Design = "Design"))

setValidity("SandwichLayer", function(object) {
  psl <- new("PreSandwichLayer",
             object@.Data,
             fitted_covariance_model = object@fitted_covariance_model,
             prediction_gradient = object@prediction_gradient)
  validObject(psl)
  
  if (nrow(object@keys) != nrow(model.matrix(object@fitted_covariance_model))) {
    return(paste0("Keys does not have the same number of rows as the dataset used ",
                  "to fit the covariance model"))
  }

  if (any(is.na(object))) {
    msg <- paste("Some covariance adjustments are NA; be careful of dropping",
                 "these observations when fitting the design model")
    warning(msg)
  }
  TRUE
})

show_layer <- function(object) {
  print(object@.Data)
  invisible(object)
}

##' @title Show a PreSandwichLayer
##' @param object PreSandwichLayer object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "PreSandwichLayer", show_layer)

##' @title Show a SandwichLayer
##' @param object SandwichLayer object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "SandwichLayer", show_layer)

##' (Internal) Get a vector of "response" predictions from a covariance model
##' for a certain dataframe and its gradient with respect to the parameters of
##' the covariance model
##' @param model Any model of class \code{glm} or \code{lm} (excluding those from
##' the \CRANpkg{gam} package) that supports \code{predict} and \code{model.matrix}
##' methods
##' @param newdata Optional; a data.frame of new data
##' @return Covariate adjusted outcomes and their gradient with respect to the
##' parameters of the covariance model (a list of a numeric vector and a matrix)
.get_ca_and_prediction_gradient <- function(model, newdata = NULL) {
  if (!is.null(newdata) & !is.data.frame(newdata)) {
    stop("If supplied, `newdata` must be a dataframe")
  }
  
  X <- tryCatch(
    if (is.null(newdata)) {
      stats::model.matrix(model)
    } else {
      form <- as.formula(model$call$formula[-2])
      stats::model.matrix(form,
                          stats::model.frame(form, data = newdata, na.action = 'na.pass'))
    }, error = function(e) {
      stop("`model` must have a `call` object and `model.matrix` method")
    })
  
  # TODO: support predict(..., type = "response"/"link"/other?)
  ca <- tryCatch(stats::predict(model, type = "response", newdata = newdata),
                 error = function(e) {
                   stop(paste("covariate adjustment model",
                              "must support predict function"))
                 })
  
  # this branch applies to `glm`, `surveyglm`, `robustbase::glmrob` models
  if (is(model, "glm") & !is(model, "gam")) {
    pred_gradient <- model$family$mu.eta(ca) * X
  } else if (is(model, "lm") | is(model, "lmrob")) {
    # `lm` doesn't have a `family` object, but we know its prediction gradient
    pred_gradient <- X
  } else {
    stop("`model` must be a `glm` or `lm` object (and not a `gam` object)")
  }
  
  return(list("ca" = ca,
              "prediction_gradient" = pred_gradient))
}


##' @title Convert a PreSandwichLayer to a SandwichLayer via a Design Object
##' @param x a \code{PreSandwichLayer} object.
##' @param design a \code{Design} object created by one of \code{rct_design()},
##' \code{rd_design()}, or \code{obs_design()}.
##' @param by optional; named vector or list connecting names of cluster/unit of
##' assignment variables in \code{design} to cluster/unit of assignment
##' variables in \code{data}. Names represent variables in the Design; values
##' represent variables in the data. Only needed if variable names differ.
##' @return a \code{SandwichLayer} object
##' @export
as.SandwichLayer <- function(x, design, by = NULL) {
  if (!is(x, "PreSandwichLayer")) {
    stop("x must be a `PreSandwichLayer` object")
  }
  
  data_call <- x@fitted_covariance_model$call$data
  if (is.null(data_call)) {
    stop("The fitted covariance model for x must be fit using a `data` argument")
  }

  covmoddata <- eval(data_call,
                     envir = environment(formula(x@fitted_covariance_model)))

  if (!is.null(by)) {
    # .update_by handles checking input
    design <- .update_by(design, covmoddata, by)
  }

  desvars <- var_names(design, "u")
  wide_frame <- tryCatch(
    stats::expand.model.frame(x@fitted_covariance_model, desvars, na.expand = TRUE)[desvars],
    error = function(e) {
      stop(paste("The",
                 gsub("_", " ", design@unit_of_assignment_type),
                 "columns",
                 paste(setdiff(desvars, colnames(covmoddata)), collapse = ", "),
                 "are missing from the covariance model dataset"),
           call. = FALSE)
    })
  keys <- .merge_preserve_order(wide_frame, design@structure, all.x = TRUE, sort = FALSE)
  keys[is.na(keys[, var_names(design, "t")]), desvars] <- NA
  keys <- keys[, desvars, drop = FALSE]
  
  return(new("SandwichLayer",
             x@.Data,
             fitted_covariance_model = x@fitted_covariance_model,
             prediction_gradient = x@prediction_gradient,
             keys = keys,
             Design = design))
}
