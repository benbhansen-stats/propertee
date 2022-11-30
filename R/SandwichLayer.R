#' @include Design.R DesignAccessors.R
NULL

setClass("PreSandwichLayer",
         contains = "numeric",
         slots = c(fitted_covariance_model = "ANY",
                   prediction_gradient = "matrix"))

setValidity("PreSandwichLayer", function(object) {
  if (!inherits(object@fitted_covariance_model$terms, "terms")) {
    return("Fitted covariance model must have a valid 'terms' attribute")
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
                  "covariance adjustment vector (", length(object), ")")
    return(msg)
  }
  if (dim(object@prediction_gradient)[2] !=
      length(object@fitted_covariance_model$coefficients)) {
    return(paste("Prediction gradient does not have the same number of columns",
                 "predictors in the covariance model"))
  }
  TRUE
})

setClass("SandwichLayer",
         contains = "PreSandwichLayer",
         slots = c(keys = "data.frame",
                   Design = "Design"))

setValidity("SandwichLayer", function(object) {
  if (nrow(object@keys) != nrow(stats::model.matrix(object@fitted_covariance_model))) {
    return(paste("Keys does not have the same number of rows as the dataset",
                 "used to fit the covariance model"))
  }

  if (any(is.na(object))) {
    msg <- paste("Some covariance adjustments are NA; be careful of dropping",
                 "these observations when fitting the design model")
    warning(msg)
  }
  TRUE
})

##' @title (Internal) \code{show} helper for
##'   \code{PreSandwichLayer}/\code{SandwichLayer}
##' @param object \code{PreSandwichLayer} or \code{SandwichLayer}
##' @return \code{object}, invisibly
##' @keywords internal
.show_layer <- function(object) {
  print(object@.Data)
  invisible(object)
}

##' @title Show a PreSandwichLayer or SandwichLayer
##' @param object PreSandwichLayer object
##' @return an invisible copy of `object`
##' @export
##' @rdname PreSandwichLayer.show
setMethod("show", "PreSandwichLayer", .show_layer)

setGeneric("subset")

##' @title PreSandwichLayer and SandwichLayer subsetting
##' @param subset Logical vector identifying values to keep or drop
##' @return \code{x} subset by \code{i}
##' @export
##' @rdname PreSandwichLayer.subset
setMethod("subset", "PreSandwichLayer", function(x, subset) {
  x@.Data <- subset(x@.Data, subset = subset)
  x@prediction_gradient <- subset(x@prediction_gradient, subset = subset)
  return(x)
})

setGeneric("[")

##' @param x \code{PreSandwichLayer} object
##' @param i indices specifying elements to extract or replace. See
##'   \code{help("[")} for further details.
##' @return \code{x} subset by \code{i}
##' @export
##' @importFrom methods callNextMethod
##' @rdname PreSandwichLayer.subset
setMethod("[", "PreSandwichLayer",
          function(x, i) {
            idx <- rep(FALSE, length(x@.Data))
            idx[i] <- TRUE
            dat <- methods::callNextMethod()
            x@.Data <- dat
            x@prediction_gradient <- subset(x@prediction_gradient, idx)
            return(x)

          })

##' (Internal) Get the a vector of "response" predictions from a covariance model
##' and its gradient with respect to the fitted coefficients
##' @param model Any model that inherits from a \code{glm}, \code{lm}, or \code{
##' robustbase::lmrob} object
##' @param newdata Optional; a data.frame of new data
##' @return Covariate adjusted outcomes and their gradient with respect to the
##' parameters of the covariance model (a list of a numeric vector and a matrix)
.get_ca_and_prediction_gradient <- function(model, newdata = NULL) {
  if (is.null(newdata)) {
    newdata <- stats::model.frame(model)
  }

  if (!is.data.frame(newdata)) {
    stop("If supplied, `newdata` must be a dataframe")
  }

  model_terms <- tryCatch(
    terms(model),
    error = function(e) stop("`model` must have `terms` method", call. = FALSE)
  )
  newdata <- tryCatch({
    try_terms <- terms(newdata) # terms will be the same if newdata generated with model's model.frame
    if (!is.null(try_terms)) newdata
    }, error = function(e) {
      stats::model.frame(stats::delete.response(model_terms),
                         data = newdata,
                         na.action = na.pass)
  })
  X <- stats::model.matrix(stats::delete.response(model_terms), data = newdata)

  # TODO: support predict(..., type = "response"/"link"/other?)
  xb <- drop(X %*% model$coefficients)

  # this branch applies to (at least) `glm`, `survey::surveyglm`,
  # `robustbase::glmrob` and `gam` models
  if (inherits(model, "glm")) {
    ca <- model$family$linkinv(xb)
    pred_gradient <- model$family$mu.eta(xb) * X
  } else if (inherits(model, "lm") | inherits(model, "lmrob")) {
    # `lm` doesn't have a `family` object, but we know its prediction gradient
    ca <- xb
    pred_gradient <- X
  } else {
    stop("`model` must inherit from a `glm`, `lm`, or `lmrob` object")
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
  if (!inherits(x, "PreSandwichLayer")) {
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

  uoanames <- var_names(design, "u")
  wide_frame <- tryCatch(
    stats::expand.model.frame(x@fitted_covariance_model, uoanames, na.expand = TRUE)[uoanames],
    error = function(e) {
      stop(paste("The",
                 gsub("_", " ", design@unit_of_assignment_type),
                 "columns",
                 paste(setdiff(uoanames, colnames(covmoddata)), collapse = ", "),
                 "are missing from the covariance model dataset"),
           call. = FALSE)
    })
  keys <- .merge_preserve_order(wide_frame, design@structure, all.x = TRUE, sort = FALSE)
  keys[is.na(keys[, var_names(design, "t")]), uoanames] <- NA
  keys <- keys[, uoanames, drop = FALSE]

  return(new("SandwichLayer",
             x,
             keys = keys,
             Design = design))
}
