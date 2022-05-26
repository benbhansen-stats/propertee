#' @include DesignAccessors.R
NULL

SandwichLayer <- setClass("SandwichLayer",
                          contains = "numeric",
                          slots = c(fitted_covariance_model = "ANY",
                                    prediction_gradient = "matrix",
                                    keys = "data.frame"))

setValidity("SandwichLayer", function(object) {
  if (!("terms" %in% attr(object@fitted_covariance_model, "names"))) {
    return("Fitted covariance model must have a 'terms' attribute")
  }
  
  sandwich_env <- new.env()
  assign("bread", sandwich::bread, envir=sandwich_env)
  assign("estfun", sandwich::estfun, envir=sandwich_env)
  if (is.null(utils::getS3method("bread",
                                 class(object@fitted_covariance_model),
                                 optional = T,
                                 envir = sandwich_env)) |
      is.null(utils::getS3method("estfun",
                                 class(object@fitted_covariance_model),
                                 optional = T,
                                 envir = sandwich_env))) {
    return("Functions for extracting vcov elements not applicable to fitted covariance model")
  }
  
  if (!is.numeric(object@prediction_gradient)) {
    return("Prediction gradient must be a numeric matrix")
  }
  if (dim(object@prediction_gradient)[1] != length(object)) {
    return("Prediction gradient does not have the same number of rows as the offset")
  }
  if (dim(object@prediction_gradient)[2] != length(object@fitted_covariance_model$coefficients)) {
    return(paste0("Prediction gradient does not have the same number of columns as ",
                  "predictors in the covariance model"))
  }
  
  if (nrow(object@keys) != nrow(model.matrix(object@fitted_covariance_model))) {
    return(paste0("Keys does not have the same number of rows as the dataset used ",
                  "to fit the covariance model"))
  }

  if (any(is.na(object))) {
    warning(paste0("Offset has NA values; be careful of dropping these observations ",
                   "when fitting the design model"))
  }
  TRUE
})

PreSandwichLayer <- setClass("PreSandwichLayer",
                             contains = "numeric",
                             slots = c(fitted_covariance_model = "ANY",
                                       prediction_gradient = "matrix"))

setValidity("PreSandwichLayer", function(object) {
  if (!("terms" %in% attr(object@fitted_covariance_model, "names"))) {
    return("Fitted covariance model must have a 'terms' attribute")
  }
  
  sandwich_env <- new.env()
  assign("bread", sandwich::bread, envir=sandwich_env)
  assign("estfun", sandwich::estfun, envir=sandwich_env)
  if (is.null(utils::getS3method("bread",
                                 class(object@fitted_covariance_model),
                                 optional = T,
                                 envir = sandwich_env)) |
      is.null(utils::getS3method("estfun",
                                 class(object@fitted_covariance_model),
                                 optional = T,
                                 envir = sandwich_env))) {
    return("Functions for extracting vcov elements not applicable to fitted covariance model")
  }
  
  if (!is.numeric(object@prediction_gradient)) {
    return("Prediction gradient must be a numeric matrix")
  }
  if (dim(object@prediction_gradient)[1] != length(object)) {
    return("Prediction gradient does not have the same number of rows as the offset")
  }
  if (dim(object@prediction_gradient)[2] != length(object@fitted_covariance_model$coefficients)) {
    return(paste0("Prediction gradient does not have the same number of columns as ",
                  "predictors in the covariance model"))
  }
  TRUE
})

##' @title Convert a PreSandwichLayer to a SandwichLayer via a Design Object
##' @param x PreSandwichLayer
##' @param design Design
##' @param by character; named vector with entries corresponding to Design's
##' group assignment columns and names corresponding to the new desired column
##' names
##' @param envir Environment; the environment in which to find the covariance
##' model data to instantiate, default is `parent.frame()`
##' @return SandwichLayer
##' @export
as.SandwichLayer <- function(x, design, by = NULL, envir = parent.frame()) {
  check_desvar_cols <- function(cols, cov_mod) {
    missing_desvar_cols <- setdiff(cols,
                                   colnames(eval(cov_mod[["call"]][["data"]],
                                                 envir = envir)))
    if (length(missing_desvar_cols) > 0) {
      stop(paste0("The treatment assignment columns \"",
                  paste(missing_desvar_cols, collapse = "\", \""),
                  "\" are missing from the dataset used to fit the covariance model"))
    }
  }

  if (!is(x, "PreSandwichLayer")) {
    stop("x must be a `PreSandwichLayer` object")
  }
  if (is.null(x@fitted_covariance_model[["call"]][["data"]])) {
    stop("The fitted covariance model for x must be fit using a `data` argument")
  }

  if (is.null(by)) {
    by <- desvars <- var_names(design, "u")
  } else {
    desvars <- names(by)
    if (is.null(desvars)) {
      stop("`by` must be a named vector")
    }
  }

  check_desvar_cols(by, x@fitted_covariance_model)
  wide_frame <- stats::expand.model.frame(x@fitted_covariance_model,
                                          by,
                                          na.expand = T)[by]
  wide_frame$idx <- 1:nrow(wide_frame)  # add idx to keep merge results in place
  keys <- merge(wide_frame, design@structure, all.x = T)
  keys[is.na(keys[, var_names(design, "t")]), by] <- NA
  keys <- keys[order(keys$idx), by]
  colnames(keys) <- desvars
  row.names(keys) <- wide_frame$idx
  
  return(new("SandwichLayer",
             x@.Data,
             fitted_covariance_model = x@fitted_covariance_model,
             prediction_gradient = x@prediction_gradient,
             keys = keys))
}
