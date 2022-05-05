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
  if (dim(object@prediction_gradient)[2] != length(coef(object@fitted_covariance_model))) {
    return(paste0("Prediction gradient does not have the same number of columns as ",
                  "predictors in the covariance model"))
  }
  
  if (is.null(dim(object@keys))) {
    return("Keys must be a valid dataframe")
  }
  if (nrow(object@keys) != nrow(model.matrix(object@fitted_covariance_model))) {
    return("Keys does not have the same number of rows as the experiment design matrix")
  }
  
  if (any(is.na(object))) {
    warning("Offset has NA values; these observations will be dropped in the design model")
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
  if (dim(object@prediction_gradient)[2] != length(coef(object@fitted_covariance_model))) {
    return(paste0("Prediction gradient does not have the same number of columns as ",
                  "predictors in the covariance model"))
  }
  TRUE
})

setMethod("as.SandwichLayer", "PreSandwichLayer", function(object, design, by=NULL) {
  if (is.null(by)) {
    desvars <- colnames(design@structure)
    missing_desvars <- setdiff(desvars,
                               colnames(eval(object@fitted_covariance_model$call$data)))
    for (col in missing_desvars) {
      object@fitted_covariance_model$call$data[col] <- NA_integer_
    }
    
    keys <- expand.model.frame(object@fitted_covariance_model,
                               desvars,
                               na.expand = T)[colnames(design@structure)]
  } else {
    if (is.null(names(by))) {
      stop("`by` must be a named vector")
    }
    
    desvars <- names(by)
    keys <- expand.model.frame(object@fitted_covariance_model,
                               by[desvars],
                               na.expand = T)[by[desvars]]
    colnames(keys) <- desvars
  }
  
  return(SandwichLayer(object@.Data,
                       fitted_covariance_model = object@fitted_covariance_model,
                       prediction_gradient = object@prediction_gradient,
                       keys = keys))
})
