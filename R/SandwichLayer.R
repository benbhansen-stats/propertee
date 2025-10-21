#' @include StudySpecification.R StudySpecificationAccessors.R
NULL

#' (Internal) model predictions with some model artifacts, as S4 object
#' @slot .Data numeric vector of predictions
#' @slot fitted_covariance_model model standing behind the predictions
#' @slot prediction_gradient matrix, predictions gradient w/r/t model params
#' @keywords internal
setClass("PreSandwichLayer",
         contains = "numeric",
         slots = c(fitted_covariance_model = "ANY",
                   prediction_gradient = "matrix"))

setValidity("PreSandwichLayer", function(object) {
  if (!inherits(object@fitted_covariance_model$terms, "terms")) {
    return("Fitted covariance adjustment model must have a valid 'terms' attribute")
  }

  tryCatch({
    sandwich::bread(object@fitted_covariance_model)
    sandwich::estfun(object@fitted_covariance_model)
    }, error = function(e) {
      stop(paste("Functions for extracting vcov elements not applicable to",
                 "fitted covariance adjustment model"))
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
      qr(stats::model.matrix(object@fitted_covariance_model))$rank) {
    return(paste("Prediction gradient does not have the same number of columns",
                 "predictors in the covariance adjustment model"))
  }
  TRUE
})

#' (Internal) model predictions with more model artifacts, as S4 object
#'
#' Contains PreSandwichLayer class.  Only additional slots listed here.
#' @slot keys a data.frame
#' @slot StudySpecification a StudySpecification
#' @keywords internal
setClass("SandwichLayer",
         contains = "PreSandwichLayer",
         slots = c(keys = "data.frame",
                   StudySpecification = "StudySpecification"))

setValidity("SandwichLayer", function(object) {
  if (nrow(object@keys) != nrow(stats::model.matrix(object@fitted_covariance_model))) {
    return(paste("Keys does not have the same number of rows as the dataset",
                 "used to fit the covariance adjustment model"))
  }

  if (any(is.na(object))) {
    msg <- paste("Some covariance adjustments are NA; be careful of dropping",
                 "these observations when fitting the ITT effect model")
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

##' @title Show a \code{PreSandwichLayer} or \code{SandwichLayer}
##'
##' @description Display information about a \code{PreSandwichLayer} or
##'   \code{SandwichLayer} object
##'
##' @param object \code{PreSandwichLayer} or \code{SandwichLayer} object
##' @return an invisible copy of \code{object}
##' @export
##' @rdname PreSandwichLayer.show
setMethod("show", "PreSandwichLayer", .show_layer)

setGeneric("subset")

##' @title \code{PreSandwichLayer} and \code{SandwichLayer} subsetting
##'
##' @description Return subset of a \code{PreSandwichLayer} or
##'   \code{SandwichLayer} which meets conditions.
##'
##' @param subset Logical vector identifying values to keep or drop
##' @return \code{x} subset by \code{subset} or \code{i}
##' @export
##' @rdname PreSandwichLayer.subset
setMethod("subset", "PreSandwichLayer", function(x, subset) {
  x@.Data <- subset(x@.Data, subset = subset)
  x@prediction_gradient <- subset(x@prediction_gradient, subset = subset)
  return(x)
})

setGeneric("[")

##' @param x \code{PreSandwichLayer} or \code{SandwichLayer} object
##' @param i indices specifying elements to extract or replace. See
##'   \code{help("[")} for further details.
##' @export
##' @importFrom methods callNextMethod
##' @rdname PreSandwichLayer.subset
setMethod("[", "PreSandwichLayer",
          function(x, i) {
            idx <- rep(FALSE, length(x@.Data))
            idx[i] <- TRUE
            dat <- methods::callNextMethod()
            x@.Data <- dat
            x@prediction_gradient <- subset(x@prediction_gradient, subset = idx)
            return(x)

          })

##' @title (Internal) Get covariance adjustments and their gradient with respect
##'   to covariance adjustment model coefficients
##' @description \code{.make_PreSandwichLayer()} takes a fitted covariance
##'   adjustment model passed to the \code{model} argument and generates
##'   adjustments to outcomes for observations in the \code{newdata} argument.
##'   It also evaluates the gradient of the adjustments taken with respect to
##'   the coefficients at the coefficient estimates.
##' @param model a fitted model to use for generating covariance adjustment
##'   values
##' @param newdata a dataframe with columns called for in \code{model}
##' @param ... additional arguments to pass on to \code{model.frame} and
##'   \code{model.matrix}. These cannot include \code{na.action}, \code{xlev},
##'   or \code{contrasts.arg}: the former is fixed to be \code{na.pass}, while
##'   the latter two are provided by elements of the \code{model} argument.
##' @return A \code{PreSandwichLayer} object
##' @keywords internal
.make_PreSandwichLayer <- function(model, newdata = NULL, ...) {
  UseMethod(".make_PreSandwichLayer")
}

##' @keywords internal
.make_PreSandwichLayer.default <- function(model, newdata = NULL, ...) {
  model_terms <- tryCatch(
    terms(model),
    error = function(e) stop("`model` must have `terms` method", call. = FALSE)
  )

  if (is.null(newdata)) {
    mf <- .get_data_from_model("cov_adj", formula(model))
  } else if (!is.null(term <- attr(newdata, "terms"))) {
    mf <- newdata
  } else {
    mf <- stats::model.frame(stats::delete.response(terms(model)),
                             data = newdata,
                             na.action = na.pass,
                             xlev = model$xlevels,
                             ...)
  }

  X <- stats::model.matrix(stats::delete.response(terms(model)),
                           data = mf,
                           contrasts.arg = model$contrasts, # contrasts.arg also follows `predict.lm`
                           ...)

  # use the `stats` package's method for handling model fits not of full rank
  p <- model$rank
  p1 <- seq_len(p)
  piv <- if(p) model$qr$pivot[p1]
  if(p < ncol(X)) {
    X <- X[, piv, drop = FALSE]
  }

  # TODO: support predict(..., type = "response"/"link"/other?)
  xb <- drop(X %*% model$coefficients[piv])

  return(new("PreSandwichLayer",
             xb,
             fitted_covariance_model = model,
             prediction_gradient = X))
}
.S3method(".make_PreSandwichLayer", "default", .make_PreSandwichLayer.default)

##' @keywords internal
.make_PreSandwichLayer.glm <- function(model, newdata = NULL, ...) {
  model_terms <- tryCatch(
    terms(model),
    error = function(e) stop("`model` must have `terms` method", call. = FALSE)
  )

  if (is.null(newdata)) {
    mf <- .get_data_from_model("cov_adj", formula(model))
  } else if (!is.null(term <- attr(newdata, "terms"))) {
    mf <- newdata
  } else {
    mf <- stats::model.frame(stats::delete.response(terms(model)),
                             data = newdata,
                             na.action = na.pass,
                             xlev = model$xlevels,
                             ...)
  }

  X <- stats::model.matrix(stats::delete.response(terms(model)),
                           data = mf,
                           contrasts.arg = model$contrasts,
                           ...)

  p <- model$rank
  p1 <- seq_len(p)
  piv <- if(p) model$qr$pivot[p1]
  if(p < ncol(X)) {
    X <- X[, piv, drop = FALSE]
  }

  # TODO: support predict(..., type = "response"/"link"/other?)
  xb <- drop(X %*% model$coefficients[piv])

  return(new("PreSandwichLayer",
             model$family$linkinv(xb),
             fitted_covariance_model = model,
             prediction_gradient = model$family$mu.eta(xb) * X))
}
.S3method(".make_PreSandwichLayer", "glm", .make_PreSandwichLayer.glm)

##' @keywords internal
.make_PreSandwichLayer.glmrob <- function(model, newdata = NULL, ...) {
  model_terms <- tryCatch(
    terms(model),
    error = function(e) stop("`model` must have `terms` method", call. = FALSE)
  )

  if (is.null(newdata)) {
    mf <- .get_data_from_model("cov_adj", formula(model))
  } else if (!is.null(term <- attr(newdata, "terms"))) {
    mf <- newdata
  } else {
    mf <- stats::model.frame(stats::delete.response(terms(model)),
                             data = newdata,
                             na.action = na.pass,
                             xlev = model$xlevels,
                             ...)
  }

  X <- stats::model.matrix(stats::delete.response(terms(model)),
                           data = mf,
                           contrasts.arg = model$contrasts,
                           ...)

  QR <- qr(model$w.r * stats::model.matrix(model))
  p <- QR$rank
  p1 <- seq_len(p)
  piv <- if(p) QR$pivot[p1]
  if(p < ncol(X)) {
    X <- X[, piv, drop = FALSE]
  }

  # TODO: support predict(..., type = "response"/"link"/other?)
  xb <- drop(X %*% model$coefficients[piv])

  return(new("PreSandwichLayer",
             model$family$linkinv(xb),
             fitted_covariance_model = model,
             prediction_gradient = model$family$mu.eta(xb) * X))
}
.S3method(".make_PreSandwichLayer", "glmrob", .make_PreSandwichLayer.glmrob)


##' @title Convert a \code{PreSandwichLayer} to a \code{SandwichLayer} with a
##'   \code{StudySpecification} object
##' @description \code{as.SandwichLayer()} uses the \code{StudySpecification}
##'   object passed to the \code{specification} argument to populate the slots
##'   in a \code{SandwichLayer} object that a \code{PreSandwichLayer} does not
##'   have sufficient information for.
##' @param x a \code{PreSandwichLayer} object
##' @param specification a \code{StudySpecification} object
##' @param by optional; a string or named vector of unique identifier columns in
##'   the data used to create \code{specification} and the data used to fit the
##'   covariance adjustment model. Default is NULL, in which case unit of
##'   assignment columns are used for identification (even if they do not
##'   uniquely identify units of observation). If a named vector is provided,
##'   names should represent variables in the data used to create
##'   \code{specification}, while values should represent variables in the
##'   covariance adjustment data.
##' @param Q_data dataframe of direct adjustment sample, which is needed to
##'   generate the \code{keys} slot of the \code{SandwichLayer} object. Defaults
##'   to NULL, in which case if \code{by} is NULL, the data used to create
##'   \code{specification} is used, and if \code{by} is not NULL, appropriate
##'   data further up the call stack (passed as arguments to \code{cov_adj()} or
##'   \code{lmitt.formula()}, for example) is used.
##' @return a \code{SandwichLayer} object
##' @export
as.SandwichLayer <- function(x, specification, by = NULL, Q_data = NULL) {
  if (!inherits(x, "PreSandwichLayer")) {
    stop("x must be a `PreSandwichLayer` object")
  }

  if (is.null(by)) {
    by <- stats::setNames(var_names(specification, "u"), var_names(specification, "u"))
    if (is.null(Q_data)) Q_data <- specification@structure
  } else {
    if (is.null(names(by))) {
      names(by) <- by
    } else if (!is.null(names(by))) {
      names(by)[names(by) == ""] <- by[names(by) == ""]
    }

    if (is.null(Q_data)) {
      form <- .update_ca_model_formula(x@fitted_covariance_model, by, specification)
      Q_data <- tryCatch(
        .get_data_from_model("cov_adj", form),
        error = function(e) {
          warning(paste("Could not find direct adjustment data in the call stack,",
                        "or it did not contain the columns specified in `by`.",
                        "Searching for `names(by)` in `specification@structure`.",
                        "Supply the direct adjustment data to the `Q_data`",
                        "argument when using `by` to avoid this error."),
                  call. = FALSE)
          tryCatch({
            specification@structure[, names(by), drop = FALSE]
          }, error = function(e) {
            stop("Could not find columns specified in `by` in specification@structure",
                 call. = FALSE)
          })

        })
    }
  }

  data_call <- x@fitted_covariance_model$call$data
  if (is.null(data_call)) {
    stop("The fitted covariance adjustment model for x must be fit using a `data` argument")
  }

  if (specification@unit_of_assignment_type == "none") {
    keys <- data.frame(..uoa.. = rownames(x@fitted_covariance_model$model))
  } else {
    keys <- tryCatch(
      stats::expand.model.frame(x@fitted_covariance_model, by, na.expand = TRUE)[by],
      error = function(e) {
        covmoddata <- eval(data_call,
                           envir = environment(formula(x@fitted_covariance_model)))
        stop(paste("Columns",
                   paste(setdiff(by, colnames(covmoddata)), collapse = ", "),
                   "are missing from the covariance adjustment model dataset"),
           call. = FALSE)
      })
  }

  keys$in_Q <- apply(
    mapply(
      function(covmod_col, spec_col) {
        unique(Q_data[[spec_col]])[
          match(keys[[covmod_col]], unique(Q_data[[spec_col]]), incomparables = NA)
        ]
      },
      by, names(by)
    ),
    1,
    function(uids) all(!is.na(uids))
  )

  return(new("SandwichLayer",
             x,
             keys = keys,
             StudySpecification = specification))
}

#' @title (Internal) Return ID's used to order observations in the covariance
#'   adjustment sample
#' @param x a \code{SandwichLayer} object
#' @param id_col character vector or list; optional. Specifies column names that
#'   appear in botn the covariance adjustment and direct adjustmet samples.
#'   Defaults to NULL, in which case unit of assignment columns in the
#'   \code{SandwichLayer}'s \code{StudySpecification} slot will be used to
#'   generate ID's.
#' @param ... arguments passed to methods
#' @return A vector of length equal to the number of units of observation used
#'   to fit the covariance adjustment model
#' @keywords internal
.sanitize_C_ids <- function(x, id_col = NULL, sorted = FALSE, ...) {
  if (!inherits(x, "SandwichLayer")) {
    stop("x must be a `SandwichLayer` object")
  }
  if (is.null(id_col)) {
    id_col <- var_names(x@StudySpecification, "u")
  }

  C_ids <- tryCatch(
    x@keys[, id_col, drop = FALSE],
    error = function(e) {
      tryCatch({
        camod <- x@fitted_covariance_model
        stats::expand.model.frame(camod, id_col, na.expand = TRUE)[id_col]
      },
      error = function(e) {
        camod <- x@fitted_covariance_model
        data <- eval(camod$call$data, envir = environment(formula(camod)))
        stop(paste("The columns",
                   paste(setdiff(id_col, colnames(data)), collapse = ", "),
                   "could not be found in either the SandwichLayer's `keys`",
                   "slot or the covariance adjustment model dataset"),
             call. = FALSE)
      })
    })

  nas <- apply(is.na(C_ids), 1, all)
  out <- apply(C_ids, 1, function(...) paste(..., collapse = "_"))
  out[nas] <- NA_character_
  names(out) <- NULL

  if (sorted) {
    if (suppressWarnings(any(is.na(as.numeric(out))))) {
      out <- sort(out, index.return = TRUE)
    } else {
      out <- sort(as.numeric(out), index.return = TRUE)
    }
  }

  return(out)
}
