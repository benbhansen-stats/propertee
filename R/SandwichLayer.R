#' @include Design.R DesignAccessors.R
NULL

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
  if (dim(object@prediction_gradient)[2] != object@fitted_covariance_model$rank) {
    return(paste("Prediction gradient does not have the same number of columns",
                 "predictors in the covariance adjustment model"))
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
            x@prediction_gradient <- subset(x@prediction_gradient, subset = idx)
            return(x)

          })

##' (Internal) Get the a vector of "response" predictions from a covariance adjustment model
##' and its gradient with respect to the fitted coefficients
##' @param model Any model that inherits from a \code{glm}, \code{lm}, or \code{
##' robustbase::lmrob} object
##' @param newdata Optional; a data.frame of new data
##' @return Covariate adjusted outcomes and their gradient with respect to the
##' parameters of the covariance adjustment model (a list of a numeric vector and a matrix)
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
  # use the `stats` package's method for handling model fits not of full rank
  p <- model$rank
  p1 <- seq_len(p)
  piv <- if(p) model$qr$pivot[p1]
  if(p < ncol(X)) {
    warning("prediction from a rank-deficient fit may be misleading")
    X <- X[, piv, drop = FALSE]
  }

  # TODO: support predict(..., type = "response"/"link"/other?)
  xb <- drop(X %*% model$coefficients[piv])

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
##' variables in the covariance adjustment data. Names represent variables in the
##' Design; values represent variables in the data. Only needed if variable names
##' differ.
##' @param Q_data optional; quasiexperimental dataframe to allow for merging
##' to the covariance adjustment data via the `by` argument. If `by` is not NULL
##' and `Q_data` is NULL, the function will search through the call stack to find
##' the quasiexperimental data; if the search is unsuccessful, the function will
##' use the `design` data to attempt to merge.
##' @return a \code{SandwichLayer} object
##' @export
as.SandwichLayer <- function(x, design, by = NULL, Q_data = NULL) {
  if (!inherits(x, "PreSandwichLayer")) {
    stop("x must be a `PreSandwichLayer` object")
  }
  
  if (is.null(by)) {
    by <- stats::setNames(var_names(design, "u"), var_names(design, "u"))
    if (is.null(Q_data)) Q_data <- design@structure
  } else {
    if (is.null(names(by))) {
      names(by) <- by
    } else if (!is.null(names(by))) {
      names(by)[names(by) == ""] <- by[names(by) == ""]
    }

    if (is.null(Q_data)) {
      form <- .update_ca_model_formula(x@fitted_covariance_model, by, design)
      Q_data <- tryCatch(
        .get_data_from_model("cov_adj", form),
        error = function(e) {
          warning(paste("Could not find quasiexperimental data in the call stack,",
                        "or it did not contain the columns specified in `by`.",
                        "Searching for `names(by)` in `design@structure`.",
                        "Supply the quasiexperimental data to the `Q_data`",
                        "argument when using `by` to avoid this error."),
                  call. = FALSE)
          tryCatch({
            design@structure[, names(by), drop = FALSE]
          }, error = function(e) {
            stop("Could not find columns specified in `by` in design@structure",
                 call. = FALSE)
          })
          
        })
    }
    by <- c(by, stats::setNames(var_names(design, "u"), var_names(design, "u")))
    by <- by[!duplicated(by)]
  }

  data_call <- x@fitted_covariance_model$call$data
  if (is.null(data_call)) {
    stop("The fitted covariance adjustment model for x must be fit using a `data` argument")
  }

  covmoddata <- eval(data_call,
                     envir = environment(formula(x@fitted_covariance_model)))

  keys <- tryCatch(
    stats::expand.model.frame(x@fitted_covariance_model, by, na.expand = TRUE)[by],
    error = function(e) {
      stop(paste("Columns",
                 paste(setdiff(by, colnames(covmoddata)), collapse = ", "),
                 "are missing from the covariance adjustment model dataset"),
           call. = FALSE)
    })

  keys$in_Q <- apply(
    mapply(
      function(covmod_col, des_col) {
        unique(Q_data[[des_col]])[
          match(keys[[covmod_col]], unique(Q_data[[des_col]]), incomparables = NA)
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
             Design = design))
}

#' Generate a list of sanitized units of assignment from C
#' @param x a \code{SandwichLayer} object.
#' @param by Defaults to NULL, which means unit of assignment columns
#' indicated in the Design will be used to align contributions to estimating
#' equations or generate clustered covariance estimates. A non-NULL argument to
#' `by` specifies a string or character vector of column names appearing in both
#' the covariance adjustment and quasiexperimental sa,ples that should be used
#' for clustering covariance estimates.
#' @param verbose Boolean defaulting to TRUE, which will produce rather than
#' swallow any warnings about the coding of the units of assignment in the
#' covariance adjustment model data
#' @param ... arguments passed to methods
#' @return A vector of length \eqn{|C|}, where C is the covariance adjustment sample
#' @keywords internal
.sanitize_C_ids <- function(x, by = NULL, verbose = TRUE, sorted = FALSE, ...) {
  if (!inherits(x, "SandwichLayer")) {
    stop("x must be a `SandwichLayer` object")
  }
  if (is.null(by)) {
    by <- var_names(x@Design, "u")
  }

  C_ids <- tryCatch(
    x@keys[, by, drop = FALSE],
    error = function(e) {
      tryCatch({
        camod <- x@fitted_covariance_model
        stats::expand.model.frame(camod, by, na.expand = TRUE)[by]
      },
      error = function(e) {
        camod <- x@fitted_covariance_model
        data <- eval(camod$call$data, envir = environment(formula(camod)))
        stop(paste("The columns",
                   paste(setdiff(by, colnames(data)), collapse = ", "),
                   "could not be found in either the SandwichLayer's `keys`",
                   "slot or the covariance adjustment model dataset"),
             call. = FALSE)
      })
    })
  
  ## NOTE: This code has been updated because of the lingering issue that we
  ## can't cluster on partially known unit of assignment information without
  ## being provided a notion of nesting levels for the units of assignment.
  ## Commented out code in the remaining lines pertain to this issue. For the time
  ## being, we treat any observations that couldn't be fully matched to the
  ## quasiexperimental sample as independent.
  # check_nas_funcs <- list(all = all, any_not_all = function(row) any(row) & !all(row))
  check_nas_funcs <- list(any = any)
  nas <- lapply(check_nas_funcs,
                function(f) which(apply(is.na(C_ids), 1, f)))
  out <- apply(C_ids, 1, function(...) paste(..., collapse = "_"))
  names(out) <- NULL

  # warn if verbose
  if (verbose) {
    # if (length(nas[["any_not_all"]]) > 0) {
    #   warning(paste("Some rows in the covariance adjustment model dataset have",
    #                 "NA's for some but not all clustering columns. Rows sharing",
    #                 "the same non-NA cluster ID's will be clustered together.",
    #                 "If this is not intended, provide unique non-NA cluster ID's",
    #                 "for these rows."))
    # }
    # if (length(nas[["all"]]) > 0) {
    if (length(nas[["any"]]) > 0) {
      warning(paste("Some or all rows in the covariance adjustment model could",
                    # "are found to have NA's for the given clustering columns.",
                    "not be matched to the quasiexperimental sample.",
                    "This is taken to mean these observations should be treated",
                    "as independent. To avoid this warning, provide unique non-NA cluster",
                    "ID's for each row."))
    }
  }
  
  # function for creating new uoa ID's for observations in C but not Q: create
  # "_"-separated randomly generated collections of 4 upper- or lower-case letters,
  # with the number of separations given by the number of uoa columns - 1. This
  # will maintain structure if future strsplit calls are desired
  create_new_ids <- function(n_ids, n_id_cols) {
    new_ids <- split(sample(c(letters, LETTERS), n_ids * n_id_cols * 4, replace = TRUE),
                     seq_len(n_ids))
    mapply(
      function(id) {
        split_ids <- tapply(id, rep(seq_len(n_id_cols), each = 4),
                            function(...) paste(..., collapse = ""))
        paste(split_ids, collapse = "_")
      },
      new_ids,
      USE.NAMES = FALSE
    )
  }

  # since we use available unit of assignment information, we make sure to assign
  # new uoa ID's to all observations that fall in the same uoa
  # if (length(nas$any_not_all) > 0) {
  #   replace_ids <- unique(out[nas$any_not_all])
  #   n_uoa_cols <- length(cluster)
  #   new_ids <- create_new_ids(length(replace_ids), n_uoa_cols)
  #   names(new_ids) <- replace_ids
  #   out[nas$any_not_all] <- new_ids[out[nas$any_not_all]]
  # }
  # for observations with no unit of assignment information, we create unique uoa
  # ID's
  # if (length(nas$all) > 0) {
  if (length(nas[["any"]]) > 0) {
    # n_replace_ids <- length(nas$all)
    n_replace_ids <- length(nas$any)
    n_id_cols <- length(by)
    new_ids <- create_new_ids(n_replace_ids, n_id_cols)
    # out[nas$all] <- new_ids
    out[nas$any] <- new_ids
  }

  if (sorted) {
    if (suppressWarnings(any(is.na(as.numeric(out))))) {
      out <- sort(out, index.return = TRUE)
    } else {
      out <- sort(as.numeric(out), index.return = TRUE)
    }
  }
  
  return(out)
}
