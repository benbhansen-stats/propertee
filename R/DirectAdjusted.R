#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R
NULL
# The above ensures that `Design`, `WeightedDesign`, and `vcovDA` are defined
# prior to `DirectAdjusted`

setClass("DirectAdjusted",
         contains = "lm",
         slots = c(Design = "Design",
                   lmitt_fitted = "logical",
                   absorbed_intercepts = "logical",
                   absorbed_moderators = "character"))

setValidity("DirectAdjusted", function(object) {
  return(TRUE)
})

##' @title Show an DirectAdjusted
##' @param object DirectAdjusted object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "DirectAdjusted", function(object) {
  coeffs <- object$coefficients
  print(coeffs)
  invisible(object)
})

##' @title Variance-Covariance matrix of \code{DirectAdjusted} object
##' @details If a \code{DirectAdjusted} object is fit with a \code{SandwichLayer}
##' offset, then its \code{vcov()} method provides a sandwich estimate of the
##' covariance-adjusted variance-covariance matrix. Otherwise, it provides
##' the default OLS estimate of the matrix.
##' @param object DirectAdjusted
##' @param ... Additional arguments to \code{vcovDA()} or \code{stats:::vcov.lm()}.
##' @return Variance-Covariance matrix
##' @exportS3Method
vcov.DirectAdjusted <- function(object, ...) {
  cl <- match.call()

  if (is.null(cl[["type"]])) {
    confint_calls <- grepl("confint.DirectAdjusted", lapply(sys.calls(), "[[", 1))
    if (any(confint_calls)) {
      type <- tryCatch(get("call", sys.frame(which(confint_calls)[1]))$type,
                       error = function(e) NULL)
      cl$type <- type # will not append if type is NULL
    }
  }

  cl$x <- cl$object
  argmatch <- match(c("x", "type", "cluster"), names(cl), nomatch = 0L)
  new_cl <- cl[c(1L, argmatch)]
  new_cl[[1L]] <-  quote(vcovDA)

  vmat <- eval(new_cl, parent.frame())
  return(vmat)
}

##' @title Variance-Covariance matrix
##' @param object DirectAdjusted
##' @param parm a specification of which parameters are to be given confidence
##'   intervals, either a vector of numbers or a vector of names. If missing,
##'   all parameters are considered.
##' @param level the confidence level required.
##' @param ... Add'l arguments
##' @return Variance-Covariance matrix
##' @exportS3Method
confint.DirectAdjusted <- function(object, parm, level = 0.95, ...) {
  call <- match.call()
  call[[1L]] <- quote(stats::confint.lm)

  ci <- eval(call, parent.frame())
  return(ci)
}

##' @title Extract empirical estimating functions from a \code{DirectAdjusted} model fit
##' @param x a fitted \code{DirectAdjusted} object
##' @param ... arguments passed to methods
##' @return An \eqn{n\times k} matrix containing the empirical estimating
##' functions of the directly adjusted ITT model. \eqn{n} represents the
##' combined number of units at the observation level in the samples used to
##' fit the covariance adjustment and ITT effect models. \eqn{k} represents the
##' number of parameters in the ITT effect model.\cr\cr
##' Each row represents an observation's contribution to the stacked estimating
##' equations. This contribution, denoted \eqn{\tilde{\psi}_{i}} for the \eqn{i}th
##' observation, is given by \deqn{\tilde{\psi}_{i} = \psi_{i} +
##' \phi_{i}A_{11}^{-1}A_{21}^{T}} where \eqn{\psi_{i}} is the observation's
##' contribution to the ITT effect model fit, \eqn{\phi_{i}} is the observation's
##' contribution to the covariance adjustment model fit, and the \eqn{A} matrices
##' are given by typical sandwich calculations.\cr\cr
##' Note that this formulation requires a row for each observation used to fit
##' either the covariance adjustment model or the ITT effect model. The output
##' matrix is ordered such that the units used to fit the latter comprise the
##' initial rows, and those used to fit the former follow. When there is overlap
##' of the rows used to fit the two models, there is no guarantee this method
##' aligns each observation's contributions. However, given clustering information
##' in a \code{DirectAdjusted}'s \code{Design} object, the contributions can be
##' aggregated such that the resulting sandwich variance estimates are correct.
##' @exportS3Method
estfun.DirectAdjusted <- function(x, ...) {
  ## If user has passed a custom `cluster` argument to use in `sandwich::meatCL`,
  ## use it here, otherwise use the default unit of assignment columns from the
  ## Design (validity of custom `cluster` argument will be checked in )
  args <- list(...)
  if (!inherits(cluster <- args$cluster, "character")) {
    cluster <- var_names(x@Design, "u")
  }

  ## this vector indicates the hierarchy of `sandwich::estfun` methods to use
  ## to extract estimating equations for ITT model
  valid_classes <- c("glm", "lmrob", "svyglm", "lm")
  base_class <- match(x@.S3Class, valid_classes)
  if (all(is.na(base_class))) {
    stop(paste("ITT effect model must have been fitted using a function from the",
               "`flexida`, `stats`, `robustbase`, or `survey` package"))
  }
  psi <- getS3method("estfun", valid_classes[min(base_class, na.rm = TRUE)])(x)

  ## if ITT model offset doesn't contain info about covariance model, psi should
  ## be the matrix of estimating equations returned
  ca <- x$model$`(offset)`
  if (is.null(ca) | !inherits(ca, "SandwichLayer")) {
    return(psi)
  }

  ## otherwise, extract/compute the rest of the relevant matrices/quantities
  cmod <- ca@fitted_covariance_model
  phi <- estfun(cmod)
  a11_inv <- .get_a11_inverse(x)
  a21 <- .get_a21(x)

  # add rows to estimating equations based on overlap of C and Q
  uoas <- .sanitize_uoas(x, cluster, verbose = FALSE, ...)
  if (any(uoas == "C")) {
    psi <- rbind(psi, matrix(0, nrow = sum(uoas == "C"), ncol = ncol(psi)))
  }
  if (any(uoas == "Q")) {
    phi <- rbind(matrix(0, nrow = sum(uoas == "Q"), ncol = ncol(phi)), phi)
  }

  ## form matrix of estimating equations
  mat <- psi - phi %*% a11_inv %*% t(a21)

  return(mat)
}

##' @title Extract bread matrix from a \code{DirectAdjusted} model fit
##' @param x a fitted \code{DirectAdjusted} object
##' @param ... arguments passed to methods
##' @return A \eqn{k\times k} matrix where k denotes the number of parameters
##' in the ITT effect model. This corresponds to the Hessian of the ITT effect
##' model estimating equations defined in our accompanying documentation.
##' @exportS3Method
bread.DirectAdjusted <- function(x, ...) {
  if (!inherits(ca <- x$model$`(offset)`, "SandwichLayer")) {
    return(utils::getS3method("bread", "lm")(x))
  }
  if (is.null(x$qr)) {
    stop(paste("Cannot compute the Hessian of the ITT effect model estimating",
               "equations if the model fit does not have a `qr` element."))
  }

  mm <- stats::model.matrix(x)
  # compute scaling factor
  nq <- nrow(mm)
  nc_not_q <- sum(apply(is.na(ca@keys), 1, any))
  n <- nq + nc_not_q

  out <- n * chol2inv(x$qr$qr)
  dimnames(out) <- list(colnames(mm), colnames(mm))

  return(out)
}

#' Make unit of assignment ID's that align with the output of
#' \code{estfun.DirectAdjusted}
#' @details \code{estfun.DirectAdjusted} stacks the rows from Q, the
#' quasiexperimental sample, atop the rows in C that don't overlap with Q, C
#' being the covariance adjustment sample. Thus, the number of rows in the
#' estimating equations matrix is equal to \eqn{|Q| + |C \ Q|}, so
#' \code{.make_uoa_ids} returns a vector of that length with corresponding
#' units of assignment.
#' @param x a fitted \code{DirectAdjusted} object
#' @param cluster Defaults to NULL, which means unit of assignment columns
#' indicated in the Design will be used to generate clustered covariance estimates.
#' A non-NULL argument to `cluster` specifies a string or character vector of
#' column names appearing in both the covariance adjustment and quasiexperimental
#' samples that should be used for clustering covariance estimates.
#' @param ... arguments passed to methods
#' @return A vector of length \eqn{|Q| + |C \ Q|}
#' @keywords internal
.make_uoa_ids <- function(x, cluster = NULL, ...) {
  if (!inherits(cluster, "character") & !inherits(x, "DirectAdjusted")) {
    stop(paste("Cannot deduce units of assignment for clustering without a",
               "Design object (stored in a DirectAdjusted object) or a `cluster`",
               "argument specifying the columns with the units of assignment"))
  }

  # Must be a DirectAdjusted object for this logic to occur
  if (!inherits(cluster, "character")) {
    cluster <- var_names(x@Design, "u")
  }

  # If there's no covariance adjustment info, return the ID's found in Q
  if (!inherits(x, "DirectAdjusted") | !inherits(x$model$`(offset)`, "SandwichLayer")) {
    Q_uoas <- .sanitize_Q_uoas(x, cluster, ...)
    return(factor(Q_uoas, levels = unique(Q_uoas)))
  }

  uoas <- .sanitize_uoas(x, cluster, ...)

  return(factor(names(uoas), levels = unique(names(uoas))))
}

#' Generate a vector of sanitized units of assignment from C and Q
#' @param x a fitted \code{DirectAdjusted} object
#' @param cluster Defaults to NULL, which means unit of assignment columns
#' indicated in the Design will be used to generate clustered covariance estimates.
#' A non-NULL argument to `cluster` specifies a string or character vector of
#' column names appearing in both the covariance adjustment and quasiexperimental
#' samples that should be used for clustering covariance estimates.
#' @param ca SandwichLayer object storing information about the covariance
#' adjustment model; usually stored as the `offset` of a \code{DirectAdjusted}
#' object when covariance adjustment is performed
#' @param verbose Boolean defaulting to TRUE, which will produce rather than
#' swallow any warnings about the coding of the units of assignment in the
#' covariance adjustment model data
#' @return A named vector of length \eqn{|Q| + |C \ Q|}, where Q and C represent
#' the sets of observations in the quasiexperimental sample and covariance
#' adjustment sample, respectively. Values can be "Q", for observations that only
#' appear in Q, "C", for those that only appear in C, and "QC" for those that
#' appear in both. Names correspond to the unit of assignment ID.
#' @param ... arguments passed to methods
#' @keywords internal
.sanitize_uoas <- function(x, cluster = NULL, ca = x$model$`(offset)`, verbose = TRUE, ...) {
  if (!inherits(x, "DirectAdjusted") | !inherits(ca, "SandwichLayer")) {
    stop(paste("x must be a DirectAdjusted object with a SandwichLayer offset or",
               "ca must be a SandwichLayer object to retrieve information about",
               "the covariance adjustment model"))
  }
  if (is.null(cluster)) {
    cluster <- var_names(x@Design, "u")
  }

  Q_uoas <- .sanitize_Q_uoas(x, cluster, ...)
  C_uoas <- .sanitize_C_uoas(ca, cluster, verbose = verbose, ...)

  # return a vector of the sample the observations pertain to named by their
  # unit of assignment
  not_Q_idx <- !(C_uoas %in% unique(Q_uoas))
  uoas <- vector("character", length(Q_uoas) + sum(not_Q_idx))
  uoas[which(!(Q_uoas %in% unique(C_uoas)))] <- "Q"
  uoas[which(Q_uoas %in% unique(C_uoas))] <- "Q_C"
  uoas[which(uoas == "")] <- "C"
  names(uoas) <- c(Q_uoas, C_uoas[not_Q_idx])

  return(uoas)
}

#' Generate a list of sanitized units of assignment from Q
#' @param x a fitted \code{DirectAdjusted} object
#' @param cluster Defaults to NULL, which means unit of assignment columns
#' indicated in the Design will be used to generate clustered covariance estimates.
#' A non-NULL argument to `cluster` specifies a string or character vector of
#' column names appearing in both the covariance adjustment and quasiexperimental
#' samples that should be used for clustering covariance estimates.
#' @param ... arguments passed to methods
#' @return A vector of length \eqn{|Q|}, where Q is the quasiexperimental sample
#' @keywords internal
.sanitize_Q_uoas <- function(x, cluster = NULL, ...) {
  if (is.null(cluster)) {
    cluster <- var_names(x@Design, "u")
  }
  uoas <- tryCatch(
    stats::expand.model.frame(x, cluster, na.expand = TRUE)[, cluster, drop = FALSE],
    error = function(e) {
      mf <- eval(x$call$data, envir = environment(x))
      missing_cols <- setdiff(cluster, colnames(mf))
      stop(paste("Could not find unit of assignment columns",
                 paste(missing_cols, collapse = ", "), "in ITT effect model data"),
           call. = FALSE)
    })
  out <- apply(uoas, 1, function(...) paste(..., collapse = "_"))
  names(out) <- NULL

  return(out)
}

#' @exportS3Method getCall DirectAdjusted
#' @importFrom stats getCall
getCall.DirectAdjusted <- function(x, ...) {
  return(x$call)
}

#' @export
update.DirectAdjusted <- function(object, ...) {
  stop(paste("DirectAdjusted objects do not support an `update` method.\n",
             "You can use `update` on the formula object passed to `lmitt`",
             "instead."))
}
