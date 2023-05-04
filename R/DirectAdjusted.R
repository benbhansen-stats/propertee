#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R
NULL
# The above ensures that `Design`, `WeightedDesign`, and `vcovDA` are defined
# prior to `DirectAdjusted`

setClass("DirectAdjusted",
         contains = "lm",
         slots = c(Design = "Design",
                   lmitt_fitted = "logical",
                   absorbed_intercepts = "logical",
                   absorbed_moderators = "character",
                   lmitt_call = "call"))

setValidity("DirectAdjusted", function(object) {
  if (length(object@lmitt_fitted) != 1) {
    return("@lmitt_fitted slot must be a single logical")
  }
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

##' @title Extract empirical estimating equations from a \code{DirectAdjusted} model fit
##' @param x a fitted \code{DirectAdjusted} object
##' @param ... arguments passed to methods
##' @return An \eqn{n\times k} matrix of empirical estimating equations for the
##' covariance-adjusted ITT effect regression. \eqn{n} represents the number of
##' observations in the union of the samples used to fit the two regressions.
##' \eqn{k} represents the number of parameters in the latter model.\cr\cr
##' Each row represents an observation's contribution to the stacked estimating
##' equations. This contribution, denoted \eqn{\tilde{\psi}_{i}} for the \eqn{i}th
##' observation, is given by \deqn{\tilde{\psi}_{i} = \psi_{i} +
##' \phi_{i}A_{11}^{-1}A_{21}^{T}} where \eqn{\psi_{i}} is the observation's
##' contribution to the ITT effect model fit, \eqn{\phi_{i}} is the observation's
##' contribution to the covariance adjustment model fit, and the \eqn{A} matrices
##' are given by typical sandwich calculations. The output matrix is orderded
##' such that the units used to fit the latter comprise the initial rows, and
##' any additional observations used to fit the former are stacked below.\cr\cr
##' Note that the formulation of the output matrix \eqn{\tilde{\Psi}} requires
##' information about each observation's contributions to both the covariance
##' adjustment and ITT effect models (where some observations may not contribute
##' to both models). Estimating equations are taken from \code{sandwich::estfun}
##' calls on both fitted models and aligned as closely as possible.
##' The `by` argument in `cov_adj()` can be used to specify a column unrelated to
##' the design that will allow for exact alignment of such matrices. If no `by`
##' argument is provided, clustering information given in the \code{DirectAdjusted}'s
##' \code{Design} object will be used to align rows by unit of assignment, even
##' though no guarantees can be made about aligning the matrices within units
##' of assignment. Regardless of the eventual alignment and initial ordering of the
##' observations in the two matrices, however, when using \code{vcovDA},
##' variance estimates will ultimately be the same due to the clustering passed
##' to any \code{sandwich::meatCL} calls.
##' @exportS3Method
estfun.DirectAdjusted <- function(x, ...) {
  ## if ITT model offset doesn't contain info about covariance model, estimating
  ## equations should be the ITT model estimating equations
  if (is.null(x$model$`(offset)`) | !inherits(x$model$`(offset)`, "SandwichLayer")) {
    return(.base_S3class_estfun(x))
  }

  ## otherwise, extract/compute the rest of the relevant matrices/quantities
  estmats <- .align_and_extend_estfuns(x, ...)
  a11_inv <- .get_a11_inverse(x)
  a21 <- .get_a21(x)

  ## form matrix of estimating equations
  mat <- estmats[["psi"]] - estmats[["phi"]] %*% t(a11_inv) %*% t(a21)

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
  nc_not_q <- sum(!ca@keys$in_Q)
  n <- nq + nc_not_q

  out <- n * chol2inv(x$qr$qr)
  dimnames(out) <- list(colnames(mm), colnames(mm))

  return(out)
}

##' (Internal) Align the dimensions and rows of estimating equations matrices
##' from the ITT effect and covariance adjustment models
##' @inheritParams estfun.DirectAdjusted
##' @return list of two matrices, one being the aligned contributions to the
##' estimating equations for the ITT effect model, and the other being the
##' aligned contributions to the covariance adjustment model
##' @keywords internal
.align_and_extend_estfuns <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted") | !inherits(x$model$`(offset)`, "SandwichLayer")) {
    stop("`x` must be a fitted DirectAdjusted object with a SandwichLayer offset")
  }

  # get the original estimating equations
  psi <- .base_S3class_estfun(x)
  phi <- estfun(x$model$`(offset)`@fitted_covariance_model)

  # define the row ordering using .order_samples() and insert rows of 0's into
  # the matrices where necessary while maintaining observation alignment
  ids <- .order_samples(x, verbose = FALSE)
  aligned_psi <- rbind(
    psi[ids$Q_order, , drop = FALSE],
    matrix(0, nrow = sum(!(names(ids$C_order) %in% names(ids$Q_order))), ncol = ncol(psi))
  )
  rownames(aligned_psi) <- NULL
  aligned_phi <- matrix(0, nrow = nrow(aligned_psi), ncol = ncol(phi),
                        dimnames = list(seq_len(nrow(aligned_psi)), colnames(phi)))
  aligned_phi[which(!is.na(ids$C_order)),] <- phi[ids$C_order[!is.na(ids$C_order)], ]

  return(list(psi = aligned_psi, phi = aligned_phi))
}

##' (Internal) Call \code{sandwich::estfun} method for a fitted \code{DirectAdjusted}
##' object based on its base S3 class
##' @inheritParams estfun.DirectAdjusted
##' @return S3 method
##' @keywords internal
.base_S3class_estfun <- function(x) {
  ## this vector indicates the hierarchy of `sandwich::estfun` methods to use
  ## to extract ITT model's estimating equations
  valid_classes <- c("glm", "lmrob", "svyglm", "lm")
  base_class <- match(x@.S3Class, valid_classes)
  if (all(is.na(base_class))) {
    stop(paste("ITT effect model must have been fitted using a function from the",
               "`flexida`, `stats`, `robustbase`, or `survey` package"))
  }
  return(getS3method("estfun", valid_classes[min(base_class, na.rm = TRUE)])(x))
}

#' Make unit of assignment ID's to pass to \code{sandwich::meatCL} `cluster`
#' argument
#' @details These ID's should align with the output of \code{estfun.DirectAdjusted},
#' which stacks the rows from Q atop the rows in C that don't overlap with Q.
#' @param x a fitted \code{DirectAdjusted} object
#' @param cluster character vector or list; optional. Specifies column names that appear
#' in both the covariance adjustment dataframe C and quasiexperimental dataframe
#' Q. Defaults to NULL, in which case unit of assignment columns indicated in
#' the Design will be used to generate clustered covariance estimates.
#' @param ... arguments passed to methods
#' @return A vector of length \eqn{|Q| + |C} \ \eqn{Q|}
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
    Q_uoas <- .sanitize_Q_ids(x, cluster, sorted = FALSE, ...)
    return(factor(Q_uoas, levels = unique(Q_uoas)))
  }

  ids <- .order_samples(x, ...)
  Q_uoas <- .sanitize_Q_ids(x, cluster, ...)[ids$Q_order]
  C_uoas <- .sanitize_C_ids(x$model$`(offset)`, cluster, verbose = FALSE, ...)[
    ids$C_order[!is.na(ids$C_order)]]
  uoas <- c(Q_uoas, C_uoas[!(C_uoas %in% Q_uoas)])

  return(factor(uoas, levels = unique(uoas)))
}

#' Order observations used to fit a \code{DirectAdjusted} model and its
#' covariance adjustment model
#' @details \code{.order_samples} underpins the ordering for both \code{.make_uoa_ids}
#' and \code{estfun.DirectAdjusted}, which need to be aligned for proper
#' clustering to occur in \code{vcovDA} calls. Since \code{estfun.DirectAdjusted}
#' returns a matrix with a row count equal to \eqn{|Q| + |C} \ \eqn{Q|},
#' \code{.order_samples} must not only order the rows in Q and C, but also
#' provide information about which observations appear in both samples. How this
#' manifests is explained below.
#' @param x a fitted \code{DirectAdjusted} object
#' @param verbose boolean; optional. Defaults to TRUE, which will produce rather
#' than swallow any warnings about the coding of the units of assignment in the
#' covariance adjustment model data.
#' @return a list of two named vectors. The first element `Q_order` corresponds
#' to the order of the observations in Q with names corresponding to their ID's.
#' The second, `C_order`, corresponds to the order of observations in C but with
#' any observations in Q that do not appear in C appended as NA's to the front
#' (names still correspond to the observations' ID's). This vector has length
#' \eqn{|Q| + |C} \ \eqn{Q|}.
#' @param ... arguments passed to methods
#' @keywords internal
.order_samples <- function(x, verbose = TRUE, ...) {
  if (!inherits(x, "DirectAdjusted") | !inherits(ca <- x$model$`(offset)`, "SandwichLayer")) {
    stop(paste("x must be a DirectAdjusted object with a SandwichLayer offset or",
               "ca must be a SandwichLayer object to retrieve information about",
               "the covariance adjustment model"))
  }
  ## `keys` may have additional columns beyond "in_Q" and the uoa columns if a
  ## `by` argument was specified in `cov_adj()` or `as.SandwichLayer()`. If it
  ## does, use the columns exclusively specified in `by` to merge.
  by <- setdiff(colnames(ca@keys), "in_Q")
  if (length(setdiff(by, var_names(x@Design, "u"))) > 0) {
    by <- setdiff(by, var_names(x@Design, "u"))
  }

  # first sort the ID's in Q
  Q_ids <- .sanitize_Q_ids(x, by, sorted = TRUE, ...)
  Q_ix <- Q_ids$ix
  Q_ix <- stats::setNames(Q_ix, Q_ids$x)

  # find the overlap of C and Q and dedupe the indices returned by `match()`
  C_ids <- .sanitize_C_ids(ca, by, verbose = verbose, sorted = TRUE, ...)
  Q_C_ix <- match(Q_ids$x, C_ids$x)
  for (start_loc in unique(Q_C_ix[!is.na(Q_C_ix)])) {
    Q_C_ix[!is.na(Q_C_ix) & Q_C_ix == start_loc] <- start_loc - 1 + seq_along(
      Q_C_ix[!is.na(Q_C_ix) & Q_C_ix == start_loc])
  }
  Q_C_ix <- C_ids$ix[Q_C_ix]

  # append the remaining ID's in C if necessary
  C_ix <- C_ids$ix[!(C_ids$x %in% Q_ids$x)]
  C_ix <- c(Q_C_ix, C_ix)
  C_ix <- stats::setNames(C_ix, c(Q_ids$x, C_ids$x[!(C_ids$x %in% Q_ids$x)]))

  return(list(Q_order = Q_ix, C_order = C_ix))
}

#' Return ID's for observations in the quasiexperimental sample Q
#' @param x a fitted \code{DirectAdjusted} object
#' @param by character vector or list; optional. Specifies column names that appear in
#' botn the covariance adjustment dataframe C and quasiexperimental dataframe Q.
#' Defaults to NULL, in which case unit of assignment columns indicated in the
#' Design will be used to generate ID's.
#' @param sorted boolean defaulting to FALSE, which provides ID's in the same
#' order as the model frame. If TRUE, ID's will be sorted alphanumerically and
#' returned as a list, as given by \code{sort.int} with `index.return = TRUE`.
#' @param ... arguments passed to methods
#' @return If not `sorted`, a vector of length \eqn{|Q|}, where Q is the
#' quasiexperimental sample. If `sorted`, a list whose elements are vectors of
#' length \eqn{\{|i : i \in Q|\}}. The sorted output will be used to align observations'
#' contributions to the ITT effect model with their contributions to the covariance
#' adjustment model in \code{estfun.DirectAdjusted}.
#' @keywords internal
.sanitize_Q_ids <- function(x, by = NULL, sorted = FALSE, ...) {
  if (is.null(by)) {
    by <- var_names(x@Design, "u")
  }
  ids <- tryCatch(
    stats::expand.model.frame(x, by, na.expand = TRUE)[, by, drop = FALSE],
    error = function(e) {
      mf <- eval(x$call$data, envir = environment(x))
      missing_cols <- setdiff(by, colnames(mf))
      stop(paste("Could not find unit of assignment columns",
                 paste(missing_cols, collapse = ", "), "in ITT effect model data"),
           call. = FALSE)
    })
  out <- apply(ids, 1, function(...) paste(..., collapse = "_"))
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
