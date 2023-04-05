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
  C_uoas <- ca@keys
  phi <- estfun(cmod)
  uoa_cols <- colnames(C_uoas)
  a11_inv <- .get_a11_inverse(x)
  a21 <- .get_a21(x)

  ## figure out if rows need to be added to the matrix of estimating equations
  Q_uoas <- stats::expand.model.frame(x, uoa_cols, na.expand = TRUE)[, uoa_cols, drop = FALSE]
  Q_uoas <- apply(Q_uoas, 1, function(...) paste(..., collapse = "_"))

  C_uoas <- apply(C_uoas, 1, function(...) paste(..., collapse = "_"))
  C_uoas[vapply(strsplit(C_uoas, "_"), function(x) all(x == "NA"), logical(1))] <- NA_character_

  nas <- is.na(C_uoas)
  if (any(nas)) {
    # give unique ID's to uoa's in C but not Q
    n_Q_uoas <- length(unique(Q_uoas))
    C_uoas[nas] <- paste0(n_Q_uoas + seq_len(sum(nas)), "*")
  }

  # add rows if necessary
  add_C_uoas <- setdiff(unique(C_uoas), unique(Q_uoas))
  add_Q_uoas <- setdiff(unique(Q_uoas), unique(C_uoas))
  if (length(add_C_uoas) > 0) {
    psi <- rbind(psi, matrix(0, nrow = sum(C_uoas %in% add_C_uoas), ncol = ncol(psi)))
  }
  if (length(add_Q_uoas) > 0) {
    phi <- rbind(matrix(0, nrow = sum(Q_uoas %in% add_Q_uoas), ncol = ncol(phi)), phi)
  }

  ## form matrix of estimating equations
  mat <- psi - phi %*% a11_inv %*% t(a21)

  return(mat)
}

##' @title Extract bread matrix from a \code{DirectAdjusted} model fit
##' @param x a fitted \code{DirectAdjusted} object
##' @return A \eqn{k\times k} matrix where k denotes the number of parameters
##' in the ITT effect model. This corresponds to the Hessian of the ITT effect
##' model estimating equations defined in our accompanying documentation.
##' @exportS3Method
bread.DirectAdjusted <- function(x) {
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

##' (Internal) Obtain which variaiton of \code{assigned()}, \code{a.()} or
##' \code{z.()} is used in the model
##' @title Treatment variable name
##' @param object \code{DirectAdjusted} object
##' @return Character string identifying the name (e.g "assigned()" or
##'   "a.(des)")
##' @keywords internal
.txt_fn <- function(object) {
  cs <- names(object$coefficients)
  found_assigned <- FALSE
  if (any(grepl("assigned\\(.*\\)", cs))) {
    found_assigned <- TRUE
    assigned_form <- regmatches(cs, regexpr("assigned\\(.*\\)", cs))
    if (length(unique(assigned_form)) > 1) {
      stop("Differing forms of `assigned()` found. Keep form consistent.")
    }
    assigned_form <- assigned_form[1]
  }
  found_a. <- FALSE
  if (any(grepl("a\\.\\(.*\\)", cs))) {
    found_a. <- TRUE
    a._form <- regmatches(cs, regexpr("a\\.\\(.*\\)", cs))
    if (length(unique(a._form)) > 1) {
      stop("Differing forms of `a.()` found. Keep form consistent.")
    }
    a._form <- a._form[1]
  }
  found_z. <- FALSE
  if (any(grepl("z\\.\\(.*\\)", cs))) {
    found_z. <- TRUE
    z._form <- regmatches(cs, regexpr("z\\.\\(.*\\)", cs))
    if (length(unique(z._form)) > 1) {
      stop("Differing forms of `z.()` found. Keep form consistent.")
    }
    z._form <- z._form[1]

  }

  if (sum(found_assigned, found_a., found_z.) > 1) {
    stop(paste("Differing treatment identification (`assigned()`, `a.()`",
               " or `z.()`) found. Only one can be used."))
  }
  if (sum(found_assigned, found_a., found_z.) == 0) {
    stop("No treatment variable found")
  }

  if (found_assigned) {
    return(assigned_form)
  }
  if (found_a.) {
    return(a._form)
  }
  if (found_z.) {
    return(z._form)
  }
  stop("This error should never be hit!")
  return(NULL)


}
