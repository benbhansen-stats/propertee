#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R
NULL
# The above ensures that `Design`, `WeightedDesign`, and `vcovDA` are defined
# prior to `DirectAdjusted`

setClass("DirectAdjusted",
         contains = "lm",
         slots = c(Design = "Design"))

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
  call <- match.call()

  if (is.null(call[["type"]])) {
    confint_calls <- grepl("confint.DirectAdjusted", lapply(sys.calls(), "[[", 1))
    if (any(confint_calls)) {
      type <- tryCatch(get("call", sys.frame(which(confint_calls)[1]))$type,
                       error = function(e) NULL)
      call$type <- type # will not append if type is NULL
    }
  }

  call[[1L]] <- if (inherits(object$model$`(offset)`, "SandwichLayer")) vcovDA else getS3method("vcov", "lm")
  vmat <- eval(call, parent.frame())

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
