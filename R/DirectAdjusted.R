DirectAdjusted <- setClass("DirectAdjusted",
                           contains = "lm",
                           slots = c(Design = "Design",
                                     target = "character"))

setValidity("DirectAdjusted", function(object) {
  if (length(object@target) != 1 || !object@target %in% c("ett", "ate")) {
    return(paste("@target must be one of [ett, ate]. unknown @target:",
                 paste(object@target, collapse = " ")))
  }
  if (!is(object$model$"(weights)", "WeightedDesign")) {
    return("weight must be WeightedDesign created by `ate()` or `ett()`")
  }
  if (ncol(object$model) < 3) {
    return("model must contain a binary treatment predictor")
  }
  if (!all(object$model[,2] %in% 0:1)) {
    return("first predictor must be binary treatment indicator")
  }
  if (sd(object$model[,2]) == 0) {
    return("treatment indicator is constant")
  }
  if (is.na(object$coef[2])) {
    return("treatment effect failed to estimate")
  }
  TRUE
})

##' @title Show an DirectAdjusted
##' @param object DirectAdjusted object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "DirectAdjusted", function(object) {
  print(as(object, "lm"))
  invisible(object)
})

##' @title Convert lm object into DirectAdjusted
##' @param x lm object with weights containing a WeightedDesign
##' @param ... Add'l arguments
##' @return DirectAdjusted
##' @export
as.DirectAdjusted <- function(x, ...) {
  if (!is(x, "lm")) {
    stop("input must be lm object")
  }

  if (!is(x$model$"(weights)", "WeightedDesign")) {
    stop("input model must contain WeightedDesign weights")
  }

  DirectAdjusted(x,
                 Design = x$model$"(weights)"@Design,
                 target = x$model$"(weights)"@target)
}

setGeneric("vcov")

##' @title Variance-Covariance matrix
##' @param object DirectAdjusted
##' @param ... Add'l arguments
##' @return Variance-Covariance matrix
##' @export
setMethod("vcov", "DirectAdjusted", function(object, ...) {
  vcov(as(object, "lm"), ...)
})


setGeneric("confint")

##' @title Variance-Covariance matrix
##' @param object DirectAdjusted
##' @param parm a specification of which parameters are to be given confidence
##'   intervals, either a vector of numbers or a vector of names. If missing,
##'   all parameters are considered.
##' @param level the confidence level required.
##' @param ... Add'l arguments
##' @return Variance-Covariance matrix
##' @export
setMethod("confint", "DirectAdjusted", function(object, parm, level = 0.95, ...) {
  confint(as(object, "lm"), parm, level = level, ...)
})
