#' @include Design.R WeightedDesign.R DesignAccessors.R
NULL
# The above ensures that `Design` and `WeightedDesign` are defined prior to
# `DirectAdjusted`

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
  if (!is_binary_or_dichotomized(object@Design)) {
    return("Treatment must be binary or have a dichotomy.")
  }
  if (all(object$model[1, 2] == object$model[, 2])) {
    return("treatment variable must not be constant")
  }
  if (!treatment(object) %in% names(coef(object))) {
    return("treatment not found in model")
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
##' @param design Optional, explicitly specify the \code{Design} to be used. If
##'   the \code{Design} is specified elsewhere in \code{x} (e.g. passed as an
##'   argument to any of \code{ate()}, \code{ett()}, \code{cov_adj()} or
##'   \code{adopters()}) it will be found automatically and does not need to be
##'   passed here as well. (If the \code{Design} is found in the model, this
##'   argument is ignored.)
##' @param target Optional, explicitly specify the estimand. If the model in
##'   \code{x} does not contain a \code{weights} argument generated using either
##'   \code{ate()} or \code{ett()}, specify whether the goal is estimating ATE
##'   ("ate") or ETT ("ett").
##' @param ... Add'l arguments
##' @return DirectAdjusted
##' @export
as.DirectAdjusted <- function(x, design = NULL, target = NULL, ...) {
  if (!is(x, "lm")) {
    stop("input must be lm object")
  }

  hasweights <- is(x$model$"(weights)", "WeightedDesign")
  hascov_adj <- is(x$model$"(offset)", "CovAdjPrediction")

  # Check if we can find a design in either Weights (preferred) or cov_adj
  found_design <- NULL
  if (hasweights) {
    found_design <- x$model$"(weights)"@Design
  } else if (hascov_adj) {
    found_design <- x$model$"(offset)"@Design
  }

  if (is(found_design, "Design")) {
    design <- found_design
  }
  if (!is(design, "Design")) {
    stop("Cannot locate `Design`, pass via `design=` argument")
  }

  # Check if we can find "target" in the weights
  found_target <- NULL
  if (hasweights) {
    found_target <- x$model$"(weights)"@target
  }

  if (found_target %in% c("ate", "ett")) {
    target <- found_target
  }

  if (!(target %in% c("ate", "ett"))) {
    stop('Cannot locate `target`, pass via `target=` argument ("ate" or "ette")')
  }

  new("DirectAdjusted",
      x,
      Design = design,
      target = target)
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
setMethod("confint", "DirectAdjusted",
          function(object, parm, level = 0.95, ...) {
  confint(as(object, "lm"), parm, level = level, ...)
})

##' Identify treatment variable in DirectAdjusted object
##' @param x DirectAdjusted model
##' @param ... Ignored
##' @return Name of treatment in model.
##' @author Josh Errickson
##' @importFrom stats coef
##' @examples
##' data(simdata)
##' des <- rct_design(z ~ unitid(cid1, cid2), data = simdata)
##' mod <- lm(y ~ z, data = simdata, weights = ett(des))
##' damod <- as.DirectAdjusted(mod)
##' damod$coef[treatment(damod)]
##' des2 <- rct_design(dose ~ unitid(cid1, cid2), data = simdata,
##'                    dichotomy = dose > 200 ~ . )
##' mod2 <- lm(y ~ adopters(), data = simdata, weights = ett(des2))
##' damod2 <- as.DirectAdjusted(mod2)
##' damod2$coef[treatment(damod2)]
setMethod("treatment", "DirectAdjusted", function(x, ...) {
  if (has_binary_treatment(x@Design)) {
    # Only if Design has a truly binary treatment variable...
    zname <- var_names(x@Design, "t")
    if (zname %in% names(stats::coef(x))) {
      # If treatment variable name is found in coefficients, return it
      return(zname)
    }
    # Otherwise return adopters
    return("adopters()")
  }
  # If Design's treatment variable isn't exactly binary, always require
  # adopters().
  return("adopters()")
})
