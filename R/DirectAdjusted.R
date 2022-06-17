#' @include Design.R WeightedDesign.R DesignAccessors.R
NULL
# The above ensures that `Design` and `WeightedDesign` are defined prior to
# `DirectAdjusted`

setClass("DirectAdjusted",
         contains = "lm",
         slots = c(Design = "Design",
                   target = "character"))

setValidity("DirectAdjusted", function(object) {
  if (length(object@target) != 1 || !object@target %in% c("ett", "ate")) {
    return(paste("@target must be one of [ett, ate]. unknown @target:",
                 paste(object@target, collapse = " ")))
  }
  if (!is_binary_or_dichotomized(object@Design)) {
    return("Treatment must be binary or have a dichotomy.")
  }
  if (!treatment_name(object) %in%
        rownames(attr(terms(object), "factors"))[-1]) {
    # Ensures that treatment variable appears somewhere in RHS (-1 removes
    # outcome) of formula. If created by `lmitt()`, treatment will be
    # "adopters()"; but if passed from lm to as.DA, it could be a variable name.
    return("treatment not found in model")
  }
  return(TRUE)
})

##' @title Show an DirectAdjusted
##' @param object DirectAdjusted object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "DirectAdjusted", function(object) {
  print(as(object, "lm"))
  invisible(object)
})

##' @title Convert \code{lm} object into \code{DirectAdjusted}
##' @param x \code{lm} object with weights containing a \code{WeightedDesign}
##' @param design Optional, explicitly specify the \code{Design} to be used. If
##'   the \code{Design} is specified elsewhere in \code{x} (e.g. passed as an
##'   argument to any of \code{ate()}, \code{ett()}, \code{cov_adj()} or
##'   \code{adopters()}) it will be found automatically and does not need to be
##'   passed here as well. (If the \code{Design} is found in the model, this
##'   argument is ignored.)
##' @param target Optional, explicitly specify the estimand. If the model in
##'   \code{x} does not contain a \code{weights} argument generated using either
##'   \code{ate()} or \code{ett()}, specify whether the goal is estimating ATE
##'   ("ate") or ETT ("ett"). (If weights are specified, this argument is
##'   ignored.)
##' @return \code{DirectAdjusted} object
##' @export
as.DirectAdjusted <- function(x, design = NULL, target = NULL) {
  if (!is(x, "lm")) {
    stop("input must be lm object")
  }

  hasweights <- is(x$model$"(weights)", "WeightedDesign")

  # Check if we can find a design in either Weights (preferred) or cov_adj
  found_design <- NULL
  if (hasweights) {
    found_design <- x$model$"(weights)"@Design
  } else {
    capred <- .get_cov_adj(x)
    if (is(capred, "CovAdjPrediction")) {
      found_design <- capred@Design
    }
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

  if (!is.null(found_target) && found_target %in% c("ate", "ett")) {
    target <- found_target
  }

  if (is.null(target) || !(target %in% c("ate", "ett"))) {
    stop(paste('Cannot locate `target`, pass via `target=` argument',
               '("ate" or "ett")'))
  }

  return(new("DirectAdjusted",
             x,
             Design = design,
             target = target))
}

setGeneric("vcov")

##' @title Variance-Covariance matrix
##' @param object DirectAdjusted
##' @param ... Add'l arguments
##' @return Variance-Covariance matrix
##' @export
setMethod("vcov", "DirectAdjusted", function(object, ...) {
  return(vcov(as(object, "lm"), ...))
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
  return(confint(as(object, "lm"), parm, level = level, ...))
})

##' Identify treatment variable in \code{DirectAdjusted} object
##'
##' @param x \code{DirectAdjusted} model
##' @return Name of treatment in model.
##' @export
##' @examples
##' data(simdata)
##' des <- rct_design(z ~ unitid(cid1, cid2), data = simdata)
##' mod <- lm(y ~ z, data = simdata, weights = ett(des))
##' damod <- as.DirectAdjusted(mod)
##' damod$coef[treatment_name(damod)]
##' des2 <- rct_design(dose ~ unitid(cid1, cid2), data = simdata,
##'                    dichotomy = dose > 200 ~ . )
##' mod2 <- lm(y ~ adopters(), data = simdata, weights = ett(des2))
##' damod2 <- as.DirectAdjusted(mod2)
##' damod2$coef[treatment_name(damod2)]
treatment_name <- function(x) {

  cnames <- names(x$coefficients)
  adopters_regexp <- "adopters\\([^)]*\\)"
  adopters_match <- regmatches(cnames, regexpr(adopters_regexp, cnames))
  if (length(unique(adopters_match)) > 1) {
    stop(paste("Differing `adopters()` calls found;",
               " all `adopters()` in formula must be identical."))
  }
  if (length(adopters_match) > 0) {
    return(adopters_match[1])
  }

  if (has_binary_treatment(x@Design)) {
    # Only if Design has a truly binary treatment variable...
    zname <- var_names(x@Design, "t")
    if (zname %in% cnames) {
      # If treatment variable name is found in coefficients, return it
      return(zname)
    }
    stop(paste("Treatment", zname, "or `adopters()` must be found in formula"))
  }
  # If we hit this point, there's no adopter and non-binary treatment, so we
  # must have non-binary treatment specified
  stop("With non-binary treatment, `adopters()` must be found in formula")
}
