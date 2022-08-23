#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R
NULL
# The above ensures that `Design`, `WeightedDesign`, and `vcovDA` are defined
# prior to `DirectAdjusted`

setClass("DirectAdjusted",
         contains = "lm",
         slots = c(Design = "Design"))

setValidity("DirectAdjusted", function(object) {
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
##'   passed here as well. (If different \code{Design} objects are passed
##'   (either through the \code{lm} in weights or covariance adjustment, or
##'   through this argument), an error will be produced.)
##' @return \code{DirectAdjusted} object
##' @rdname as_lmitt
##' @export
as.lmitt <- function(x, design = NULL) {
  if (!inherits(x, "lm")) {
    stop("input must be lm object")
  }

  # Check if we can find a design in either Weights (preferred) or cov_adj
  design_weights <- tryCatch(x$model$"(weights)"@Design,
                             error = function(e) NULL)
  design_cov_adj <- tryCatch(.get_cov_adj(x)@Design,
                             error = function(e) NULL)

  # The list contains all designs possible found (one passed in, and one in each
  # of weights and cov_adj). Passing `unique` removes any duplicates (since
  # duplicates are OK).
  unique_designs <- unique(list(design, design_weights, design_cov_adj))
  # Drop any designs which aren' `Design`. Mostly NULL hopefully.
  unique_designs <- unique_designs[vapply(unique_designs,
                                          is, logical(1), "Design")]
  # At this point, if the lenght of `unique_designs` is 1, we're done. More than
  # one is an error.
  if (length(unique_designs) == 1) {
    design <- unique_designs[[1]]
  } else if (length(unique_designs) > 1) {
    stop("Multiple differing `Design` found in object.")
  } else {
    stop("Cannot locate a `Design`, pass via it `design=` argument")
  }

  return(new("DirectAdjusted",
             x,
             Design = design))
}

##' @rdname as_lmitt
##' @export
as.DirectAdjusted <- as.lmitt

setGeneric("summary")

##' @title Summary of \code{DirectAdjusted} object
##' @details The usual \code{stats::summary.lm()} output is enhanced by
##' covariance-adjusted sandwich standard errors, with t-test values
##' recalculated to reflect this adjustment.
##' @param object DirectAdjusted
##' @parm ... Add'l aguments
##' @return summary object based from \code{stats::summary.lm()}
##' @export
setMethod("summary", "DirectAdjusted", function(object, ...) {
  ans <- summary(as(object, "lm"))
  if (inherits(object$model$`(offset)`, "SandwichLayer")) {
    ans$coefficients[, 2L] <- sqrt(diag(vcovDA(object, ...)))
    ans$coefficients[, 3L] <- ans$coefficients[, 1L] / ans$coefficients[, 2L]
    ans$coefficients[, 4L] <- 2*pt(abs(ans$coefficients[, 3L]), ans$df[2L],
                                   lower.tail = FALSE)
  }

  return(ans)
})

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
##' damod <- as.lmitt(mod)
##' damod$coef[treatment_name(damod)]
##' des2 <- rct_design(dose ~ unitid(cid1, cid2), data = simdata,
##'                    dichotomy = dose > 200 ~ . )
##' mod2 <- lm(y ~ adopters(), data = simdata, weights = ett(des2))
##' damod2 <- as.lmitt(mod2)
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
