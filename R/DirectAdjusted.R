#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R
NULL
# The above ensures that `Design`, `WeightedDesign`, and `vcovDA` are defined
# prior to `DirectAdjusted`

setClass("DirectAdjusted",
         contains = "lm",
         slots = c(Design = "Design"))

setValidity("DirectAdjusted", function(object) {
 if (!treatment_name(object) %in%
        rownames(attr(terms(object), "factors"))[-1]) {
    # Ensures that treatment variable appears somewhere in RHS (-1 removes
    # outcome) of formula. If created by `lmitt()`, treatment will be
    # "assigned()"; but if passed from lm to as.DA, it could be a variable name.
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
##'   \code{assigned()}) it will be found automatically and does not need to be
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

  eval_env <- new.env(parent = environment(formula(x)))
  data <- eval(x$call$data, eval_env)
  x$call$data <- quote(data)
  assign("data", data, envir = eval_env)
  assign("design", design, envir = eval_env)
  environment(x$terms) <- eval_env
  if (inherits(x, "glm")) {
    x$formula <- as.formula(x, env = eval_env)
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
##' @details If a \code{DirectAdjusted} object is fit with a \code{SandwichLayer}
##' offset, then the usual \code{stats::summary.lm()} output is enhanced by
##' the use of covariance-adjusted sandwich standard errors, with t-test values
##' recalculated to reflect the new standard errors.
##' @param object DirectAdjusted
##' @param ... Additional arguments to \code{vcovDA()}, such as the desired
##' finite sample heteroskedasticity-robust standard error adjustment.
##' @return summary object based from \code{stats:::summary.lm()}
##' @export
setMethod("summary", "DirectAdjusted", function(object, ...) {
  ans <- summary(as(object, "lm"))
  if (inherits(object$model$`(offset)`, "SandwichLayer")) {
    ans$coefficients[, 2L] <- sqrt(diag(vcovDA(object, ...)))
    ans$coefficients[, 3L] <- ans$coefficients[, 1L] / ans$coefficients[, 2L]
    ans$coefficients[, 4L] <- 2*stats::pt(abs(ans$coefficients[, 3L]), ans$df[2L],
                                          lower.tail = FALSE)
  }

  return(ans)
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
##' mod2 <- lm(y ~ assigned(), data = simdata, weights = ett(des2))
##' damod2 <- as.lmitt(mod2)
##' damod2$coef[treatment_name(damod2)]
treatment_name <- function(x) {

  cnames <- names(x$coefficients)
  assigned_regexp <- "assigned\\([^)]*\\)"
  assigned_match <- regmatches(cnames, regexpr(assigned_regexp, cnames))
  if (length(unique(assigned_match)) > 1) {
    stop(paste("Differing `assigned()` calls found;",
               " all `assigned()` in formula must be identical."))
  }
  if (length(assigned_match) > 0) {
    return(assigned_match[1])
  }

  if (has_binary_treatment(x@Design)) {
    # Only if Design has a truly binary treatment variable...
    zname <- var_names(x@Design, "t")
    if (zname %in% cnames) {
      # If treatment variable name is found in coefficients, return it
      return(zname)
    }
    stop(paste("Treatment", zname, "or `assigned()` must be found in formula"))
  }
  # If we hit this point, there's no adopter and non-binary treatment, so we
  # must have non-binary treatment specified
  stop("With non-binary treatment, `assigned()` must be found in formula")
}
