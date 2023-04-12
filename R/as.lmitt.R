#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R
#' @include DirectAdjusted.R
NULL

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
##' @importFrom stats formula
##' @export
as.lmitt <- function(x, design = NULL) {
  if (!inherits(x, "lm")) {
    stop("input must be lm object")
  }


  ### 10/31/22 JE - The below causes a ton of errors in the test suite that do
  ### NOT show up interactively and thus are extremely challenging to debug. It
  ### is supposed to replace the treatment variable name with `assigned()` and
  ### refit the model. I've spent a few days trying to debug this and failing,
  ### so going to leave it as unfixable for the moment, and instead simply error
  ### if the user does not include `assigned()` in the `lm`.
  ## # Update formula to use `assigned()` if needed
  ## ff <- stats::formula(x)
  ## newff <- formula(gsub(var_names(design, "t"), "assigned()", deparse(ff)))
  ## environment(newff) <- environment(ff)

  ## # Ensure updated model will be the same.

  ## newx <- update(x, newff)
  ## if (!isTRUE(all.equal(newx$coefficients, x$coefficients,
  ##                       check.attributes = FALSE))) {
  ##   stop(paste("Treatment variable found in model formula.",
  ##              "Updating model to use `assigned()` instead produces",
  ##              "different results. Please refit original model using",
  ##              "`assigned()` in place of treatment variable name."))
  ## }
  ## x <- newx

  tt <- terms(stats::formula(x), specials = c("assigned", "a.", "z."))

  if (all(vapply(attr(tt, "specials"), is.null, logical(1)))) {
    stop(paste("`assigned()` or its aliases are not found in the model formula.",
               "`assigned()` needs to be found in place of the treatment",
               "variable name. The `lmitt()` function may be used to",
               "avoid explicitly indicating `assigned()`."))
  }


  return(.convert_to_lmitt(x,
                           design,
                           lmitt_fitted = FALSE,
                           absorbed_intercepts = FALSE,
                           absorbed_moderators = vector("character")))
}

##' @rdname as_lmitt
##' @export
as.DirectAdjusted <- as.lmitt

.convert_to_lmitt <- function(lm_model,
                              design,
                              lmitt_fitted,
                              absorbed_intercepts,
                              absorbed_moderators) {
  if (!inherits(lm_model, "lm")) {
    stop("input must be lm object")
  }

  # Ensure `design=` is a proper object
  if (!is.null(design) & !is(design, "Design")) {
    # Allow WeightedDesign just in case
    if (is(design, "WeightedDesign")) {
      design <- design@Design
    } else {
      stop(paste("If provided, `design` must be a `Design` or",
                 "`WeightedDesign` object"))
    }
  }

  # Check if we can find a design in either Weights (preferred) or cov_adj
  design_weights <- tryCatch(lm_model$model$"(weights)"@Design,
                             error = function(e) NULL)
  design_cov_adj <- tryCatch(.get_cov_adj(lm_model)@Design,
                             error = function(e) NULL)

  # The list contains all designs possible found (one passed in, and one in each
  # of weights and cov_adj). Passing `unique` removes any duplicates (since
  # duplicates are OK).
  unique_designs <- unique(list(design, design_weights, design_cov_adj))
  # Drop any NULL designs
  unique_designs <- unique_designs[!vapply(unique_designs,
                                           is.null, logical(1))]
  # At this point, if the length of `unique_designs` is 1, we're done. More than
  # one is an error.
  if (length(unique_designs) == 1) {
    design <- unique_designs[[1]]
  } else if (length(unique_designs) > 1) {
    stop("Multiple differing `Design` found in object.")
  } else {
    # This should never be hit
    stop("Cannot locate a `Design`, pass via it `design=` argument")
  }

  eval_env <- new.env(parent = environment(formula(lm_model)))
  # Find data
  if (lmitt_fitted) {
    # If `lmitt.formula` is called, get the data from there directly (since
    # inside `lmitt.formula`, we pass in the data directly after appenindg on
    # the updated RHS and LHS).
    data <- lm_model$call$data
  } else {
    # If `as.lmitt` (or `lmitt.lm`), evaluate the lm call's data
    data <- eval(lm_model$call$data, envir = eval_env)
  }
  lm_model$call$data <- data
  assign("data", data, envir = eval_env)
  assign("design", design, envir = eval_env)
  environment(lm_model$terms) <- eval_env
  if (inherits(lm_model, "glm")) {
    lm_model$formula <- as.formula(lm_model, env = eval_env)
  }


  return(new("DirectAdjusted",
             lm_model,
             Design = design,
             absorbed_intercepts = absorbed_intercepts,
             absorbed_moderators = absorbed_moderators))

}
