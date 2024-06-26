#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R
#' @include teeMod.R
NULL

##' @title Convert \code{lm} object into \code{teeMod}
##'
##' @description Converts the output of [lm()] into a \code{teeMod}
##'   object so that it can be used to properly account for standard errors.
##'
##' @details The formula with which \code{x} was created must include a
##'   treatment identifier (e.g. [assigned()]).
##'
##' @param x \code{lm} object with weights containing a \code{WeightedDesign},
##'   or an offset from [cov_adj()].
##' @param design Optional, explicitly specify the \code{Design} to be used. If
##'   the \code{Design} is specified elsewhere in \code{x} (e.g. passed as an
##'   argument to any of \code{ate()}, \code{ett()}, \code{cov_adj()} or
##'   \code{assigned()}) it will be found automatically and does not need to be
##'   passed here as well. (If different \code{Design} objects are passed
##'   (either through the \code{lm} in weights or covariance adjustment, or
##'   through this argument), an error will be produced.)
##' @return \code{teeMod} object
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

  tt <- terms(stats::formula(x),
              specials = c("assigned", "a.", "z.", "adopters"))

  if (all(vapply(attr(tt, "specials"), is.null, logical(1)))) {
    stop(paste("`assigned()` or its aliases are not found in the model formula.",
               "`assigned()` needs to be found in place of the treatment",
               "variable name. The `lmitt()` function may be used to",
               "avoid explicitly indicating `assigned()`."))
  }

  #### Obtain the proper call.
  # The first conditional is the user calls `lmitt()` and passes an `lm` object.
  # The second conditional if the user calls `lmitt.lm()` directly
  # The third conditional is if the user calls `as.lmitt()`.

  # `&&` necessary to return FALSE immediately if not enough frames on stack
  if (sys.nframe() >=3 &&
      !is.null(sys.call(-2)) &&
      sys.call(-2)[[1]] == as.name("lmitt")) {
    # If we're in `lmitt.lm()` via `lmitt()`, save the call to `lmitt`.
    lmitt_call <- sys.call(-2)
  } else if (sys.nframe() >= 2 &&
             is.name(sys.call(-1)[[1]]) && # Catches mapply issue
             !is.null(sys.call(-1)) &&
             sys.call(-1)[[1]] == as.name("lmitt.lm")) {
    # If we're in `lmitt.lm()` directly, save that call.
    lmitt_call <- sys.call(-1)
  } else {
    # Otherwise save the `as.lmitt()` call
    lmitt_call <- sys.call()
  }
  return(.convert_to_lmitt(x,
                           design,
                           lmitt_fitted = FALSE,
                           absorbed_intercepts = FALSE,
                           moderator = vector("character"),
                           lmitt_call = lmitt_call))
}

##' @rdname as_lmitt
##' @export
as.teeMod <- as.lmitt

.convert_to_lmitt <- function(lm_model,
                              design,
                              lmitt_fitted,
                              absorbed_intercepts,
                              moderator,
                              lmitt_call) {
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

  designs <- list(design, design_weights, design_cov_adj)
  designs <- designs[!vapply(designs, is.null, logical(1))]

  # The list contains all designs possible found (one passed in, and one in each
  # of weights and cov_adj). Passing `unique` removes any duplicates (since
  # duplicates are OK).
  unique_designs <- unique(designs)

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

  # #123 Make PreSandwichLayer to SandwichLayer if necessary
  ca <- .get_cov_adj(lm_model)
  if (is(ca, "PreSandwichLayer") & !is(ca, "SandwichLayer")) {
    lm_model$model$`(offset)` <- as.SandwichLayer(.get_cov_adj(lm_model),
                                                  design = design)
  }

  eval_env <- new.env(parent = environment(formula(lm_model)))
  # Find data
  if (lmitt_fitted) {
    # If `lmitt.formula` is called, get the data from there directly (since
    # inside `lmitt.formula`, we pass in the data directly after appenindg on
    # the updated RHS and LHS).
    data <- lm_model$call$data
  } else {
    # If `as.lmitt` (or `lmitt.lm`), evaluate the lm call's data using a
    # `fallback_data_search`-esque approach
    data_call <- lm_model$call$data
    #for (i in c(-1L, seq_len(sys.nframe()))) {
    for (i in seq_len(length(sys.calls()))) {
      try(data <- eval(data_call, envir = parent.frame(i)),
          silent = TRUE)
      if (!is.null(data) && is.data.frame(data)) {
        break()
      }
    }
    if (!inherits(data, "data.frame")) {
      stop("Could not determine appropriate data")
    }
    if (!is.null(lm_model$call$subset)) {
      sbst <- eval(lm_model$call$subset, data)
      data <- subset(data, sbst)
    }
  }
  lm_model$call$data <- data
  assign("data", data, envir = eval_env)
  assign("design", design, envir = eval_env)
  environment(lm_model$terms) <- eval_env
  if (inherits(lm_model, "glm")) {
    lm_model$formula <- as.formula(lm_model, env = eval_env)
  }

  return(new("teeMod",
             lm_model,
             Design = design,
             absorbed_intercepts = absorbed_intercepts,
             moderator = moderator,
             lmitt_call = call("quote", lmitt_call),
             lmitt_fitted = lmitt_fitted))

}
