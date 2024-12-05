#' @include StudySpecification.R WeightedStudySpecification.R StudySpecificationAccessors.R SandwichLayerVariance.R
#' @include teeMod.R
NULL

##' @title Convert \code{lm} object into \code{teeMod}
##'
##' @description Converts the output of [lm()] into a \code{teeMod}
##'   object, for standard errors that account for block and cluster
##'   information carried with the \code{lm}'s weights, and/or an
##'   offset incorporating predictions of the outcome from a
##'   separate model.
##'
##' @details The formula with which \code{x} was created must include
##'   a treatment identifier (e.g. [assigned()]).  If a model-based
##'   offset is incorportated, the model's predictions would have to
##'   have been extracted using [cov_adj()] (as opposed to
##'   \code{predict{}} in order for \code{teeMod} standard error
##'   calculations to reflect propagation of error from these
##'   predictions. This mechanism only supports treatment main effects:
##'   to estimate interactions of treatment assignment with a moderator
##'   variable, use [lmitt()] instead of \code{lm()} and
##'   \code{as.lmitt()}.
##'
##' @param x \code{lm} object with weights containing a
##'   \code{WeightedStudySpecification}, or an offset from [cov_adj()].
##' @param specification Optional, explicitly specify the
##'   \code{StudySpecification} to be used. If the \code{StudySpecification} is
##'   specified elsewhere in \code{x} (e.g. passed as an argument to any of
##'   \code{ate()}, \code{ett()}, \code{cov_adj()} or \code{assigned()}) it will
##'   be found automatically and does not need to be passed here as well. (If
##'   different \code{StudySpecification} objects are passed (either through the
##'   \code{lm} in weights or covariance adjustment, or through this argument),
##'   an error will be produced.)
##' @return \code{teeMod} object
##' @rdname as_lmitt
##' @importFrom stats formula
##' @export
as.lmitt <- function(x, specification = NULL) {
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
  ## newff <- formula(gsub(var_names(specification, "t"), "assigned()", deparse(ff)))
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
                           specification,
                           lmitt_fitted = FALSE,
                           absorbed_intercepts = FALSE,
                           moderator = vector("character"),
                           lmitt_call = lmitt_call))
}

##' @rdname as_lmitt
##' @export
as.teeMod <- as.lmitt

.convert_to_lmitt <- function(lm_model,
                              specification,
                              lmitt_fitted,
                              absorbed_intercepts,
                              moderator,
                              lmitt_call) {
  if (!inherits(lm_model, "lm")) {
    stop("input must be lm object")
  }

  # Ensure `specification=` is a proper object
  if (!is.null(specification) & !is(specification, "StudySpecification")) {
    # Allow WeightedStudySpecification just in case
    if (is(specification, "WeightedStudySpecification")) {
      specification <- specification@StudySpecification
    } else {
      stop(paste("If provided, `specification` must be a `StudySpecification` or",
                 "`WeightedStudySpecification` object"))
    }
  }

  # Check if we can find a specification in either Weights (preferred) or cov_adj
  specification_weights <- tryCatch(lm_model$model$"(weights)"@StudySpecification,
                             error = function(e) NULL)
  specification_cov_adj <- tryCatch(.get_cov_adj(lm_model)@StudySpecification,
                             error = function(e) NULL)

  specifications <- list(specification, specification_weights, specification_cov_adj)
  specifications <- specifications[!vapply(specifications, is.null, logical(1))]

  # The list contains all specifications possible found (one passed in, and one in each
  # of weights and cov_adj). Passing `unique` removes any duplicates (since
  # duplicates are OK).
  unique_specs <- unique(specifications)

  # At this point, if the length of `unique_specs` is 1, we're done. More than
  # one is an error.
  if (length(unique_specs) == 1) {
    specification <- unique_specs[[1]]
  } else if (length(unique_specs) > 1) {
    stop("Multiple differing `StudySpecification` found in object.")
  } else {
    # This should never be hit
    stop("Cannot locate a `StudySpecification`, pass via it `specification=` argument")
  }

  # #123 Make PreSandwichLayer to SandwichLayer if necessary
  ca <- .get_cov_adj(lm_model)
  if (is(ca, "PreSandwichLayer") & !is(ca, "SandwichLayer")) {
    lm_model$model$`(offset)` <- as.SandwichLayer(.get_cov_adj(lm_model),
                                                  specification = specification)
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
  # set call's na.action to na.pass so expand.model.frame includes NA rows
  lm_model$call$na.action <- "na.pass"
  assign("data", data, envir = eval_env)
  assign("specification", specification, envir = eval_env)
  environment(lm_model$terms) <- eval_env
  if (inherits(lm_model, "glm")) {
    lm_model$formula <- as.formula(lm_model, env = eval_env)
  }

  return(new("teeMod",
             lm_model,
             StudySpecification = specification,
             absorbed_intercepts = absorbed_intercepts,
             moderator = moderator,
             lmitt_call = call("quote", lmitt_call),
             lmitt_fitted = lmitt_fitted))

}
