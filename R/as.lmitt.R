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
