##' Linear Model for Intention To Treat
##'
##' Generates a linear model object which allows proper estimation of variances
##' accounting for the study design.
##'
##' The first argument to \code{lmitt}, \code{obj}, specifies a linear
##' regression model. The formula is specified exactly how it would be in
##' \code{lm()}, with one deviation. Rather than including the treatment
##' variable by name, the helper function \code{assigned()} should be used in
##' its place. For example, if your outcome is \code{y}, you could call
##' \code{lmitt(y ~ assigned(), design = ...)}.
##'
##' If the formula does \emph{not} contain at least one \code{assigned()}, the
##' right hand side of the formula will be treated as a series of
##' strataification variables. This means that a formula such as \code{y ~ x +
##' q} will be replaced with \code{y ~ x:assigned() + q:assigned()} such that
##' \code{lmitt} will estimate the treatment effect within each level of both
##' \code{x} and \code{q}. If \code{assigned()} is found anywhere on the right
##' hand side of \code{obj}, this modification is \emph{not} made.
##'
##' Any additional arguments to \code{lm()} can be passed into \code{...}. See
##' the help for \code{lm()} for a list of possible argument. Most commonly used
##' in this scenario would be \code{weights=} to pass \code{ate()} or
##' \code{ett()}, or \code{offset=} to pass \code{cov_adj()}.
##'
##' Alternatively, \code{obj} can be a pre-created \code{lm} object. No
##' modification is made to the formula of the object. See the help for
##' \code{as.lmitt()} for details of this conversion.
##'
##' Note that although the \code{Design} creation functions (e.g.
##' \code{rct_design()}) take an optional \code{subset=} argument used in the
##' creation of the \code{Design}, this is \bold{not} the same as the
##' \code{subset=} argument passed to \code{lm()} or \code{lmitt()}. The
##' \code{subset=} argument when creating a \code{Design} restricts the data
##' used to generate the \code{Design}, but has no direct impact on the future
##' \code{lm()} or \code{lmitt()} calls using that \code{Design}. (It can
##' indirectly have an impact by excluding particular clusters/units of
##' assignment from recieving a treatment assignment and thus complete case
##' analysis removes them from the model.)
##'
##' On the other hand, the \code{subset=} argument in \code{lm()} or
##' \code{lmitt()} refers only to subsetting the \code{data} argument passed
##' into \code{lm()} or \code{lmitt()}.
##' @param obj A \code{formula} or a \code{lm} object. See details.
##' @param design Optional, explicitly specify the \code{Design} to be used. If
##'   the \code{Design} is specified elsewhere in the model (e.g. passed as an
##'   argument to any of \code{ate()}, \code{ett()}, \code{cov_adj()} or
##'   \code{assigned()}) it will be found automatically and does not need to be
##'   passed here as well. (If different \code{Design} objects are passed
##'   (either through the \code{lm} in weights or covariance adjustment, or
##'   through this argument), an error will be produced.) Alternatively, a
##'   formula creating a design (of the type of that would be passed as the
##'   first argument to \code{rd_design()}, \code{rct_design()}, or
##'   \code{obs_design()}.
##' @param absorb If \code{TRUE}, fixed effects are included for units of
##'   assignemnt/clusters identified in the \code{Design}. Excluded in
##'   \code{FALSE}. Default is \code{FALSE}.
##' @param ... Additional arguments passed to \code{lm()}. Ignored if \code{obj}
##'   is already an \code{lm} object.
##' @return \code{DirectAdjusted} model.
##' @export
##' @importFrom stats lm predict weights
##' @rdname lmitt
lmitt <- function(obj,
                  design = NULL,
                  absorb = FALSE,
                  ...) {
  UseMethod("lmitt")
}

##' @export
##' @rdname lmitt
lmitt.formula <- function(obj,
                          design = NULL,
                          absorb = FALSE,
                          ...) {
  mf <- match.call()

  if (!is.null(attr(terms(obj, specials = ".absorbed"),
                    "specials")$.absorbed)) {
    stop(paste("`.absorbed()` is an internal function",
               "and should not be used by end-users"))
  }


  # If there are no assigned() in the formula, assume all RHS variables are
  # stratified and add interaction with `assigned()`
  no_assigned <- is.null(attr(terms(obj, specials = "assigned"),
                              "specials")$assigned)
  if (no_assigned) {
      obj <- update(obj, . ~ . : assigned())
  }

  if (absorb) {
    fixed_eff_term <- paste(paste0(".absorbed(", var_names(design, "b"), ")"),
                          collapse = "*")
    obj <- update(obj, paste0(". ~ . + ", fixed_eff_term))
  }

  if (inherits(design, "formula")) {
    # If there's a `forcing()`, user wants RDD. If not, force Obs. To do RCT,
    # must create Design manually.
    if (!is.null(attr(terms(design, specials = "forcing"),
                      "specials")$forcing)) {
      des_call <- "rd_design"
    } else {
      des_call <- "obs_design"
    }

    # Build new call. All calls must include obj and data
    new_d_call <- paste0(des_call, "(",
                         "formula = ", deparse(design),
                         ", data = ", deparse(mf$data))
    # If user passed dichotomy, include it. We do this so the
    # `design@call` will be in agreement.
    if (!is.null(mf$dichotomy)) {
      new_d_call <- paste0(new_d_call, ", dichotomy = ", deparse(mf$dichotomy))
    }
    new_d_call <- paste0(new_d_call, ")")
    # str2lang converts character into call
    design <- eval(str2lang(new_d_call))
  }

  m <- match(c("obj", "data", "subset", "weights", "na.action",
               "method", "model", "x", "y", "qr", "singular.ok",
               "contrasts", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::lm)

  # Reset the formula  for the `lm`, giving  it a proper name (since  we take in
  # the  generic "obj"  name), and  replacing it  with the  updated `obj`  if we
  # modified it above due to lack of `assigned()
  names(mf)[2] <- "formula"
  mf[[2L]] <- obj

  model <- eval(mf, parent.frame())

  model$call[[1]] <- as.name("lmitt")

  return(as.lmitt(model, design))

}

##' @export
##' @rdname lmitt
lmitt.lm <- function(obj,
                     design = NULL,
                     ...) {
  return(as.lmitt(obj, design))
}

##' For internal use only
##' @export
##' @param x x
.absorbed <- as.factor
