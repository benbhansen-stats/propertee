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
##'   is already an \code{lm} object. If \code{weights=} is passed, it can be
##'   either a call to \code{ate()} or \code{ett()}, or if no additional
##'   arguments or manipulations of those weights are needed, the strings
##'   \code{"ate"} or \code{"ett"}.
##' @return \code{DirectAdjusted} model.
##' @export
##' @importFrom stats lm predict weights
##' @rdname lmitt
lmitt <- function(obj,
                  design,
                  absorb = FALSE,
                  ...) {
  UseMethod("lmitt")
}

##' @export
##' @rdname lmitt
lmitt.formula <- function(obj,
                          design,
                          absorb = FALSE,
                          ...) {
  mf <- match.call()

  if (!is.null(attr(terms(obj, specials = ".absorbed"),
                    "specials")$.absorbed)) {
    stop(paste("`.absorbed()` is an internal function",
               "and should not be used by end-users"))
  }

  # First, make sure we have a valid `design=` - if given a formula, make a new
  # `Design`, otherwise ensure `design=` is `Design` class.
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
  } else if (!is(design, "Design")) {
    stop(paste("`design=` must be an object created by `*_design`",
               "function, or a formula specifying such a design"))
  }

  ### Next, update main formula

  # Try to ensure proper formula. Valid forms are:
  # ~ 1
  # ~ sbgrp
  # ~ 1 + sbgrp <- NYI, same as ~ sbgrp right now
  rhs <- trimws(strsplit(deparse(obj[[3]]), "+", fixed = TRUE)[[1]])
  rhs <- rhs[rhs != "+"]
  if (length(rhs) > 2) {
    stop("Too many terms on right hand side")
  }
  if ("0" %in% rhs) {
    stop("'0' is not a valid entry on right hand side")
  }
  if (length(rhs[rhs != "1"]) > 1) {
    stop("Too many variable on right hand ide")
  }
  if (var_names(design, "t") %in% rhs) {
    stop("Treatment variable should not be manually entered.")
  }

  #### This blocks on `~ 1 + sbgrp` and `~ sbgrp + 1`. When we eventually
  #### enable this functionality, remove this block
  #### See #73
  if (length(rhs) == 2) {
    stop("To estimate subgroup effects, do not include '1'")
  }

  # About the modify the formula; doing so strips off envir, so saving to
  # restore below.
  saveenv <- environment(obj)

  # Replace alises for assigned with assigned
  obj <- formula(gsub("adopters(", "assigned(", deparse(obj), fixed = TRUE))
  obj <- formula(gsub("a.(",       "assigned(", deparse(obj), fixed = TRUE))
  obj <- formula(gsub("z.(",       "assigned(", deparse(obj), fixed = TRUE))
  # Only using opening parans to catch things like `a.(des)`.
  environment(obj) <- saveenv


  # If there are no assigned() in the formula, assume all RHS variables are
  # stratified and add interaction with `assigned()`


  no_assigned <- is.null(attr(terms(obj, specials = "assigned"),
                              "specials")$assigned)
  if (no_assigned) {
    if (length(attr(terms(obj), "term.labels")) == 0) {
      # If user passes `y ~ 1`, there will be no `terms`, and `1:assigned`
      # is ignored, so add it instead
      obj <- update(obj, . ~ . + assigned())
    } else {
      # There's a sbgrp (as checked above) so interact it with assigned
      obj <- update(obj, . ~ . : assigned())
    }
  }

  # Handle the `absorb=` argument
  if (absorb) {
    fixed_eff_term <- paste(paste0(".absorbed(", var_names(design, "b"), ")"),
                          collapse = "*")
    obj <- update(obj, paste0(". ~ . + ", fixed_eff_term))
  }


  m <- match(c("obj", "data", "subset", "weights", "na.action",
               "method", "model", "x", "y", "qr", "singular.ok",
               "contrasts", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::lm)

  ### Allow users to pass in "ate" and "ett" rather than functions if they have
  ### no special modifications/additional arguments
  wt <- substitute(list(...))$weights
  if (is(wt, "character")) {
    if (tolower(wt) == "ate") {
      mf$weights <- quote(ate())
    } else if (tolower(wt) == "ett") {
      mf$weights <- quote(ett())
    } else {
      warning(paste("Character other than \"ate\" or \"ett\" passed to",
                    "`weights=` argument.\nIf you are trying to pass a",
                    "character to the internal `lm` you can disregard this",
                    "warning.\nIf you are attemping to use `flexida`\'s",
                    "weight generation, only \"ate\" and \"ett\" are",
                    "accepted."))
    }
  }

  # Reset the formula  for the `lm`, giving  it a proper name (since  we take in
  # the  generic "obj"  name), and  replacing it  with the  updated `obj`  if we
  # modified it above due to lack of `assigned()
  names(mf)[2] <- "formula"
  mf[[2L]] <- obj

  model <- eval(mf, parent.frame())

  return(as.lmitt(model, design))

}

##' @export
##' @rdname lmitt
lmitt.lm <- function(obj,
                     design = NULL,
                     ...) {
  return(as.lmitt(obj, design))
}

##' (Internal) Include blocks as absorbed effects
##'
##' Not to be used interactively. Alias for \code{as.factor()} to identify fixed
##' effects generated through the \code{absorb=TRUE} option to \code{lmitt()}.
##' @param x Block variables
##' @return Identifies block variables as categorical for the model.
##' @export
.absorbed <- as.factor
