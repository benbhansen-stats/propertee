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
##' @importFrom stats lm predict weights weighted.mean reformulate
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
  lmitt.call <- match.call()
  m <- match(c("obj", "data", "subset", "weights", "na.action",
               "method", "model", "x", "y", "qr", "singular.ok",
               "contrasts", "offset"), names(lmitt.call), 0L)
  lm.call <- lmitt.call[c(1L, m)]
  names(lm.call)[2] <- "formula"
  lm.call[[1]] <- quote(stats::lm)
  #lmitt.call contains design, absorb. lm.call contains only things that get
  #passed to lm, model.matrix

  ### Allow users to pass in "ate" and "ett" rather than functions if they have
  ### no special modifications/additional arguments
  wt <- substitute(list(...))$weights
  if (is(wt, "character")) {
    if (tolower(wt) == "ate") {
      lm.call$weights <- quote(ate())
      lmitt.call$weights <- quote(ate())
    } else if (tolower(wt) == "ett") {
      lm.call$weights <- quote(ett())
      lmitt.call$weights <- quote(ett())
    } else {
      warning(paste("Character other than \"ate\" or \"ett\" passed to",
                    "`weights=` argument.\nIf you are trying to pass a",
                    "character to the internal `lm` you can disregard this",
                    "warning.\nIf you are attemping to use `flexida`\'s",
                    "weight generation, only \"ate\" and \"ett\" are",
                    "accepted."))
    }
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
                         ", data = ", deparse(lmitt.call$data))
    # If user passed dichotomy, include it. We do this so the
    # `design@call` will be in agreement.
    if (!is.null(lmitt.call$dichotomy)) {
      new_d_call <- paste0(new_d_call, ", dichotomy = ", deparse(lmitt.call$dichotomy))
    }
    new_d_call <- paste0(new_d_call, ")")
    # str2lang converts character into call
    design <- eval(str2lang(new_d_call))
  } else if (!is(design, "Design")) {
    stop(paste("`design=` must be an object created by `*_design`",
               "function, or a formula specifying such a design"))
  }

  # Extract formula bits
  rhs <- trimws(strsplit(deparse(obj[[3]]), "+", fixed = TRUE)[[1]])
  rhs <- rhs[rhs != "+"]
  if (length(rhs) > 1) {
    # At max two terms ("sbgrp" and "1")
    stop(paste("Too many terms on right hand side. To estimate subgroup",
               "effects and main effect simultaneously, fit two",
               "separate models."))
  }
  if (rhs == "0") {
    # We don't support y ~ 0
    stop("'0' is not a valid entry on right hand side")
  }
  if (grepl("assigned\\(", rhs)) {
    stop(paste("Do not specify `assigned()` in the right hand side of `lmitt()`.\n",
               "To estimate only a treatment effect, pass `~ 1` as the right",
               "hand side."))
  }
  # `rhs` is now either "1" or subgrouping variable


  # Evaluate model.frame
  mf.call <- lm.call
  mf.call[[1]] <- quote(stats::model.frame)
  # Add assigned() to the model so it utilizes Design characteristics (primarly
  # concerned about subset)
  mf.call[[2]] <- stats::update(eval(mf.call[[2]]), . ~ . + flexida::assigned())
  mf <- eval(mf.call, parent.frame())

  # Get weights
  lm.call$weights <- mf$"(weights)"

  # get response
  mr <- stats::model.response(mf)

  areg.center <- function(var, grp, wts = NULL) {
    if (!is.null(wts)) {
          df <- data.frame(var = var, wts = wts)
          var - sapply(split(df, grp), function(x) {
            stats::weighted.mean(x$var, x$wts)
          })[as.character(grp)] +
            stats::weighted.mean(var, w = wts, na.rm = TRUE)
    } else {
      var - tapply(var, grp, mean, na.rm = TRUE)[as.character(grp)] +
        mean(var, na.rm = TRUE)
    }
  }

  if (absorb) {
    blocks <- .get_col_from_new_data(design,
                                     eval(lmitt.call$data, parent.frame()),
                                     "b")[, 1]
    # To be used below
  }



  if (!is.null(lm.call$subset)) {
    # If provided, subset the data
    lm.call$data <- subset(eval(lm.call$data),
                              eval(lm.call$subset, eval(lm.call$data)))
  }


  if (rhs == "1") {
    # Define new RHS and obtain model.matrix
    new.form <- formula(~ flexida::assigned()) # need flexida:: or assigned()
                                               # can't be found
#    environment(new.form) <- saveenv # Do I need this?
    mm.call <- lm.call
    mm.call[[2]] <- str2lang(deparse(new.form))
    mm.call[[1]] <- quote(stats::model.matrix.default)
    names(mm.call)[2] <- "object"
    mm <- eval(mm.call, parent.frame())

    if (absorb) {
      mm <- apply(mm, 2, areg.center, as.factor(blocks), lm.call$weights)
    }

  } else {

    # Create model.matrix with subgroup main effects (to BE residualized out)
    sbgrp.form <- stats::reformulate(paste0(rhs, "+ 0"))
    sbgrp.call <- lm.call
    sbgrp.call[[2]] <- str2lang(deparse(sbgrp.form))
    sbgrp.call[[1]] <- quote(stats::model.matrix.default)
    names(sbgrp.call)[2] <- "object"
    sbgrp.mm <- eval(sbgrp.call, parent.frame())

    # Create model.matrix with treatment:subgroup interaction (to be kept in)
    effect.form <- stats::reformulate(paste0("flexida::assigned():", rhs, "+0"))
    effect.call <- lm.call
    effect.call[[2]] <- str2lang(deparse(effect.form))
    effect.call[[1]] <- quote(stats::model.matrix.default)
    names(effect.call)[2] <- "object"
    effect.mm <- eval(effect.call, parent.frame())

    if (absorb) {
      sbgrp.mm <- apply(sbgrp.mm, 2, areg.center, as.factor(blocks),
                        lm.call$weights)
      effect.mm <- apply(effect.mm, 2, areg.center, as.factor(blocks),
                         lm.call$weights)
    }

    # Using `__xx__` to try and ensure no collision with variable names
    mm <- apply(effect.mm, 2, function(xx__) {
      resid.call <- lm.call
      resid.call$offset <- NULL #see issue #101
      resid.call$formula <- stats::reformulate("sbgrp.mm", "xx__")
      eval(resid.call, parent.frame())$resid
    })

  }


  # Strip intercept from model
  mm <- mm[, !grepl("(Intercept)", colnames(mm)), drop = FALSE]

  # Center (weighted if needed)
  .center <- function(x, wts) {
    if (is.null(wts)) {
      x - mean(x, na.rm = TRUE)
    } else {
      x - weighted.mean(x, wts, na.rm = TRUE)
    }
  }

  # Center variables to remove intercept
  mm <- apply(mm, 2, .center, lm.call$weights)

  if (is.matrix(mr)) {
    # if somehow user passes in a matrix outcome, handle centering appropriately
    mr <- apply(mr, 2, .center, lm.call$weights)
  } else {
    mr <- .center(mr, lm.call$weights)
  }

  # Scale terms to remove intercept
  lm.call[[2]] <- reformulate("mm + 0", "mr")

  model <- eval(lm.call, parent.frame())

  return(.convert_to_lmitt(model, design, lmitt_fitted = TRUE))
}

##' @export
##' @rdname lmitt
lmitt.lm <- function(obj,
                     design = NULL,
                     ...) {
  return(as.lmitt(obj, design))
}
