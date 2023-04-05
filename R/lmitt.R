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
##' @param design The \code{Design} to be used. Alternatively, a formula
##'   creating a design (of the type of that would be passed as the first
##'   argument to \code{rd_design()}, \code{rct_design()}, or
##'   \code{obs_design()}).
##' @param data Data frame such as would be passed into \code{lm()}.
##' @param absorb If \code{TRUE}, fixed effects are included for units of
##'   assignemnt/clusters identified in the \code{Design}. Excluded in
##'   \code{FALSE}. Default is \code{FALSE}.
##' @param offset Offset of the kind which would be passed into \code{lm()}. To
##'   utilize flexida's functionality, the output of \code{cov_adj()}.
##' @param weights Weight of the kind which would be passed into \code{lm()}. To
##'   utilize flexida's functionality, the output of \code{ate()} or
##'   \code{ett()}, or the strings \code{"ate"}/\code{"ett"}.
##' @param ... Additional arguments passed to \code{lm()}.
##' @return \code{DirectAdjusted} model.
##' @export
##' @importFrom stats lm predict weights weighted.mean reformulate residuals
##' @rdname lmitt
lmitt <- function(obj,
                  design,
                  data,
                  ...) {
  UseMethod("lmitt")
}

##' @export
##' @rdname lmitt
lmitt.formula <- function(obj,
                          design,
                          data,
                          absorb = FALSE,
                          offset = NULL,
                          weights = NULL,
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

  # Save the subset and remove it from `lm.call`, to be re-added at the enc
  savedsubset <- lm.call$subset
  logicalsubset <- eval(savedsubset, data)
  lm.call$subset <- NULL

  ### Allow users to pass in "ate" and "ett" rather than functions if they have
  ### no special modifications/additional arguments
  wt <- lmitt.call$weights
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
      new_d_call <- paste0(new_d_call, ", dichotomy = ",
                           deparse(lmitt.call$dichotomy))
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
    stop(paste("Do not specify `assigned()` in the right hand side of",
               "`lmitt()`.\nTo estimate only a treatment effect, pass",
               "`~ 1` as the right hand side."))
  }
  # `rhs` is now either "1" or subgrouping variable

  # Get weights and offset if passed
  lm.call$weights <- eval.parent(lm.call$weights)
  lm.call$offset <- eval.parent(lm.call$offset)

  # Ensure same design is used in weights and offset, if they're there
  wtdes <- try(lm.call$weights@Design, silent = TRUE)
  ofdes <- try(lm.call$offset@Design, silent = TRUE)

  if (is(wtdes, "Design")) {
    if (!identical(design, wtdes)) {
      stop(paste("Multiple differing `Design` found (`design` argument to",
                 " `lmitt` and `design` object inside the weights differ)."))
    }
  }

  if (is(ofdes, "Design")) {
    if (!identical(design, ofdes)) {
      stop(paste("Multiple differing `Design` found (`design` argument to",
                 " `lmitt` and `design` object inside the offset differ)."))
    }
  }


  # Evaluate model.frame
  mf.call <- lm.call
  mf.call[[1]] <- quote(stats::model.frame)
  # Add assigned() to the model so it utilizes Design characteristics (primarly
  # concerned about subset)
  mf.call[[2]] <- stats::update(eval(mf.call[[2]]), . ~ . + flexida::assigned())
  mf.call$na.action <- "na.pass"
  mf <- eval(mf.call, parent.frame())

  areg.center <- function(var, grp, wts = NULL) {
    if (!is.null(wts)) {
      # weighted.mean produces NA if any weights are NA
      if (any(is.na(wts))) {
        var2 <- var[!is.na(wts)]
      } else {
        var2 <- var
      }

      df <- data.frame(var = var, wts = wts)
      var - sapply(split(df, grp), function(x) {
        stats::weighted.mean(x$var, x$wts, na.rm = TRUE)
      })[as.character(grp)] +
        stats::weighted.mean(var2, w = wts[!is.na(wts)], na.rm = TRUE)
    } else {
      var - tapply(var, grp, mean, na.rm = TRUE)[as.character(grp)] +
        mean(var, na.rm = TRUE)
    }
  }

  if (absorb) {
    if (length(var_names(design, "b")) == 0) {
      stop("No blocks found in Design, cannot absorb")
    }
    blocks <- .get_col_from_new_data(design,
                                     eval(lm.call$data, parent.frame()),
                                     "b", all.x = TRUE)[, 1]
    # To be used below
  }

  if (rhs == "1") {
    # Define new RHS and obtain model.matrix
    new.form <- formula(~ flexida::assigned()) # need flexida:: or assigned()
                                               # can't be found
#    environment(new.form) <- saveenv # Do I need this?
    mm.call <- lm.call
    mm.call[[2]] <- str2lang(deparse(new.form))
    # model.matrix.lm supports as `na.action` argument where
    # model.matrix.default doesn't
    mm.call[[1]] <- quote(stats::model.matrix.lm)
    mm.call$na.action <- "na.pass"
    names(mm.call)[2] <- "object"
    mm <- eval(mm.call, parent.frame())

    if (absorb) {
      mm <- apply(mm, 2, areg.center, as.factor(blocks), lm.call$weights)
    }

    absorbed_moderators = character()

  } else {

    # Create model.matrix with subgroup main effects (to BE residualized out)
    sbgrp.form <- stats::reformulate(paste0(paste(rhs,
                                                  "flexida::assigned()",
                                                  sep = "+"),
                                            "+ 0"))
    sbgrp.call <- lm.call
    sbgrp.call[[2]] <- str2lang(deparse(sbgrp.form))
    sbgrp.call[[1]] <- quote(stats::model.matrix.lm)
    names(sbgrp.call)[2] <- "object"
    sbgrp.call$na.action <- "na.pass"
    sbgrp.mm <- eval(sbgrp.call, parent.frame())
    sbgrp.mm <- sbgrp.mm[, !grepl("assigned\\(", colnames(sbgrp.mm)),
                         drop = FALSE]

    # Create model.matrix with treatment:subgroup interaction (to be kept in)
    effect.form <- stats::reformulate(paste0("flexida::assigned():", rhs, "+0"))
    effect.call <- lm.call
    effect.call[[2]] <- str2lang(deparse(effect.form))
    effect.call[[1]] <- quote(stats::model.matrix.lm)
    effect.call$na.action <- "na.pass"
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
      # By switching from `na.omit` to `na.exclude`, `residuals()` includes NAs
      resid.call$na.action <- "na.exclude"
      stats::residuals(eval(resid.call, parent.frame()))
    })

    absorbed_moderator <- as.character(rhs)
  }


  # Strip intercept from model
  mm <- mm[, !grepl("(Intercept)", colnames(mm)), drop = FALSE]
  # Make sure to keep it as a named matrix

  # Center (weighted if needed)
  .center <- function(x, wts, sbst) {
    # If the user passed in `$subset`, we'll want to center by the mean of the
    # subset, not the whole data.
    if (!is.null(sbst)) {
      xs <- x[sbst]
    } else {
      xs <- x
    }
    # `weighted.mean` with a NULL weights argument calls `mean`, but due to
    # the below issue with NAs in the weights, don't try combining these two
    # calls into one.
    if (is.null(wts)) {
      x - mean(xs, na.rm = TRUE)
    } else {
      # `weighted.mean` with `na.rm = TRUE` only drops `x` as NA, any NA weights
      # will return an NA mean
      x - weighted.mean(xs[!is.na(lm.call$weights)],
                        lm.call$weights[!is.na(lm.call$weights)],
                        na.rm = TRUE)
    }
  }

  # Center variables to remove intercept
  mm <- apply(mm, 2, .center, lm.call$weights, logicalsubset)

  # get response
  mr <- stats::model.response(mf)

  if (is.matrix(mr)) {
    # if somehow user passes in a matrix outcome, handle centering appropriately
    flexida_y <- apply(mr, 2, .center, lm.call$weights, logicalsubset)
  } else {
    flexida_y <- .center(mr, lm.call$weights, logicalsubset)
  }

  toreturn <- as.lmitt(model, design)
  toreturn@lmitt_fitted <- TRUE
  toreturn@absorbed_intercepts <- absorbed_intercepts
  toreturn@absorbed_moderators <- absorbed_moderators
  return(toreturn)

  # Rebuild formula. LHS is updated outcome (`flexida_y`), RHS is 0 (everything
  # is centered) and `flexida::assigned()` or `flexida::assigned():sbrp`
  lm.call$formula <- formula(paste("flexida_y ~ 0 + ",
                                   paste(
                                     paste0("`", colnames(mm), "`"),
                                     collapse = "+"),
                                   collapse = ""))

  # Data for model should be original data, plus updated outcome (flexida_y) and
  # RHS (mm)
  lm.call$data <- cbind(data, flexida_y, mm)

  # restore subset
  lm.call$subset <- savedsubset

  model <- eval(lm.call, parent.frame())

  return(.convert_to_lmitt(model,
                           design,
                           lmitt_fitted = TRUE,
                           absorbed_moderators = absorbed_moderators,
                           absorbed_intercepts = TRUE))
}

##' @export
##' @rdname lmitt
lmitt.lm <- function(obj,
                     design = NULL,
                     ...) {
  return(as.lmitt(obj, design))
}
