##' @title Linear Model for Intention To Treat
##'
##' @description Generates a linear model object to estimate a treatment effect,
##'   with proper estimation of variances accounting for the study design.
##'
##' @details The first argument to [lmitt()] should be a formula specifying the
##'   outcome on the left hand side. The right hand side of the formula can be
##'   any of the following:
##'
##' \itemize{
##'   \item \code{1}: Estimates a main treatment effect.
##'   \item a subgroup variable: Estimates a treatment effect within each level
##'        of your subgrouping variable.
##'   \item a continuous moderator: Estimates a main treatment effect as well as
##'        a treatment by moderator interaction. The moderator is not
##'        automatically centered.
##' }
##'
##' Alternatively, \code{obj} can be a pre-created \code{lm} object. No
##' modification is made to the formula of the object. See the help for
##' \code{as.lmitt()} for details of this conversion.
##'
##' Note that although the \code{Design} creation functions (e.g.
##' [rct_design()]) take an optional \code{subset=} argument used in the
##' creation of the \code{Design}, this is \bold{not} the same as the
##' \code{subset=} argument passed to [lm()] or [lmitt()]. The \code{subset=}
##' argument when creating a \code{Design} restricts the data used to generate
##' the \code{Design}, but has no direct impact on the future \code{lm()} or
##' \code{lmitt()} calls using that \code{Design}. (It can indirectly have an
##' impact by excluding particular units of assignment from receiving a
##' treatment assignment and thus complete case analysis removes them from the
##' model.)
##'
##' On the other hand, the \code{subset=} argument in [lm()] or [lmitt()] refers
##' only to subsetting the \code{data} argument passed into [lm()] or [lmitt()].
##'
##' To avoid variable name collision, the treatment variable defined in the
##' \code{design} will have a "\code{.}" appended to it. For example, if you
##' request a main treatment effect (with a formula of \code{~ 1}) with a
##' treatment variable named "txt", you can obtain it's estimate from the
##' returned \code{teeMod} object via \code{$coefficients["txt."]}.
##'
##' [lmitt()] will produce a message if the \code{Design} passed in has block
##' information that is not being utilized in the model. Note that this is
##' \emph{not} an error, but could be an oversight. To disable this message, run
##' \code{options("propertee_message_on_unused_blocks" = FALSE)}.
##'
##' @param obj A \code{formula} or a \code{lm} object. See Details.
##' @param design The \code{Design} to be used. Alternatively, a formula
##'   creating a design (of the type of that would be passed as the first
##'   argument to [rd_design()], [rct_design()], or [obs_design()]). If the
##'   formula includes a [forcing()] element, an RD design is created. Otherwise
##'   an observational design is created. An RCT design must be created manually
##'   using [rct_design()].
##' @param data A \code{data.frame} such as would be passed into [lm()].
##' @param absorb If \code{TRUE}, fixed effects are included for blocks
##'   identified in the \code{Design}. Excluded in \code{FALSE}. Default is
##'   \code{FALSE}. The estimates of these fixed effects are suppressed from the
##'   returned object.
##' @param offset Offset of the kind which would be passed into [lm()]. Ideally,
##'   this should be the output of [cov_adj()].
##' @param weights Which weights should be generated? Options are \code{"ate"}
##'   or \code{"ett"}. Alternatively, the output of a manually run \code{ate()}
##'   or \code{ett()} can be used.
##' @param ... Additional arguments passed to [lm()]. One such argument is
##' \code{dichotomy}, which can be used to dichotomize a non-binary treatment
##' variable in \code{design}. See the Details section of the \code{ett()} or
##' \code{att()} help pages for information on specifying this formula.
##' @return \code{teeMod} model.
##' @export
##' @importFrom stats lm predict weights weighted.mean reformulate residuals
##' @rdname lmitt
##' @examples
##' data(simdata)
##' des <- rct_design(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
##' mod1 <- lmitt(y ~ 1, data = simdata, design = des, weights = "ate")
##' mod2 <- lmitt(y ~ as.factor(o), data = simdata, design = des, weights = "ate")
##' mod3 <- lmitt(y ~ 1, data = simdata,
##'               design = z ~ uoa(uoa1, uoa2) + forcing(force))
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
  lm.call$subset <- NULL

  ### Allow users to pass in "ate" and "ett" rather than functions if they have
  ### no special modifications/additional arguments
  dichotomy <- lmitt.call$dichotomy
  wt <- lmitt.call$weights
  if (is(wt, "character")) {
    if ((wt_call <- tolower(wt)) %in% c("ate", "ett")) {
      lm.call$weights <- lmitt.call$weights <- call(wt_call)
    } else {
      warning(paste("Character other than \"ate\" or \"ett\" passed to",
                    "`weights=` argument.\nIf you are trying to pass a",
                    "character to the internal `lm` you can disregard this",
                    "warning.\nIf you are attemping to use `propertee`\'s",
                    "weight generation, only \"ate\" and \"ett\" are",
                    "accepted."))
    }
  } else if (is.call(wt) & is.null(dichotomy)) dichotomy <- wt$dichotomy

  # First, make sure we have a valid `design=` - if given a formula, make a new
  # `Design`, otherwise ensure `design=` is `Design` class.
  if (!missing(design)) {
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
                           ", data = ", deparse(lmitt.call$data), ")")
      # str2lang converts character into call
      design <- eval(str2lang(new_d_call))
    } else if (is(design, "WeightedDesign")) {
      design <- design@Design
    } else if (!is(design, "Design")) {
      stop(paste("`design=` must be an object created by `*_design`",
                 "function, or a formula specifying such a design"))
    }

    # #126 block on factor treatments
    if (is.factor(treatment(design)[, 1])) {
      if (is.ordered(treatment(design)[, 1])) {
        fact_or_ord <- "Ordered"
      } else {
        fact_or_ord <- "Factor"
      }
      stop(paste(fact_or_ord, "treatment variables are not yet supported, use",
                 "`dichotomy=` to define a binary treatment."))
    }
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
  .check_for_assigned_and_aliases <- function(fn) {
    if (grepl(paste0(fn, "\\("), rhs)) {
      stop(paste("Do not specify `assigned()` or any of its aliases in",
                 "the right hand side of `lmitt()`.\nTo estimate only",
                 "a treatment effect, pass `~ 1` as the right hand side."))
    }
  }
  .check_for_assigned_and_aliases("assigned")
  .check_for_assigned_and_aliases("adopters")
  .check_for_assigned_and_aliases("a\\.")
  .check_for_assigned_and_aliases("z\\.")

  # `rhs` is now either "1" or subgrouping variable

  # Get weights and offset if passed
  lm.call$weights <- eval.parent(lm.call$weights)
  lm.call$offset <- eval.parent(lm.call$offset)


  # Ensure same design is used in weights and offset, if they're there
  wtdes <- try(lm.call$weights@Design, silent = TRUE)
  ofdes <- try(lm.call$offset@Design, silent = TRUE)

  if (missing(design) & is(wtdes, "Design")) {
    stop(paste("You've passed a `Design` into the weight function (`ate()`",
                "or `ett()`) but not the `lmitt()` call. Please pass the",
                "`Design` into the `design=` argument of `lmitt()`. It is",
                "not needed in `ate()` or `ett()` when passed as the",
                "`design=` argument to `lmitt()`."))
  }

  if (missing(design) & is(ofdes, "Design")) {
    stop(paste("You've passed a `Design` into the offset function",
                " (`cov_adj()`) but not the `lmitt()` call. Please pass the",
                "`Design` into the `design=` argument of `lmitt()`. It is",
                "not needed in `cov_adj()` when passed as the `design=`",
                "argument to `lmitt()`."))
  }

  if (is(wtdes, "Design")) {
    if (!identical_Designs(design, wtdes)) {
      stop(paste("Multiple differing `Design` found (`design` argument to",
                 " `lmitt` and `design` object inside the weights differ)."))
    }
  }

  if (is(ofdes, "Design")) {
    if (!identical_Designs(design, ofdes)) {
      stop(paste("Multiple differing `Design` found (`design` argument to",
                 " `lmitt` and `design` object inside the offset differ)."))
    }
  }

  if (!is(lm.call$weights, "WeightedDesign") & !absorb) {
    if ("b" %in% design@column_index) {
      if (options()$propertee_message_on_unused_blocks) {
        message(paste("The Design object contains block-level information,",
                      "but it is not used in this model. Block information",
                      "is used when weights are defined via `ate()` or `ett()`",
                      "or if the `absorb=TRUE` argument is passed."))
      }
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

  # Identify whether RHS is intercept, continuous moderator, or subgroup
  if (rhs != "1") {
    new.form <- reformulate(rhs, intercept = FALSE)
    modfinder.call <- lm.call
    modfinder.call[[2]] <- str2lang(deparse(new.form))
    modfinder.call[[1]] <- quote(stats::model.matrix)
    names(modfinder.call)[2] <- "object"
    numcol <- ncol(eval(modfinder.call, parent.frame()))
    if (numcol > 1) {
      rhstype <- "categorical"
    } else {
      rhstype <- "continuous"
    }
  } else {
    rhstype <- "intercept"
  }

  # Generate formula for the internal `lm`
  if (rhstype == "intercept") {
    new.form <- formula(~ a.())
    moderator <- character()
  } else if (rhstype == "categorical") {
    new.form <- stats::reformulate(paste0("a.():", rhs, "+", rhs))
    moderator <- rhs
  } else {
    new.form <- stats::reformulate(paste0("a.() + a.():", rhs, "+", rhs))
    moderator <- rhs
  }
  mm.call <- lm.call
  mm.call[[2]] <- str2lang(deparse(new.form))
  # model.matrix.lm supports as `na.action` argument where
  # model.matrix.default doesn't
  mm.call[[1]] <- quote(stats::model.matrix.lm)
  mm.call$na.action <- "na.pass"
  names(mm.call)[2] <- "object"
  mm <- eval(mm.call, parent.frame())

  if (absorb) {
    mm <- areg.center(mm, as.factor(blocks), lm.call$weights)
  }

  # Strip intercept from data if it's in there
  mm <- mm[, !grepl("(Intercept)", colnames(mm)), drop = FALSE]
  # Make sure to keep it as a named matrix

  # Rename `a.()` to treatment named, adding a "." after to avoid conflict with
  # original treatment variable.
  colnames(mm) <- gsub("a\\.\\(\\)",
                       paste0(var_names(design, "t"), "."),
                       colnames(mm))

  # Replace : with _ for interaction to try and avoid backticks
  colnames(mm) <- gsub("\\:", "_", colnames(mm))
  # This isn't foolproof as the sbgrp variable can cause them too
  # e.g. `as.factor(sbgrp)` forces backticks.

  # Rebuild RHS of formula: RHS is "txt_" or "txt_sbgrp", where "txt" has been
  # replaced with the actual treatment name above.
  lm_form <- eval(lm.call$formula, envir = parent.frame())
  lm.call$formula <- str2lang(paste(lm_form[[2]], " ~ ",
                                    paste(
                                      paste0("`", colnames(mm), "`"),
                                      collapse = "+"),
                                    collapse = ""))

  # If the moderators already exist in the passed-in data, remove them to avoid
  # either a) overloading variable names, or b) an error if data is a
  # `grouped_df` (see #137)
  data <- data[, !(names(data) %in% colnames(mm))]

  # Data for model should be original data, plus updated RHS (mm)
  lm.call$data <- cbind(data, mm)

  # restore subset
  lm.call$subset <- savedsubset

  model <- eval(lm.call, parent.frame())

  # `&&` necessary to return FALSE immediately if not enough frames on stack
  if (sys.nframe() >= 2 &&
      !is.null(sys.call(-1)) &&
      sys.call(-1)[[1]] == as.name("lmitt")) {
    # If we're in `lmitt.formula()` via `lmitt()`, save that call.
    lmitt_call <- sys.call(-1)
  } else {
    # Otherwise save the `lmitt.formula()` call
    lmitt_call <- sys.call()
  }
  return(.convert_to_lmitt(model,
                           design,
                           lmitt_fitted = TRUE,
                           moderator = moderator,
                           absorbed_intercepts = absorb,
                           lmitt_call = lmitt_call))
}

##' @export
##' @rdname lmitt
lmitt.lm <- function(obj,
                     design = NULL,
                     ...) {
  return(as.lmitt(obj, design))
}
