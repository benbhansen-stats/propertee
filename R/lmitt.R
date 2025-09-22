##' @title Linear Model for Intention To Treat
##'
##' @description Generates a linear model object to estimate a treatment effect,
##'   with proper estimation of variances accounting for the study
##'   specification.
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
##' The [lmitt()] function's \code{subset=} argument governs
##' the subsetting of \code{data} prior to model fitting, just as with [lm()].
##' Functions such as [rct_spec()] that create \code{StudySpecification}s
##' also take an optional \code{subset=} argument, but its role differs
##' from that of the \code{subset=} argument of [lm()] or [lmitt()]. The \code{subset=}
##' argument when creating a \code{StudySpecification} restricts the data used
##' to generate the \code{StudySpecification}, but has no direct impact on the
##' future \code{lm()} or \code{lmitt()} calls using that
##' \code{StudySpecification}. (It can have an indirect impact by excluding
##' particular units from receiving a treatment assignment or weight. When
##' treatment assignments or weights are reconstructed from the \code{StudySpecification},
##' these units will receive NAs, and will be excluded from the \code{lm()}
##' or \code{lmitt()} fit under typical \code{na.action} settings.)
##'
##' To avoid variable name collision, the treatment variable defined in the
##' \code{specification} will have a "\code{.}" appended to it. For example, if
##' you request a main treatment effect (with a formula of \code{~ 1}) with a
##' treatment variable named "txt", you can obtain its estimate from the
##' returned \code{teeMod} object via \code{$coefficients["txt."]}.
##'
##' [lmitt()] will produce a message if the \code{StudySpecification}
##' designates treatment assignment by block but the blocking
##' structure appears not to be reflected in the \code{weights}, nor
##' in a block fixed effect adjustment (via \code{absorb=TRUE}). While
##' not an error, this is at odds with intended uses of
##' \code{propertee}, so \code{lmitt()} flags it as a potential
##' oversight on the part of the analyst. To disable this message, run
##' \code{options("propertee_message_on_unused_blocks" = FALSE)}.
##'
##' [lmitt()] returns objects of class \sQuote{\code{teeMod}}, for
##' Treatment Effect Estimate Model, extending the lm class to add a
##' summary of the response distribution under control (the
##' coefficients of a controls-only regression of the response on an
##' intercept and any moderator variable).  \code{teeMod} objects also
##' record the underlying \code{StudySpecification} and information
##' about any externally fitted models \code{mod} that may have been
##' used for covariance adjustment by passing
##' \code{offset=cov_adj(mod)}. In the latter case, responses are
##' offsetted by predictions from \code{mod} prior to treatment effect
##' estimation, but estimates of the response variable distribution
##' under control are calculated without reference to \code{mod}.
##'
##' The response distribution under control is also characterized when
##' treatment effects are estimated with block fixed effects, i.e. for
##' [lmitt()] with a \code{formula} first argument with option
##' \code{absorb=TRUE}.  Here as otherwise, the supplementary
##' coefficients describe a regression of the response on an intercept
##' and moderator variables, to which only control observations
##' contribute; but in this case the weights are modified for this
##' supplementary regression. The treatment effect estimates adjusted
##' for block fixed effects can be seen to coincide with estimates
##' calculated without block effect but with weights multiplied by an
##' additional factor specific to the combination of block and
##' treatment condition. For block \eqn{s} containing units with weights
##' \eqn{w_i} and binary treatment assignments \eqn{z_i}, define
##' \eqn{\hat{\pi}_s} by \eqn{\hat{\pi}_s\sum_sw_i=\sum_sz_iw_i}. If
##' \eqn{\hat{\pi}_s} is 0 or 1, the block doesn't contribute to
##' effect estimation and the additional weighting factor is 0; if
##' \eqn{0 < \hat{\pi}_s < 1}, the additional weighting factor is
##' \eqn{1 - \hat{\pi}_s} for treatment group members and
##' \eqn{\hat{\pi}_s} for controls. When estimating a main effect only
##' or a main effect with continuous moderator, supplementary
##' coefficients under option \code{absorb=TRUE} reflect regressions
##' with additional weighting factor equal to 0 or \eqn{\hat{\pi}_s},
##' respectively, for treatment or control group members of block
##' \eqn{s}. With a categorical moderator and \code{absorb=TRUE},
##' this additional weighting factor determining supplementary coefficients
##' is calculated separately for each level \eqn{\ell} of the moderator
##' variable, with the sums defining \eqn{\hat{\pi}_{s\ell}} restricted
##' not only to block \eqn{s} but also to observations with moderator
##' equal to \eqn{\ell}. 
##' 
##' 
##' @param obj A \code{formula} or a \code{lm} object. See Details.
##' @param specification The \code{StudySpecification} to be used.
##'   Alternatively, a formula creating a specification. (Of the type of that
##'   would be passed as the first argument to [rd_spec()], [rct_spec()], or
##'   [obs_spec()], with the difference that [cluster()], [uoa()] 
#'    and [unit_of_assignment()] terms can be omitted when each row of 
#'    \code{data} represents a distinct unit of assignment.) If the formula includes a [forcing()] element, an RD
##'   specification is created. Otherwise an observational specification is
##'   created. An RCT specification must be created manually using [rct_spec()].
##' @param data A \code{data.frame} such as would be passed into [lm()].
##' @param absorb If \code{TRUE}, fixed effects are included for blocks
##'   identified in the \code{StudySpecification}. Excluded in \code{FALSE}.
##'   Default is \code{FALSE}. The estimates of these fixed effects are
##'   suppressed from the returned object.
##' @param offset Offset of the kind which would be passed into [lm()]. Ideally,
##'   this should be the output of [cov_adj()].
##' @param weights Which weights should be generated? Options are \code{"ate"}
##'   or \code{"ett"}. Alternatively, the output of a manually run \code{ate()}
##'   or \code{ett()} can be used.
##' @param ... Additional arguments passed to [lm()] and other functions.
##'   An example of the latter is \code{dichotomy=}, a formula passed to
##'   [assigned()] and, as appropriate, [ate()], [att()], [atc()] or
##'   [ato()]. It is used to dichotomize a non-binary treatment
##'   variable in \code{specification}. See the Details section of the
##'   \code{ate()} help page for examples.
##' @return \code{teeMod} object (see Details)
##' @export
##' @importFrom stats lm predict weights weighted.mean reformulate residuals
##' @rdname lmitt
##' @examples
##' data(simdata)
##' spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
##' mod1 <- lmitt(y ~ 1, data = simdata, specification = spec, weights = "ate")
##' mod2 <- lmitt(y ~ as.factor(o), data = simdata, specification = spec, weights = "ate")
##' ### observational study with treatment z assigned row-wise within blocks:
##' mod3 <- lmitt(y ~ 1, data=simdata, specification=z ~ block(bid), weights="att")
##' ### regression discontinuity study with units of assignment
##' ### given by combinations of uoa1, uoa2:
##' mod4 <- lmitt(y ~ 1, data = simdata,
##'               specification = z ~ uoa(uoa1, uoa2) + forcing(force))
lmitt <- function(obj,
                  specification,
                  data,
                  ...) {
  UseMethod("lmitt")
}

##' @export
##' @rdname lmitt
lmitt.formula <- function(obj,
                          specification,
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
  #lmitt.call contains specification, absorb. lm.call contains only things that get
  #passed to lm, model.matrix

  data <- .as_data_frame(data)

  # Save the subset and remove it from `lm.call`, to be re-added at the enc
  savedsubset <- lm.call$subset
  lm.call$subset <- NULL

  ### Allow users to pass in "ate" and "ett" rather than functions if they have
  ### no special modifications/additional arguments
  dichotomy <- eval.parent(lmitt.call$dichotomy)
  wt <- lmitt.call$weights
  if (is(wt, "character")) {
    if (.isValidWeightAlias(wt_call <- tolower(wt))) {
      lm.call$weights <- lmitt.call$weights <- call(wt_call, dichotomy = dichotomy)
    } else {
      warning(paste0("Character other than [",
                    .listValidWeightAliases(),
                    "] passed to `weights=` argument.\nIf you are trying ",
                    "to pass a character to the internal `lm` you can ",
                    "disregard this warning.\nIf you are attempting to ",
                    "use `propertee`\'s weight generation, only [",
                    .listValidWeightAliases(),
                    "] are accepted."))
    }
  } else if (is.call(wt) & is.null(dichotomy)) {
    for (s in seq_len(sys.nframe())) {
      dichotomy <- tryCatch(eval.parent(wt$dichotomy, s), error = function(e) NULL)
      if (!is.null(dichotomy)) break()
    }
  }

  # First, make sure we have a valid `specification=` - if given a formula, make a new
  # `StudySpecification`, otherwise ensure `specification=` is `StudySpecification` class.
  if (!missing(specification)) {
    if (inherits(specification, "formula")) {
      # If there's a `forcing()`, user wants RDD. If not, force Obs. To do RCT,
      # must create StudySpecification manually.
      if (!is.null(attr(terms(specification, specials = "forcing"),
                        "specials")$forcing)) {
        spec_call <- "rd_spec"
        type <- "RD"
      } else {
        spec_call <- "obs_spec"
        type <- "Obs"
      }

      # Build new call. All calls must include obj and data
      new_spec_call <- paste0(spec_call, "(",
                              "formula = ", deparse1(specification),
                              ", data = ", deparse1(lmitt.call$data), ")")
      new_internal_call <- paste0(".new_StudySpecification(",
                                  "form = ", deparse1(specification),
                                  ", data = ", deparse1(lmitt.call$data),
                                  ", type = ", deparse1(type),
                                  ", call = str2lang(\"", new_spec_call, "\")",
                                  ", called_from_lmitt = TRUE",
                                  ")")
      # str2lang converts character into call
      specification <- eval(str2lang(new_internal_call))
    } else if (is(specification, "WeightedStudySpecification")) {
      specification <- specification@StudySpecification
    } else if (!is(specification, "StudySpecification")) {
      stop(paste("`specification=` must be an object created by `*_spec`",
                 "function, or a formula specifying such a specification"))
    }

    # #126 block on factor treatments
    if (is.factor(treatment(specification)[, 1]) & is.null(dichotomy)) {
      if (is.ordered(treatment(specification)[, 1])) {
        fact_or_ord <- "Ordered"
      } else {
        fact_or_ord <- "Factor"
      }
      stop(paste(fact_or_ord, "treatment variables are not supported, use",
                 "`dichotomy=` to define a binary treatment."))
    }
    
    # #205 block on continuous treatments
    if (!has_binary_treatment(specification) & is.null(dichotomy)) {
      stop("Specify a dichotomy when estimating effects of a continuous treatment variable")
    }
  }


  # Extract formula bits
  rhs <- trimws(strsplit(deparse1(obj[[3]]), "+", fixed = TRUE)[[1]])
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


  # Ensure same specification is used in weights and offset, if they're there
  wtspec <- try(lm.call$weights@StudySpecification, silent = TRUE)
  offspec <- try(lm.call$offset@StudySpecification, silent = TRUE)

  if (missing(specification) & is(wtspec, "StudySpecification")) {
    stop(paste("You've passed a `StudySpecification` into the weight function (`ate()`",
                "or `ett()`) but not the `lmitt()` call. Please pass the",
                "`StudySpecification` into the `specification=` argument of `lmitt()`. It is",
                "not needed in `ate()` or `ett()` when passed as the",
                "`specification=` argument to `lmitt()`."))
  }

  if (missing(specification) & is(offspec, "StudySpecification")) {
    stop(paste("You've passed a `StudySpecification` into the offset function",
                " (`cov_adj()`) but not the `lmitt()` call. Please pass the",
                "`StudySpecification` into the `specification=` argument of `lmitt()`. It is",
                "not needed in `cov_adj()` when passed as the `specification=`",
                "argument to `lmitt()`."))
  }

  if (is(wtspec, "StudySpecification")) {
    if (!identical_StudySpecifications(specification, wtspec)) {
      stop(paste("Multiple differing `StudySpecification` found (`specification` argument to",
                 " `lmitt` and `specification` object inside the weights differ)."))
    }
  }

  if (is(offspec, "StudySpecification")) {
    if (!identical_StudySpecifications(specification, offspec)) {
      stop(paste("Multiple differing `StudySpecification` found (`specification` argument to",
                 " `lmitt` and `specification` object inside the offset differ)."))
    }
  }

  if (!is(lm.call$weights, "WeightedStudySpecification") & !absorb) {
    if ("b" %in% specification@column_index) {
      if (options()$propertee_message_on_unused_blocks) {
        message(paste("The StudySpecification object contains block-level information,",
                      "but it is not used in this model. Block information",
                      "is used when weights are defined via `ate()` or `ett()`",
                      "or if the `absorb=TRUE` argument is passed."))
      }
    }
  }

  blocks <- blocks(specification,
                  eval(lm.call$data, parent.frame()),
                  all.x = TRUE,
                  implicit = TRUE)[,1]

  # Identify whether RHS is intercept, continuous moderator, or subgroup
  if (rhs != "1") {
    new.form <- reformulate(rhs, intercept = FALSE)
    modfinder.call <- lm.call
    modfinder.call[[2]] <- str2lang(deparse1(new.form))
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
  if (rhstype == "intercept" | rhs == var_names(specification, "t")) {
    new.form <- stats::reformulate(paste0("a.(dichotomy=", deparse1(dichotomy), ")"))
    moderator <- character()
    ctrl_means_form <- stats::reformulate("1", response = "lhs")
  } else if (rhstype == "categorical") {
    new.form <- stats::reformulate(paste0("a.(dichotomy=", deparse1(dichotomy), "):",
                                          rhs, "+", rhs))
    moderator <- rhs
    ctrl_means_form <- stats::reformulate(rhs, response = "lhs", intercept = FALSE)
  } else {
    new.form <- stats::reformulate(paste0("a.(dichotomy=", deparse1(dichotomy), ") + ",
                                          "a.(dichotomy=", deparse1(dichotomy), "):",
                                          rhs, "+", rhs))
    moderator <- rhs
    ctrl_means_form <- stats::reformulate(rhs, response = "lhs")
  }
  ctrl_means_env <- new.env()
  assign("ctrl_means_data", data, envir = ctrl_means_env)
  environment(ctrl_means_form) <- ctrl_means_env
  
  mm.call <- lm.call
  mm.call[[2]] <- str2lang(paste(deparse1(new.form, width.cutoff = 500L), collapse = ""))
  # model.matrix.lm supports as `na.action` argument where
  # model.matrix.default doesn't
  mm.call[[1]] <- quote(stats::model.matrix.lm)
  mm.call$na.action <- "na.pass"
  names(mm.call)[2] <- "object"
  mm <- eval(mm.call, parent.frame())
  colnames(mm) <- gsub(paste0("dichotomy = ", deparse1(dichotomy)), "", colnames(mm),
                       fixed = TRUE)

  if (absorb) {
    mm <- areg.center(mm, as.factor(blocks), lm.call$weights)
  }

  # Strip intercept from data if it's in there
  mm <- mm[, !grepl("(Intercept)", colnames(mm)), drop = FALSE]
  # Make sure to keep it as a named matrix

  # Rename `a.()` to treatment named, adding a "." after to avoid conflict with
  # original treatment variable.
  colnames(mm) <- gsub("a\\.\\(\\)",
                       paste0(var_names(specification, "t"), "."),
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
  
  # Data for model should be original data, plus updated RHS (mm),
  lm.call$data <- cbind(data, mm)

  # restore subset
  lm.call$subset <- savedsubset

  model <- eval(lm.call, parent.frame())
  
  # set call's na.action to na.pass so expand.model.frame includes NA rows
  model$call$na.action <- "na.pass"
  
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
                           specification,
                           lmitt_fitted = TRUE,
                           moderator = moderator,
                           absorbed_intercepts = absorb,
                           ctrl_means_form = ctrl_means_form,
                           dichotomy = dichotomy,
                           lmitt_call = lmitt_call))
}

##' @export
##' @rdname lmitt
lmitt.lm <- function(obj,
                     specification = NULL,
                     ...) {
  return(as.lmitt(obj, specification))
}
