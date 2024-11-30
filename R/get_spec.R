##' \code{assigned()}/\code{ate()}/\code{ett()}/\code{cov_adj()} all need the
##' \code{StudySpecification} to operate. If any are called in the model without
##' a \code{specification=} argument, this function sees if it can find the
##' \code{StudySpecification} in another of these functions.
##'
##' Note that it will never look inside \code{assigned()} (gets complicated in
##' formulas), only in weights or \code{cov_adj()}. E.g.
##'
##' \code{lm(y ~ assigned(), weights = ate(spec), offest = cov_adj(mod1))}
##'
##' \code{lm(y ~ assigned(), weights = ate(), offest = cov_adj(mod1, specification = spec))}
##'
##' will both work, but
##'
##' \code{lm(y ~ assigned(spec), weights = ate(), offest = cov_adj(mod1))}
##'
##' will fail.
##' @title (Internal) Locate a \code{StudySpecification} in the call stack
##' @param NULL_on_error if \code{TRUE}, returns \code{NULL} if a
##'   \code{StudySpecification} object is not found rather than throwing an
##'   error.
##' @return A \code{StudySpecification}, or \code{NULL} if \code{NULL_on_error}
##'   is \code{TRUE} and the \code{StudySpecification} can't be found.
##' @keywords internal
.get_spec <- function(NULL_on_error = FALSE) {
  specification <- NULL

  # The plain `lmitt` call may contain a specification formula, only in `lmitt.formula`
  # will that be converted into an actual specification. Want to return this immediately
  # to prevent infinite recursion
  found_lmitt <- grepl("^lmitt\\.formula$",
                       lapply(sys.calls(), "[[", 1),
                       perl = TRUE)

  if (any(found_lmitt)) {
    specification <- get("specification", sys.frame(which(found_lmitt)[1]))
    if (inherits(specification, "StudySpecification")) {
      return(specification)
    } else {
      specification <- NULL
    }
  }

  # #100 Checking for recusion and exiting early.
  if (sum(lapply(sys.calls(), "[[", 1) == ".get_spec") > 2) {
    stop(paste("Unable to locate StudySpecification in call stack, please use the",
               "`specification` argument to pass a StudySpecification object.",
               "(Internal: Inf Recusion error)"))
  }


  # Searching for weights or cov_adj is basically the same, except for argument
  # type
  .find.specification <- function(type) {
    stopifnot(type %in% c("weights", "offset"))

    specification <- NULL

    # Identify all frames with the appropriate argument
    keyframes <- !vapply(lapply(sys.calls(), `[[`, type), is.null, logical(1))

    # Loop over each frame which has an `type` argument.
    # Its most likely the first frame, but perhaps not.
    for (i in which(keyframes)) {
      possible_spec_holder <- get(type, sys.frame(i))
      if (inherits(possible_spec_holder, "WeightedStudySpecification") ||
          inherits(possible_spec_holder, "SandwichLayer")) {
        # If we have a WeightedStudySpecification, save it and break
        specification <- possible_spec_holder@StudySpecification
        break()
      }
    }
    if (!inherits(specification, "StudySpecification")) {
      return(NULL)
    }
    return(specification)
  }

  weight_spec <- NULL
  covadj_spec <- NULL
  mf_spec <- NULL
  emf_spec <- NULL
  # This avoids infinite recursion; if we're in weights or in cov_adj, don't
  # look for it again. Only assigned will look for both.
  if (sys.call(-1)[[1]] != ".weights_calc") {
    weight_spec <- .find.specification("weights")
  }
  if (sys.call(-1)[[1]] != "cov_adj" && sys.call(-1)[[1]] != "propertee::cov_adj") {
    covadj_spec <- .find.specification("offset")
  }

  # `model.frame` may have a teeMod object we can extract from
  mf_calls <- grepl("model\\.frame$", lapply(sys.calls(), "[[", 1), perl = TRUE)

  mf_spec <- lapply(which(mf_calls), function(x) {
    form <- get("formula", sys.frame(x))
    found_spec <- if (inherits(form, "teeMod")) {
      form@StudySpecification
    } else if (inherits(form, "terms") | inherits(form, "formula")) {
      tryCatch(get("specification", environment(form)),
               error = function(e) NULL)
    } else {
      NULL
    }
    if (!inherits(found_spec, "StudySpecification")) {
      found_spec <- NULL
    }
    return(found_spec)
  })

  emf_calls <- grepl("expand\\.model\\.frame$",
                     lapply(sys.calls(), "[[", 1), perl = TRUE) |
    grepl("\\.expand\\.model\\.frame_teeMod$",
          lapply(sys.calls(), "[[", 1), perl = TRUE)

  emf_spec <- lapply(which(emf_calls), function(x) {
    mod <- get("model", sys.frame(x))
    if (inherits(mod, "teeMod")) mod@StudySpecification else NULL
  })

  # At this point, each *_spec is either NULL, or a StudySpecification (as enforced by
  # .find.specification() and the special lmitt case)
  potential_specs <- append(emf_spec,
                              append(mf_spec,
                                     c(weight_spec = weight_spec,
                                       covadj_spec = covadj_spec)))

  found_specs <- !vapply(potential_specs, is.null, logical(1))

  if (!any(found_specs)) {
    # Found nothing
    if (NULL_on_error) {
      return(NULL)
    }
    stop(paste("Unable to locate StudySpecification in call stack, please use the",
               " `specification` argument to pass a StudySpecification object."))
  }
  if (sum(found_specs) > 1) {
    # Found multiple StudySpecifications; ensure they're all the same
    if (sum(!duplicated(potential_specs[found_specs])) > 1) {
      stop("Multiple differing `StudySpecification` found")
    }
  }
  # At this point, we know at least one specification is non-NULL, and if there are
  # multiple, they're duplicates, so just take the first real StudySpecification.
  return(potential_specs[found_specs][[1]])

}
