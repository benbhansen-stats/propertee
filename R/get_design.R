##' \code{assigned()}/\code{ate()}/\code{ett()}/\code{cov_adj()} all need the
##' \code{Design} to operate. If any are called in the model without a
##' \code{design=} argument, this function sees if it can find the \code{Design}
##' in another of these functions.
##'
##' Note that it will never look inside \code{assigned()} (gets complicated in
##' formulas), only in weights or \code{cov_adj()}. E.g.
##'
##' \code{lm(y ~ assigned(), weights = ate(des), offest = cov_adj(mod1))}
##'
##' \code{lm(y ~ assigned(), weights = ate(), offest = cov_adj(mod1, design = des))}
##'
##' will both work, but
##'
##' \code{lm(y ~ assigned(des), weights = ate(), offest = cov_adj(mod1))}
##'
##' will fail.
##' @title (Internal) Locate a \code{Design} in the call stack
##' @param NULL_on_error if \code{TRUE}, returns \code{NULL} if a \code{Design}
##'   object is not found rather than throwing an error.
##' @return A \code{Design}, or \code{NULL} if \code{NULL_on_error} is
##'   \code{TRUE} and the \code{Design} can't be found.
##' @keywords internal
.get_design <- function(NULL_on_error = FALSE) {
  design <- NULL

  # The plain `lmitt` call may contain a design formula, only in `lmitt.formula`
  # will that be converted into an actual design. Want to return this immediately
  # to prevent infinite recursion
  found_lmitt <- grepl("^lmitt\\.formula$",
                       lapply(sys.calls(), "[[", 1),
                       perl = TRUE)

  if (any(found_lmitt)) {
    design <- get("design", sys.frame(which(found_lmitt)[1]))
    if (inherits(design, "Design")) {
      return(design)
    } else {
      design <- NULL
    }
  }

  # #100 Checking for recusion and exiting early.
  if (sum(lapply(sys.calls(), "[[", 1) == ".get_design") > 2) {
    stop(paste("Unable to locate Design in call stack, please use the",
               "`design` argument to pass a Design object.",
               "(Internal: Inf Recusion error)"))
  }


  # Searching for weights or cov_adj is basically the same, except for argument
  # type
  .find.design <- function(type) {
    stopifnot(type %in% c("weights", "offset"))

    design <- NULL

    # Identify all frames with the appropriate argument
    keyframes <- !vapply(lapply(sys.calls(), `[[`, type), is.null, logical(1))

    # Loop over each frame which has an `type` argument.
    # Its most likely the first frame, but perhaps not.
    for (i in which(keyframes)) {
      possible_design_holder <- get(type, sys.frame(i))
      if (inherits(possible_design_holder, "WeightedDesign") ||
          inherits(possible_design_holder, "SandwichLayer")) {
        # If we have a WeightedDesign, save it and break
        design <- possible_design_holder@Design
        break()
      }
    }
    if (!inherits(design, "Design")) {
      return(NULL)
    }
    return(design)
  }

  weight_design <- NULL
  covadj_design <- NULL
  mf_design <- NULL
  # This avoids infinite recursion; if we're in weights or in cov_adj, don't
  # look for it again. Only assigned will look for both.
  if (sys.call(-1)[[1]] != ".weights_calc") {
    weight_design <- .find.design("weights")
  }
  if (sys.call(-1)[[1]] != "cov_adj") {
    covadj_design <- .find.design("offset")
  }

  # `expand.model.frame` or `model.frame` may have a DirectAdjusted object from
  # which we can extract the design object
  mf_funcs <- c("expand\\.model\\.frame$", "model.frame$")
  mf_calls <- sapply(mf_funcs,
                     function(x) grepl(x, lapply(sys.calls(), "[[", 1), perl = TRUE),
                     simplify = FALSE)

  mf_design <- mapply(function(search_func, locs) {
    des <- NULL
    if (any(locs)) {
      get_func <- switch(
        search_func,
        "expand\\.model\\.frame$" = function(mf) {
          mod <- get("model", sys.frame(mf))
          if (inherits(mod, "DirectAdjusted")) mod@Design else NULL
        },
        "model.frame$" = function(mf) {
          form <- get("formula", sys.frame(mf))
          if (inherits(form, "DirectAdjusted")) {
            form@Design
          } else if (inherits(form, "terms") | inherits(form, "formula")) {
            tryCatch(get("design", environment(form)),
                     error = function(e) NULL)
          } else {
            NULL
          }
        }
      )
      des <- get_func(which(locs)[1])
      if (!inherits(des, "Design")) des <- NULL
    }
    des
  }, names(mf_calls), mf_calls, SIMPLIFY = FALSE)

  # At this point, each *_design is either NULL, or a Design (as enforced by
  # .find.design() and the special lmitt case)
  potential_designs <- append(mf_design,
                              c(weight_design = weight_design,
                                covadj_design = covadj_design))

  found_designs <- !vapply(potential_designs, is.null, logical(1))

  if (!any(found_designs)) {
    # Found nothing
    if (NULL_on_error) {
      return(NULL)
    }
    stop(paste("Unable to locate Design in call stack, please use the",
               " `design` argument to pass a Design object."))
  }
  if (sum(found_designs) > 1) {
    # Found multiple Designs; ensure they're all the same
    if (sum(!duplicated(potential_designs[found_designs])) > 1) {
      stop("Multiple differing `Design` found")
    }
  }
  # At this point, we know at least one design is non-NULL, and if there are
  # multiple, they're duplicates, so just take the first real Design.
  return(potential_designs[found_designs][[1]])

}
