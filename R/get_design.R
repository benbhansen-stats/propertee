##' \code{adopters()}/\code{ate()}/\code{ett()}/\code{cov_adj()} all need the
##' \code{Design} to operate. If any are called in the model without a
##' \code{design=} argument, this function sees if it can find the \code{Design}
##' in another of these functions.
##'
##' Note that it will never look inside \code{adopters()} (gets complicated in
##' formulas), only in weights or \code{cov_adj()}. E.g.
##'
##' \code{lm(y ~ adopters(), weights = ate(des), offest = cov_adj(mod1))}
##'
##' \code{lm(y ~ adopters(), weights = ate(), offest = cov_adj(mod1, design = des))}
##'
##' will both work, but
##'
##' \code{lm(y ~ adopters(des), weights = ate(), offest = cov_adj(mod1))}
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
  lmitt_design <- NULL
  # This avoids infinite recursion; if we're in weights or in cov_adj, don't
  # look for it again. Only adopters will look for both.
  if (sys.call(-1)[[1]] != ".weights_calc") {
    weight_design <- .find.design("weights")
  }
  if (sys.call(-1)[[1]] != "cov_adj") {
    covadj_design <- .find.design("offset")
  }
  # The plain `lmitt` call may contain a design formula, only in `lmitt.formula`
  # will that be converted into an actual design, and `expand.model.frame` may
  # have a DirectAdjusted object from which we can extract the design object
  found_lmitt <- grepl("(^lmitt\\.formula$|expand\\.model\\.frame$)",
                       lapply(sys.calls(), "[[", 1),
                       perl = TRUE)
  if (any(found_lmitt)) {
    mf <- which(found_lmitt)[1]
    lmitt_design <- tryCatch(
      get("design", sys.frame(mf)),
      error = function(e1) {
        tryCatch(
          get("model", sys.frame(mf))@Design,
          error = function(e2) {
            NULL
          })
      })
  }
  # If its not a real design, return NULL
  if (!inherits(lmitt_design, "Design")) {
    lmitt_design <- NULL
  }

  # At this point, each *_design is either NULL, or a Design (as enforced by
  # .find.design() and the special lmitt case)

  potential_designs <- list(weight_design,
                            covadj_design,
                            lmitt_design)

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
    if (!Reduce(identical, potential_designs[found_designs])) {
      stop("Multiple differing `Design` found")
    }
  }
  # At this point, we know at least one design is non-NULL, and if there are
  # multiple, they're duplicates, so just take the first real Design.
  return(potential_designs[found_designs][[1]])

}
