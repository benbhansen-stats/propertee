# (Internal) Adopters/weights/cov_adj all need the Designs to operate. If any
# are called in the model without a Design, this function sees if it can find
# the Design in another of these functions.
#
# Note that it will never look inside adopters (gets complicated in formulas),
# only in weights or cov_adj. E.g.
#   lm(y ~ adopters, weights = ate(des), offest = cov_adj(mod1))
#   lm(y ~ adopters, weights = ate(), offest = cov_adj(mod1, design = des))
# will both work, but
#   lm(y ~ adopters(des), weights = ate(), offest = cov_adj(mod1))
# will fail.
# @param NULL_on_error if `TRUE`, returns `NULL` if a Design object is not found
# rather than an exception
# @return a \code{Design} object if one can be found in the call stack,
# otherwise an error or `NULL` depending on `NULL_on_error`
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
      if (is(possible_design_holder, "WeightedDesign") ||
          is(possible_design_holder, "SandwichLayer")) {
        # If we have a WeightedDesign, save it and break
        design <- possible_design_holder@Design
        break()
      }
    }
    if (!is(design, "Design")) {
      return(NULL)
    }
    return(design)
  }


  weight_design <- NULL
  covadj_design <- NULL
  # This avoids infinite recursion; if we're in weights or in cov_adj, don't
  # look for it again. Only adopters will look for both.
  if (sys.call(-1)[[1]] != ".weights_calc") {
    weight_design <- .find.design("weights")
  }
  if (sys.call(-1)[[1]] != "cov_adj") {
    covadj_design <- .find.design("offset")
  }

  # At this point, each *_design is either NULL, or a Design (as enforced by
  # .find.design())

  if (is.null(weight_design) && is.null(covadj_design)) {
    # Found nothing
    if (NULL_on_error) {
      return(NULL)
    }
    stop(paste("Unable to locate Design in call stack, please use the",
               " `design` argument to pass a Design object."))
  }
  if (is(weight_design, "Design") && is(covadj_design, "Design")) {
    # Found both; ensure its the same Design
    if (!identical(weight_design, covadj_design)) {
      stop("`Design`s foundn in both `cov_adj` and weights but differ")
    }
    return(weight_design)
  }

  # Since we know they're not both NULL, and not both Designs, and they have to
  # be one or the other, the only possible scenario is that one is NULL and one
  # is a Design, so return the Design.
  if (is.null(weight_design)) {
    return(covadj_design)
  } else {
    return(weight_design)
  }

}
