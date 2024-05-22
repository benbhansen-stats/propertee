#' @include WeightedDesign.R
NULL
# The above ensures that `WeightedDesign` is defined prior to
# `CombinedWeightedDesign`

setClass("CombinedWeightedDesign",
         contains = "WeightedDesign",
         slots = c(dichotomies = "list",
                   keys = "list"))

setValidity("CombinedWeightedDesign", function(object) {
  if (length(object@dichotomies) != length(object@keys)) {
    return("discrepancy between length of @dichotomies and @keys")
  }
  return(TRUE)
})

##' @title Concatenate weights
##'
##' @description Given several variations of weights generated from a single
##'   \code{Design}, combine into a single weight.
##'
##' @details Concatenating \code{WeightedDesign} objects with [c()] requires
##'   both individual \code{WeightedDesign} objects to come from the same
##'   \code{Design} (except \code{dichotomy}, see below) and have the same
##'   target (e.g all created with [ate()] or all created with [ett()], no
##'   mixing-and-matching). All arguments to [c()] must be
##'   \code{WeightedDesign}.
##'
##'   One exception is when concatenting \code{WeightedDesign} objects whose
##'   \code{Design} differ only in their dichotomies. There may be cases where
##'   the treatment is continuous or has multiple levels, and there is a need to
##'   combine the weights from the same general design, but with different
##'   dichotomies. Therefore multiple \code{WeightedDesign} objects can be
##'   combined if they are identical except for their \code{@dichotomy} slots.
##'   The resulting object will be a \code{WeightedDesign}.
##'
##' @param x, .. a \code{WeightedDesign} object, typically created from [ate()]
##'   or [ett()]
##' @param ... any number of additional \code{WeightedDesign} objects with
##'   equivalent \code{Design} to \code{x} and eachother
##' @param force_dichotomy_equal if \code{FALSE} (default), \code{WeightedDesign}s are
##'   considered equivalent even if their \code{dichotomy} differs. If
##'   \code{TRUE}, \code{@dichotomy} must also be equal.
##' @return A new \code{WeightedDesign} with
##'   the weights concatenated in the input order.
##' @export
##' @importFrom methods slot
##' @importFrom stats formula
##' @examples
##' data(simdata)
##' des <- rct_design(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
##' w1 <- ate(des, data = simdata[1:30,])
##' w2 <- ate(des, data = simdata[31:40,])
##' w3 <- ate(des, data = simdata[41:50,])
##' c_w <- c(w1, w2, w3)
##' c(length(w1), length(w2), length(w3), length(c_w))
##'
##' des <- rct_design(dose ~ unit_of_assignment(uoa1, uoa2), data = simdata)
##' w1 <- ate(des, data = simdata[1:10, ], dichotomy = dose >= 300 ~ .)
##' w2 <- ate(des, data = simdata[11:30, ], dichotomy = dose >= 200 ~ .)
##' w3 <- ate(des, data = simdata[31:50, ], dichotomy = dose >= 100 ~ .)
##' c_w <- c(w1, w2, w3)
setMethod("c", signature(x = "WeightedDesign"),
          function(x, ..., force_dichotomy_equal = FALSE) {
  dots <- list(...)
  # x must be a WeightedDesign to get here; ensure all other elements
  # are as well
  if (any(vapply(dots, function(k) !inherits(k, "WeightedDesign"), TRUE))) {
    stop("WeightedDesigns can only be combined with other WeightedDesigns")
  }

  # Ensure all WeightedDesigns have the same target
  targets <- c(x@target, lapply(dots, methods::slot, "target"))
  if (length(unique(targets)) > 1) {
    stop(paste("WeightedDesigns can only be concatenated with",
               "the same target (ate or ett)"))
  }

  designs <- c(x@Design, lapply(dots, methods::slot, "Design"))
  # temporarily remove calls from the Designs to avoid weird discrepancies here
  destmp <- lapply(designs, function(x) {
    x@call <- call("c") # placeholder for unique check
    x
  })

  if (length(unique(destmp)) == 1) {
    # if all Designs are identical (sans @call), we can c together
    # just the numeric weights, and pull out Design/target from the
    # first one
    return(new("WeightedDesign",
               Reduce(c, lapply(dots, methods::slot, ".Data"), init = x@.Data),
               Design = x@Design,
               target = x@target))
  } else {
    if (force_dichotomy_equal) {
      stop("NEED AN ERROR MESSAGE HERE IF `WeightedDesign` OBECTS HAVE `dichotomy` SLOTS")
    }
    stop("Cannot combine WeightedDesigns from differing Designs")
  }
})
