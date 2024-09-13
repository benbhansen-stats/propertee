#' @include WeightedDesign.R
NULL

##' @title Concatenate weights
##'
##' @description Given several variations of weights generated from a single
##'   \code{Design}, combine into a single weight.
##'
##' @details Concatenating \code{WeightedDesign} objects with [c()] requires
##'   both individual \code{WeightedDesign} objects to come from the same
##'   \code{Design} and have the same
##'   target (e.g all created with [ate()] or all created with [ett()], no
##'   mixing-and-matching). All arguments to [c()] must be
##'   \code{WeightedDesign}.
##'
##'   \code{WeightedDesign} objects may be concatenated together even without
##'   having the same \code{@dichotomy} slot. This procedure only prompts a
##'   warning for differing dichotomies if the argument \code{warn_dichotomy_not_equal}
##'   is set to \code{TRUE}.
##' @param x, .. a \code{WeightedDesign} object, typically created from [ate()]
##'   or [ett()]
##' @param ... any number of additional \code{WeightedDesign} objects with
##'   equivalent \code{Design} to \code{x} and eachother
##' @param warn_dichotomy_not_equal if \code{FALSE} (default), \code{WeightedDesign}s are
##'   considered equivalent even if their \code{dichotomy} differs. If
##'   \code{TRUE}, a warning is produced.
##' @return A numeric \code{vector} with
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
          function(x, ..., warn_dichotomy_not_equal = FALSE) {
  dots <- list(...)
  # x must be a WeightedDesign to get here; ensure all other elements
  # are as well
  if (any_numeric <- any(vapply(dots, function(k) !inherits(k, "WeightedDesign"), TRUE))) {
    message(paste("Concatenating a WeightedDesign with a non-WeightedDesign vector.",
                  "If non-WeightedDesign vector has been formed from previous",
                  "concatenation of WeightedDesign objects, equality of `Design`",
                  "slots cannot be confirmed"))
  }

  if (!any_numeric) {
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
      if (warn_dichotomy_not_equal) {
        dichotomies <- Reduce(c, lapply(dots, methods::slot, "dichotomy"), init = x@dichotomy)
        if (length(unique(dichotomies)) > 1) {
          warning("Concatenating WeightedDesigns with different `dichotomy` slots")
        }
      }
    } else {
      stop("Cannot combine WeightedDesigns from differing Designs")
    }
  }

  return(Reduce(c, lapply(dots, methods::slot, ".Data"), init = x@.Data))
})
