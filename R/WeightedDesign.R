WeightedDesign <- setClass("WeightedDesign",
                           contains = "numeric",
                           slots = c(Design = "Design",
                                     target = "character"))

setValidity("WeightedDesign", function(object) {
  if (any(object[!is.na(object)] < 0)) {
    return("Weights must be non-negative")
  }
  if (all(is.na(object)) || all(object[!is.na(object)] == 0)) {
    return("At least one weight must be positive")
  }
  if (!object@target %in% c("ate", "ett")) {
    return(paste("@target must be one of [ate,ett]. unknown @target:",
                 object@target))
  }
  # The design in weighted design must have an accessible binary treatment.
  if (!has_binary_treatment(object@Design)) {
    return("Treatment must be binary or have a dichotomization.")
  }
  TRUE
})

##' @title Show a WeightedDesign
##' @param object WeightedDesignDesign object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "WeightedDesign", function(object) {
  print(object@.Data)
  invisible(object)
})


##' \code{WeightedDesign}s do not support addition or subtraction, but do
##' support all other reasonable operations.
##'
##' Concatenating \code{WeightedDesign}s with \code{c()} requires both
##' individual \code{WeightedDesign}s to come from the same \code{Design} and
##' have the target (\code{ate()} or \code{ett()}). Both arguments to \code{c()}
##' must be \code{WeightedDesign}.
##'
##' @title \code{WeightedDesign} Ops
##' @param e1 \code{WeightedDesign} or numeric
##' @param e2 numeric or \code{WeightedDesign}
##' @param x \code{WeightedDesign}
##' @param ... additional \code{WeightedDesign}s
##' @rdname WeightedDesignOps
##' @export
setMethod("+", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) addsubtracterror()
          )

##' @rdname WeightedDesignOps
##' @export
setMethod("+", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) addsubtracterror()
          )

##' @rdname WeightedDesignOps
##' @export
setMethod("-", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) addsubtracterror()
          )

##' @rdname WeightedDesignOps
##' @export
setMethod("-", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) addsubtracterror()
          )

##' @rdname WeightedDesignOps
##' @export
setMethod("*", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) {
            e1@.Data <- e1@.Data * e2
            validObject(e1)
            e1
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("*", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) {
            e2@.Data <- e1 * e2@.Data
            validObject(e2)
            e2
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("/", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) {
            e1@.Data <- e1@.Data/e2
            validObject(e1)
            e1
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("/", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) {
            e2@.Data <- e1/e2@.Data
            validObject(e2)
            e2
          })

addsubtracterror <- function() {
  stop("Cannot perform addition or subtraction on WeightedDesigns")
}


setGeneric("weights")

##' @title Extract Weights from WeightedDesign
##' @param object WeightedDesign object
##' @param ... Ignored
##' @return Weights
##' @export
setMethod("weights", "WeightedDesign", function(object, ...) {
  as.numeric(object)
})

##' @rdname WeightedDesignOps
##' @export
##' @importFrom methods slot
setMethod("c", signature(x = "WeightedDesign"),
          function(x, ...) {
            dots <- list(...)
            # x must be a WeightedDesign to get here; ensure all other elements
            # are as well
            if (any(vapply(dots, function(k) !is(k, "WeightedDesign"), TRUE))) {
              stop("WeightedDesigns can only be combined with other WeightedDesigns")
            }

            # Make sure all Designs are the same. Create a list of all targets,
            # then check if `unique` returns a single element
            designs <- c(x@Design, lapply(dots, methods::slot, "Design"))
            if (length(unique(designs)) > 1) {
              stop("WeightedDesigns can only be concatenated from identical Designs")
            }

            # Same with targets
            targets <- c(x@target, lapply(dots, methods::slot, "target"))
            if (length(unique(targets)) > 1) {
              stop("WeightedDesigns can only be concatenated with the same target (ate or ett)")
            }

            # The `Reduce` call c's together all the .Data in ... . The `init`
            # argument adds the .Data to it.
            new("WeightedDesign",
                Reduce(c, lapply(dots, methods::slot , ".Data"), init = x@.Data),
                Design = x@Design, target = x@target)

          })
