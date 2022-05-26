#' @include Design.R
NULL
# The above ensures that `Design` is defined prior to `WeightedDesign`

setClass("WeightedDesign",
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
  if (!is_binary_or_dichotomized(object@Design)) {
    return("Treatment must be binary or have a dichotomy.")
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
##' individual \code{WeightedDesign}s to come from the same \code{Design}
##' (except \code{dichotomy}, see below) and have the target (e.g all created
##' with \code{ate()} or all created with \code{ett()}, no mixing-and-matching).
##' All arguments to \code{c()} must be \code{WeightedDesign}.
##'
##' One exception is when concatenting \code{WeightedDesign}s with the same
##' \code{Design} but different dichotomies. There may be cases where the
##' treatment is continuous or has multiple levels, and there is a need to
##' combine the weights from the same general design, but with different
##' dichotomys. Therefore multiple \code{WeightedDesign}s can be combined if
##' they are identical except for their \code{@dichotomy} slots. The resulting
##' object will be a \code{CombinedWeightedDesign} which tracks all individual
##' \code{dichotomy}.
##'
##' @title \code{WeightedDesign} Ops
##' @param e1 \code{WeightedDesign} or \code{numeric}
##' @param e2 \code{numeric} or \code{WeightedDesign}
##' @param x \code{WeightedDesign}
##' @param ... additional \code{WeightedDesign}sx
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

##' @title Extract Weights from \code{WeightedDesign}
##' @param object \code{WeightedDesign} object
##' @param ... Ignored
##' @return vector of weights
##' @export
setMethod("weights", "WeightedDesign", function(object, ...) {
  as.numeric(object)
})
