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
  return(TRUE)
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
##' @title \code{WeightedDesign} Ops
##' @param e1 \code{WeightedDesign} or \code{numeric}
##' @param e2 \code{numeric} or \code{WeightedDesign}
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
            return(e1)
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("*", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) {
            e2@.Data <- e1 * e2@.Data
            validObject(e2)
            return(e2)
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("/", signature(e1 = "WeightedDesign", e2 = "numeric"),
          function(e1, e2) {
            e1@.Data <- e1@.Data/e2
            validObject(e1)
            return(e1)
          })

##' @rdname WeightedDesignOps
##' @export
setMethod("/", signature(e1 = "numeric", e2 = "WeightedDesign"),
          function(e1, e2) {
            e2@.Data <- e1/e2@.Data
            validObject(e2)
            return(e2)
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
  return(as.numeric(object))
})

setGeneric("subset")

##' @title \code{WeightedDesign} subsetting
##' @param subset Logical vector identifying values to keep or drop
##' @return \code{x} subset by \code{i}
##' @export
##' @rdname WeightedDesign.subset
setMethod("subset", "WeightedDesign", function(x, subset) {
  x@.Data <- subset(x@.Data, subset = subset)
  return(x)
})

setGeneric("[")

##' @param x \code{WeightedDesign} object
##' @param i indices specifying elements to extract or replace. See
##'   \code{help("[")} for further details.
##' @return \code{x} subset by \code{i}
##' @export
##' @importFrom methods callNextMethod
##' @rdname WeightedDesign.subset
setMethod("[", "WeightedDesign",
          function(x, i) {
            dat <- methods::callNextMethod()
            x@.Data <- dat
            return(x)

          })
