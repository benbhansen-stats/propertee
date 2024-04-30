#' @include Design.R
NULL
# The above ensures that `Design` is defined prior to `WeightedDesign`

setClass("WeightedDesign",
         contains = "numeric",
         slots = c(Design = "Design",
                   target = "character"))

setValidity("WeightedDesign", function(object) {
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

##' @title Show a \code{WeightedDesign}
##'
##' @description Prints out the weights from a \code{WeightedDesign}
##'
##' @param object a \code{WeightedDesign} object
##' @return an invisible copy of \code{object}
##' @export
setMethod("show", "WeightedDesign", function(object) {
  print(object@.Data)
  invisible(object)
})


##' @title \code{WeightedDesign} Operations
##'
##' @description Algebraic operators on \code{WeightedDesign} objects and
##'   numeric vectors. \code{WeightedDesign}s do not support addition or
##'   subtraction.
##'
##' @details These are primarily used to either combine weights via
##'   multiplication, or to invert weights. Addition and subtraction are not
##'   supported and will produce errors.
##'
##' @param e1,e2 \code{WeightedDesign} or \code{numeric} objects
##' @rdname WeightedDesignOps
##' @return a \code{WeightedDesign} object
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
##'
##' @description A \code{WeightedDesign} object contains a numeric vector with a
##'   few additional slots, this extracts only the numeric vector.
##'
##' @param object a \code{WeightedDesign} object
##' @param ... Ignored
##' @return A numeric \code{vector} of the weights
##' @export
setMethod("weights", "WeightedDesign", function(object, ...) {
  return(as.numeric(object))
})

setGeneric("subset")

##' @title \code{WeightedDesign} subsetting
##'
##' @description Provides functionality to subset the weights of a
##'   \code{WeightedDesign} object.
##'
##' @param subset Logical vector identifying values to keep or drop
##' @return A \code{WeightedDesign} object which is a subsetted version of
##'   \code{x}.
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
##' @export
##' @importFrom methods callNextMethod
##' @rdname WeightedDesign.subset
setMethod("[", "WeightedDesign",
          function(x, i) {
            dat <- methods::callNextMethod()
            x@.Data <- dat
            return(x)

          })
