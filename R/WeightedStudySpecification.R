#' @include StudySpecification.R
NULL
# The above ensures that `StudySpecification` is defined prior to `WeightedStudySpecification`

#' (Internal) Modeling weights with an accompanying StudySpecification
#'
#' @slot .Data numeric vector of modeling weights
#' @slot StudySpecification a StudySpecification
#' @slot target character string, e.g. "ate"
#' @slot weightAlias alias for target appearing in an originating call
#' @slot dichotomy formula describing a treatment/comparison dichotomy
#' @details
#'  `@target` is used for calculation purpose; defining what weight to calculate
#'   `@weightAlias` is only to store the alias used in creation of the 
#'     weights in case we want to report it later.
#' @keywords internal
setClass("WeightedStudySpecification",
         contains = "numeric",
         slots = c(StudySpecification = "StudySpecification",
                   target = "character",
                   weightAlias = "character",
                   dichotomy = "formula"))

setValidity("WeightedStudySpecification", function(object) {
  if (!.isValidWeightTarget(object@target)) {
    return(paste0("@target must be one of [", .listValidWeightTargets(),
                  "]. unknown @target: ", object@target))
  }
  if (!.isValidWeightAlias(object@weightAlias)) {
    return(paste0("@weightAlias must be one of [", .listValidWeightAliases(),
                 "]. unknown @weightAlias: ", object@weightAlias))
  }
  return(TRUE)
})

##' @title Show a \code{WeightedStudySpecification}
##'
##' @description Prints out the weights from a \code{WeightedStudySpecification}
##'
##' @param object a \code{WeightedStudySpecification} object
##' @return an invisible copy of \code{object}
##' @export
setMethod("show", "WeightedStudySpecification", function(object) {
  print(object@.Data)
  invisible(object)
})


##' @title \code{WeightedStudySpecification} Operations
##'
##' @description Algebraic operators on \code{WeightedStudySpecification}
##'   objects and numeric vectors. \code{WeightedStudySpecification}s do not
##'   support addition or subtraction.
##'
##' @details These are primarily used to either combine weights via
##'   multiplication, or to invert weights. Addition and subtraction are not
##'   supported and will produce errors.
##'
##' @param e1,e2 \code{WeightedStudySpecification} or \code{numeric} objects
##' @rdname WeightedStudySpecificationOps
##' @return a \code{WeightedStudySpecification} object
##' @export
setMethod("+", signature(e1 = "WeightedStudySpecification", e2 = "numeric"),
          function(e1, e2) addsubtracterror()
          )

##' @rdname WeightedStudySpecificationOps
##' @export
setMethod("+", signature(e1 = "numeric", e2 = "WeightedStudySpecification"),
          function(e1, e2) addsubtracterror()
          )

##' @rdname WeightedStudySpecificationOps
##' @export
setMethod("-", signature(e1 = "WeightedStudySpecification", e2 = "numeric"),
          function(e1, e2) addsubtracterror()
          )

##' @rdname WeightedStudySpecificationOps
##' @export
setMethod("-", signature(e1 = "numeric", e2 = "WeightedStudySpecification"),
          function(e1, e2) addsubtracterror()
          )

##' @rdname WeightedStudySpecificationOps
##' @export
setMethod("*", signature(e1 = "WeightedStudySpecification", e2 = "numeric"),
          function(e1, e2) {
            e1@.Data <- e1@.Data * e2
            validObject(e1)
            return(e1)
          })

##' @rdname WeightedStudySpecificationOps
##' @export
setMethod("*", signature(e1 = "numeric", e2 = "WeightedStudySpecification"),
          function(e1, e2) {
            e2@.Data <- e1 * e2@.Data
            validObject(e2)
            return(e2)
          })

##' @rdname WeightedStudySpecificationOps
##' @export
setMethod("/", signature(e1 = "WeightedStudySpecification", e2 = "numeric"),
          function(e1, e2) {
            e1@.Data <- e1@.Data/e2
            validObject(e1)
            return(e1)
          })

##' @rdname WeightedStudySpecificationOps
##' @export
setMethod("/", signature(e1 = "numeric", e2 = "WeightedStudySpecification"),
          function(e1, e2) {
            e2@.Data <- e1/e2@.Data
            validObject(e2)
            return(e2)
          })

addsubtracterror <- function() {
  stop("Cannot perform addition or subtraction on WeightedStudySpecifications")
}

setGeneric("weights")

##' @title Extract Weights from \code{WeightedStudySpecification}
##'
##' @description A \code{WeightedStudySpecification} object contains a numeric
##'   vector with a few additional slots, this extracts only the numeric vector.
##'
##' @param object a \code{WeightedStudySpecification} object
##' @param ... Ignored
##' @return A numeric \code{vector} of the weights
##' @export
setMethod("weights", "WeightedStudySpecification", function(object, ...) {
  return(as.numeric(object))
})

setGeneric("subset")

##' @title \code{WeightedStudySpecification} subsetting
##'
##' @description Provides functionality to subset the weights of a
##'   \code{WeightedStudySpecification} object.
##'
##' @param subset Logical vector identifying values to keep or drop
##' @return A \code{WeightedStudySpecification} object which is a subsetted
##'   version of \code{x}.
##' @export
##' @rdname WeightedStudySpecification.subset
setMethod("subset", "WeightedStudySpecification", function(x, subset) {
  x@.Data <- subset(x@.Data, subset = subset)
  return(x)
})

setGeneric("[")

##' @param x \code{WeightedStudySpecification} object
##' @param i indices specifying elements to extract or replace. See
##'   \code{help("[")} for further details.
##' @export
##' @importFrom methods callNextMethod
##' @rdname WeightedStudySpecification.subset
setMethod("[", "WeightedStudySpecification",
          function(x, i) {
            dat <- methods::callNextMethod()
            x@.Data <- dat
            return(x)

          })
