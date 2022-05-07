#' @include Design.R
NULL
# The above ensures that `Design` is defined prior to `CovAdjPrediction`

CovAdjPrediction <- setClass("CovAdjPrediction",
                             contains = "numeric",
                             slots = c(Design = "Design"))

setValidity("CovAdjPrediction", function(object) {
  TRUE
})

##' @title Show a CovAdjPrediction
##' @param object CovAdjPredictionDesign object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "CovAdjPrediction", function(object) {
  print(object@.Data)
  invisible(object)
})
