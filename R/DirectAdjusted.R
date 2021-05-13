DirectAdjusted <- setClass("DirectAdjusted",
                           contains = "lm",
                           slots = c(Design = "Design",
                                     target = "character"))

setValidity("DirectAdjusted", function(object) {
  if (length(object@target) != 1 || !object@target %in% c("ett", "ate")) {
    return(paste("@target must be one of [ett, ate]. unknown @target:",
                 paste(object@target, collapse = " ")))
  }
  TRUE
})

##' @title Show an DirectAdjusted
##' @param object DirectAdjusted object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "DirectAdjusted", function(object) {
  print(as(object, "lm"))
  invisible(object)
})

##' @title Convert lm object into DirectAdjusted
##' @param x lm object with weights containing a WeightedDesign
##' @param ... Add'l arguments
##' @return DirectAdjusted
##' @export
as.DirectAdjusted <- function(x, ...) {
  if (!is(x, "lm")) {
    stop("input must be lm object")
  }

  if (!is(x$model$"(weights)", "WeightedDesign")) {
    stop("input model must contain WeightedDesign weights")
  }

  DirectAdjusted(x,
                 Design = x$model$"(weights)"@Design,
                 target = x$model$"(weights)"@target)
}
