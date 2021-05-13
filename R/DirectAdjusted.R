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

#' @export
print.DirectAdjusted <- function(x, digits = max(3L, getOption("digits") - 3L), ...)  {
  print(as(x, "lm"), digits = digits, ...)
  invisible(x)
}

##' @title Show an DirectAdjusted
##' @param object DirectAdjusted object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "DirectAdjusted", function(object) {
  print(as(object, "lm"))
  invisible(object)
})
