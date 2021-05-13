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
  stats:::print.lm(x, digits = digits, ...)
}

#' @export
setMethod("show", "DirectAdjusted", function(object) {
  stats:::print.lm(object)
})
