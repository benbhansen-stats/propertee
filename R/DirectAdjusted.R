DirectAdjusted <- setClass("DirectAdjusted",
                           contains = "lm",
                           slots = c(Design = "Design",
                                     target = "character"))

#' @export
print.DirectAdjusted <- function(x, digits = max(3L, getOption("digits") - 3L), ...)  {
  stats:::print.lm(x, digits = digits, ...)
}

#' @export
setMethod("show", "DirectAdjusted", function(object) {
  stats:::print.lm(object)
})
