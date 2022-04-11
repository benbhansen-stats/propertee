#' @importFrom utils getS3method methods
#' @importFrom stats model.frame terms na.pass
#' @importFrom methods as is new validObject show
NULL
#> NULL


.onLoad <- function(lib, pkg) {
  options("flexida_warn_on_conditional_treatment" = TRUE)
}
