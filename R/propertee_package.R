#' @importFrom utils getS3method methods
#' @importFrom stats model.frame terms na.pass vcov confint
#' @importFrom methods as is new validObject show
#' @importFrom sandwich estfun bread
NULL
#> NULL


.onLoad <- function(lib, pkg) {
  options("propertee_warn_on_conditional_treatment" = TRUE)
  options("propertee_message_on_unused_blocks" = TRUE)
  options("propertee_warn_on_no_unit_of_assignment" = TRUE)
}
