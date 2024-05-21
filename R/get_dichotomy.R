##' @title Find and validate a \code{dichotomy}
##' @param dichotomy a formula or call evaluating to a formula
##' @details Whether a \code{dichotomy} is provided or not, \code{.get_dichotomy()}
##' finds any \code{dichotomy} argument passed to an \code{lmitt.formula()} call
##' up the stack. This argument is either used itself or used to validate a
##' provided \code{dichotomy}. 
##' @return A formula or a call that evaluates to a formula, or NULL if no
##' \code{dichotomy} is provided nor found in a \code{lmitt.formula()} call up
##' stack
##' @keywords internal
.get_dichotomy <- function(dichotomy = NULL) {
  found_lmitt <- grepl("^lmitt\\.formula$",
                       lapply(sys.calls(), "[[", 1),
                       perl = TRUE)
  lmitt.dichotomy <- NULL
  if (any(found_lmitt)) {
    lmitt.call <- get("lmitt.call", sys.frame(which(found_lmitt)[1]))
    lmitt.dichotomy <- lmitt.call$dichotomy
  }
  
  if (!is.null(dichotomy)) {
    if (!inherits(dichotomy, "formula")) {
      stop("`dichotomy` must be a `formula`")
    }
    if (!is.null(lmitt.dichotomy)) {
      if (!all.equal(lmitt.dichotomy, dichotomy)) {
        warning(paste("`dichotomy` passed to", paste0("`", sys.call(1L)[[1L]], "()`"),
                      "is not the same as the `dichotomy` passed to `lmitt()`"))
      }
    }
  } else {
    if (!is.null(dichotomy)) dichotomy <- lmitt.dichotomy
  }
  
  return(dichotomy)
}