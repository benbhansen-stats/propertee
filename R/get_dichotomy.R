##' @title (Internal) Find and validate a \code{dichotomy}
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
  found_weights <- if (any(grepl("weights_calc", lapply(sys.calls(), "[[", 1),
                             fixed = TRUE))) FALSE else {
    vapply(lapply(sys.calls(), "[[", "weights"), function(aa) {
      is.call(aa) && (aa[[1]] == quote(ate) || aa[[1]] == quote(ett) ||
                        any(vapply(aa, function(x) {
                          x == quote(ate()) || x == quote(ett())
                        }, logical(1))))
    },   logical(1))
  }

  lmitt.dichotomy <- NULL
  weights.dichotomy <- NULL
  if (any(found_lmitt)) {
    lmitt.call <- get("lmitt.call", sys.frame(which(found_lmitt)[1]))
    lmitt.dichotomy <- lmitt.call$dichotomy
  }
  if (any(found_weights)) {
    weights.calls <- lapply(which(found_weights), function(pos) sys.call(pos)$weights)
    weights.dichotomy <- unique(lapply(weights.calls, function(cl) cl$dichotomy))
  }

  if (!is.null(dichotomy)) {
    if (!is.null(lmitt.dichotomy)) {
      if (!identical(deparse(lmitt.dichotomy), deparse(dichotomy))) {
        warning(paste("`dichotomy` passed to", paste0("`", sys.call(-1L)[[1L]], "()`"),
                      "is not the same as the `dichotomy` passed to `lmitt()`"))
      }
    }
    if (!is.null(weights.dichotomy)) {
      if (!identical(deparse(weights.dichotomy[[1]]), deparse(dichotomy))) {
        warning(paste("`dichotomy` passed to", paste0("`", sys.call(-1L)[[1L]], "()`"),
                      "is not the same as the `dichotomy` passed to",
                      paste0("`", unique(weights.calls)[[1L]][[1L]], "()`")))
      }
    }
  } else {
    dichotomy <- lmitt.dichotomy
    if (is.null(dichotomy)) dichotomy <- weights.dichotomy[[1]]
  }
  
  return(dichotomy)
}