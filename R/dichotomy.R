##' @title (Internal) Find \code{dichotomy} formulas in the call stack
##' @details \code{.find_dichotomies()} searches for \code{lmitt.formula()}
##'   calls and their \code{weights} arguments for any \code{dichotomy}
##'   arguments.
##' @return A list where elements are either formulas or NULL, depending on
##'   whether a \code{dichotomy} argument was found in \code{lmitt.formula()} or
##'   its \code{weights} argument
##' @keywords internal
.find_dichotomies <- function() {
  found_lmitt <- grepl("^lmitt\\.formula$",
                       lapply(sys.calls(), "[[", 1),
                       perl = TRUE)
  # if this is called within `.weights_calc()` don't look for dichotomies in the
  # weights arguments (it already would have been available)
  weights.calls <- if (any(grepl("weights_calc", lapply(sys.calls(), "[[", 1),
                                 fixed = TRUE))) FALSE else {
    lapply(sys.calls(), "[[", "weights")
  }

  dichotomies <- list(lmitt = NULL, weights = NULL)
  if (any(found_lmitt)) {
    lmitt.call <- get("lmitt.call", sys.frame(which(found_lmitt)[1L]))
    dichotomies$lmitt <- lmitt.call$dichotomy
  }
  if (length(valid_pos <- which(!vapply(weights.calls, is.null, logical(1L)))) > 0) {
    dichotomies$weights <- mapply(
      function(cc, pos) {
        if (is.call(cc)) cc$dichotomy else if (
          inherits(wd <- eval("weights", pos), "WeightedStudySpecification")) wd@dichotomy else NULL
        },
      weights.calls[valid_pos],
      valid_pos,
      SIMPLIFY = FALSE)[[1L]]
  }

  return(dichotomies)
}

#' @title (Internal) Validate a dichotomy against other dichotomies found in the
#'   call stack
#' @param possible_dichotomy a list or a formula, the former coming from
#'   \code{.find_dichotomies()}
#' @return formula for a dichotomy or NULL
#' @keywords internal
.validate_dichotomy <- function(possible_dichotomy) {
  if (inherits(possible_dichotomy, "list")) {
    dichotomy <- possible_dichotomy$lmitt
    if (is.null(dichotomy)) {
      dichotomy <- possible_dichotomy$weights
    }
  } else if (inherits(possible_dichotomy, "call")) {
    possible_dichotomy <- as.formula(possible_dichotomy)
  }
  
  if (inherits(possible_dichotomy, "formula")) {
    dichotomy <- possible_dichotomy
    other_dichotomies <- .find_dichotomies()
    if (!is.null(other_dichotomies$lmitt)) {
      if (!identical(deparse1(other_dichotomies$lmitt), deparse1(dichotomy))) {
        warning(paste("`dichotomy` passed to", paste0("`", sys.call(-1L)[[1L]], "()`"),
                      "is not the same as the `dichotomy` passed to `lmitt()`"))
      }
    }
    if (!is.null(other_dichotomies$weights)) {
      if (!identical(deparse1(other_dichotomies$weights), deparse1(dichotomy))) {
        warning(paste("`dichotomy` passed to", paste0("`", sys.call(-1L)[[1L]], "()`"),
                      "is not the same as the `dichotomy` passed to weights"))
      }
    }
  } else if (!exists("dichotomy")) {
    stop("No formulas provided to .validate_dichotomy")
  }

  return(invisible(dichotomy))
}
