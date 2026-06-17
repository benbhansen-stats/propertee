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
  # weights arguments (no need to compare the dichotomy to itself)
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
    # possible_dichotomy is a list if it's passed the output of .find_dichotomies().
    # in that case, there's no need to compare it to other dichotomy arguments up
    # the stack, so we return a dichotomy if it exists in the list and NULL otherwise
    dichotomy <- possible_dichotomy$lmitt
    if (inherits(dichotomy, "name")) {
      # if the dichotomy passed to weights or lmitt was a stored object, find and
      # evaluate it
      for (s in seq_len(sys.nframe())) {
        eval_d <- tryCatch(eval.parent(dichotomy, s), error = function(e) NULL)
        if (!is.null(eval_d)) break
      }
      dichotomy <- eval_d
    }
    if (is.null(dichotomy)) {
      dichotomy <- possible_dichotomy$weights
    }
    if (inherits(dichotomy, "name")) {
      for (s in seq_len(sys.nframe())) {
        eval_d <- tryCatch(eval.parent(dichotomy, s), error = function(e) NULL)
        if (!is.null(eval_d)) break
      }
      dichotomy <- eval_d
    }
    if (!is.null(dichotomy)) dichotomy <- as.formula(dichotomy)
  } else if (inherits(possible_dichotomy, "call")) {
    possible_dichotomy <- as.formula(possible_dichotomy)
  } else if (inherits(possible_dichotomy, "name")) {
    # if the dichotomy passed to weights or lmitt was a stored object, find and
    # evaluate it
    for (s in seq_len(nf <- sys.nframe())) {
      eval_pd <- tryCatch(
        eval.parent(possible_dichotomy, s),
        error = function(e) if (s == nf) {
          stop(paste("Could not find", deparse1(possible_dichotomy), "in call stack"),
               call. = FALSE)
        } else return(NULL)
      )
      if (!is.null(eval_pd)) break
    }
    possible_dichotomy <- eval_pd
    if (!is.null(possible_dichotomy)) possible_dichotomy <- as.formula(possible_dichotomy)
  }

  if (inherits(possible_dichotomy, "formula")) {
    dichotomy <- possible_dichotomy
    # find any other dichotomies up the call stack to compare dichotomy against
    other_dichotomies <- .find_dichotomies()
    if (!is.null(other_dichotomies$lmitt)) {
      if (inherits(other_dichotomies$lmitt, "name")) {
        for (s in seq_len(sys.nframe())) {
          ok <- tryCatch(eval.parent(other_dichotomies$lmitt, s),
                         error = function(e) NULL)
          if (!is.null(ok)) break
        }
        if (is.null(ok)) stop(paste("Could not find", deparse1(other_dichotomies$lmitt), "in call stack"))
        other_dichotomies$lmitt <- ok
      }
      if (!identical(deparse1(other_dichotomies$lmitt),
                     do.call("deparse1", list(eval(dichotomy))))) {
        warning(paste("`dichotomy` passed to", paste0("`", sys.call(-1L)[[1L]], "()`"),
                      "is not the same as the `dichotomy` passed to `lmitt()`"))
      }
    }
    if (!is.null(other_dichotomies$weights)) {
      if (inherits(other_dichotomies$weights, "name")) {
        for (s in seq_len(sys.nframe())) {
          ok <- tryCatch(eval.parent(other_dichotomies$weights, s),
                         error = function(e) NULL)
          if (!is.null(ok)) break
        }
        if (is.null(ok)) stop(paste("Could not find", deparse1(other_dichotomies$weights), "in call stack"))
        other_dichotomies$weights <- ok
      }
      if (!identical(deparse1(other_dichotomies$weights),
                     do.call("deparse1", list(eval(dichotomy))))) {
        warning(paste("`dichotomy` passed to", paste0("`", sys.call(-1L)[[1L]], "()`"),
                      "is not the same as the `dichotomy` passed to weights"))
      }
    }
  } else if (!exists("dichotomy")) {
    stop("No formulas provided to .validate_dichotomy")
  }

  return(invisible(dichotomy))
}
