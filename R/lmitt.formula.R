#' Formulas generated from inputs to \code{lmitt()}
#' 
#' The \code{terms()} method for this class returns a \code{lmitt.terms}
#' object that helps generate artifacts for fitting \code{teeMod} objects
#' @keywords internal
setClass("lmitt.formula",
         contains = "formula")

#' Build a formula from inputs to \code{lmitt()}
#' 
#' @description
#' The formula passed to \code{lmitt()} does not always contain the extent of
#' information related to treatment assignment, nor does it express the user's
#' preference to absorb block fixed effects. \code{build_lmitt_formula()}
#' synthesizes this information into a single formula. Users need not call it
#' directly, as it is called within \code{lmitt()}. Note this is not the
#' formula ultimately used in model fitting; see \code{terms.lmitt.formula()}.
#' 
#' @details
#' The righthand side of the formula always has two terms: the first informs
#' construction of the treatment assignment variable, and the second informs
#' construction of the intercept and optional moderator variable. The
#' \code{rhstype} argument determines whether the fitted model has terms just
#' for treatment assignment and an intercept (\code{"intercept"}); an intercept,
#' columns for each level of the moderator, and interactions between treatment
#' assignment and the moderator (\code{"categorical"}); or an intercept,
#' treatment assignment, moderator, and interactions between the two
#' (\code{"continuous"}). The \code{absorb} argument dictates whether these
#' terms will be centered within blocks.
#' 
#' The \code{specification} and \code{data} arguments are left NULL when called
#' within the S3 method \code{lmitt.formula()}, but in \code{lmitt.lm()}, they
#' are set based on the \code{lm} object's call; the treatment assignment
#' variable would otherwise fail to be constructed in model artifact generation.
#' @returns A \code{lmitt.formula}.
#' @param lhs character, name of the response variable. 
#' @param rhs character, name of the moderator variable. May be "1"
#' if no moderator variable was provided in the original formula.
#' @param rhstype character, one of \code{"intercept"}, \code{"categorical"},
#' or \code{"continuous"}. 
#' @param absorb logical, \code{TRUE} if moderator variable is to be block-
#' centered.
#' @param dichotomy formula, specifies the dichotomization of a non-binary
#' treatment assignment variable. Defaults to NULL.
#' @param contrasts quoted function, specifies the contrast function for the
#' moderator variable. Defaults to NULL.
#' @param specification name, name of the specification object to look for
#' in the call stack. Defaults to NULL (see Details).
#' @param data name, name of the data to use to form the assignment column.
#' Defaults to NULL (see Details).
#' @export
#' @rdname build_lmitt_formula
#' @examples
#' build_lmitt_formula("y", "1", "intercept", FALSE)
#' build_lmitt_formula("y", "mod", "categorical", TRUE)
build_lmitt_formula <- function(lhs,
                                rhs,
                                rhstype,
                                absorb,
                                dichotomy = NULL,
                                contrasts = NULL,
                                specification = NULL,
                                data = NULL) {
  form <- stats::reformulate(
    vapply(
      list(call("a.",
                dichotomy = dichotomy,
                specification = specification,
                data = data),
           call("moderator",
                str2lang(rhs),
                rhstype = rhstype,
                absorb = absorb,
                contrasts = contrasts)),
      deparse1,
      character(1L)
    )
    , response = lhs
    , env = parent.frame(2L)
  )
  return(as(form, "lmitt.formula"))
}