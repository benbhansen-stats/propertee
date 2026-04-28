#' lmitt.formula class

#' Class for lmitt formulas
#'
#' @details
#'  This has a specific \code{model.frame()} method so that \code{lmitt()} fits
#'  the model we desire
setClass("lmitt.formula",
         contains = "formula")

#' Build the formula for \code{lmitt()}
#' @returns Returns a formula that has a term for the treatment assignment
#' variable and, if provided in \code{lmitt()}, a moderator term that will be
#' appropriately block-centered (as indicated by the \code{absorb} argument of
#' \code{lmitt()}) and interacted with the treatment assignment variable (based
#' on whether the moderator variable is continuous or categorical)
#' @param moderator_var character, name of the moderator variable. May be "1"
#' if no moderator variable was provided in the original formula.
#' @param moderator_type character, one of \code{"intercept"},
#' \code{"continuous"}, or \code{"categorical"}. 
#' @param absorb logical, \code{TRUE} if moderator variable is to be block-
#' centered.
#' @param dichotomy formula, specifies the dichotomization of a non-binary
#' treatment assignment variable. Defaults to NULL if the treatment assignment
#' variable is already binary.
#' @export
build_lmitt_formula <- function(moderator_var,
                                moderator_type,
                                absorb,
                                dichotomy = NULL) {
  if (!is.null(dichotomy)) dichotomy <- deparse1(dichotomy)
  form <- stats::reformulate(
    vapply(
      list(call("a.", dichotomy = dichotomy),
           call("areg.center",
                str2lang(moderator_var),
                rhstype = moderator_type,
                absorb = absorb)),
      deparse,
      character(1L)
    )
  )
  return(as(form, "lmitt.formula"))
}