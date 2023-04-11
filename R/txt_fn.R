##' (Internal) Obtain which variaiton of \code{assigned()}, \code{a.()} or
##' \code{z.()} is used in the model
##' @title Treatment variable name
##' @param object \code{DirectAdjusted} object
##' @return Character string identifying the name (e.g "assigned()" or
##'   "a.(des)")
##' @keywords internal
.txt_fn <- function(object) {
  cs <- names(object$coefficients)
  found_assigned <- FALSE
  if (any(grepl("assigned\\(.*\\)", cs))) {
    found_assigned <- TRUE
    assigned_form <- regmatches(cs, regexpr("assigned\\(.*\\)", cs))
    if (length(unique(assigned_form)) > 1) {
      stop("Differing forms of `assigned()` found. Keep form consistent.")
    }
    assigned_form <- assigned_form[1]
  }
  found_a. <- FALSE
  if (any(grepl("a\\.\\(.*\\)", cs))) {
    found_a. <- TRUE
    a._form <- regmatches(cs, regexpr("a\\.\\(.*\\)", cs))
    if (length(unique(a._form)) > 1) {
      stop("Differing forms of `a.()` found. Keep form consistent.")
    }
    a._form <- a._form[1]
  }
  found_z. <- FALSE
  if (any(grepl("z\\.\\(.*\\)", cs))) {
    found_z. <- TRUE
    z._form <- regmatches(cs, regexpr("z\\.\\(.*\\)", cs))
    if (length(unique(z._form)) > 1) {
      stop("Differing forms of `z.()` found. Keep form consistent.")
    }
    z._form <- z._form[1]

  }

  if (sum(found_assigned, found_a., found_z.) > 1) {
    stop(paste("Differing treatment identification (`assigned()`, `a.()`",
               " or `z.()`) found. Only one can be used."))
  }
  if (sum(found_assigned, found_a., found_z.) == 0) {
    stop("No treatment variable found")
  }

  if (found_assigned) {
    return(assigned_form)
  }
  if (found_a.) {
    return(a._form)
  }
  if (found_z.) {
    return(z._form)
  }
  stop("This error should never be hit!")
  return(NULL)


}
