##' Add dichotomization to a design
##'
##' Generally users should leave `force = TRUE`. However, if you want to avoid
##' accidentally changing the dichotomization, set `force = FALSE`.
##' @param des a Design object
##' @param dichotomize optionally, a formula defining the dichotomization of the
##'   treatment variable if it isn't already \code{0}/\code{1}. See details of
##'   help for \code{\link{rct_design()}} e.g. for details.
##' @param force If `des` is already dichotomized, should we force a new
##'   dichomization (TRUE, default) or should we instead return the original
##'   design (FALSE)?
##' @return A dichotomized design
##' @export
dichotomize_design <- function(des, dichotomize, force = TRUE) {
  if (!is(des, "Design")) {
    stop("`des` must be a Design created by *_design() functions")
  }
  if (!is(dichotomize, "formula")) {
    stop("`dichotomize` must be a formula")
  }
  if (is_dichotomized(des) & !force) {
    warning("Design is already dichotomized, no further action taken")
    return(des)
  }

  if (is_dichotomized(des) & force) {
    txt <- des@treatment_binary[[2]]
  } else {
    txt <- treatment(des)
    des@treatment_binary <- list(dichotomize = dichotomize,
                                 original_treatment = txt)
  }

  bin_txt <- .binarize_treatment(txt, dichotomize)
  treatment(des) <- data.frame(z__ = bin_txt)

  # Update original design's call for replicability
  orig_call <- des@call
  orig_call$dichotomize <- dichotomize
  des@call <- orig_call

  return(des)
}

# (Internal) Uses a dichotomization formula to create a binary version of the
# treatment variable.
.binarize_treatment <- function(trt, dichot) {

  if (!is(dichot, "formula")) {
    stop("`dichotomization` must be formula")
  }

  lhs_dot <- rhs_dot <- FALSE
  if (dichot[[3]] == ".") {
    # control group is .
    # control goes first since txt switches LHS and RHS
    dichot[[3]] <- 1
    rhs_dot <- TRUE
  }
  if (dichot[[2]] == ".") {
    # treatment group is .
    dichot[[2]] <- dichot[[3]]
    dichot[[3]] <- 1
    lhs_dot <- TRUE
  }
  if (lhs_dot & rhs_dot) {
    stop("Both dots")
  }

  m <- model.frame(dichot, trt)

  if (lhs_dot) {
    return(as.numeric(!m[,1]))
  } else if (rhs_dot) {
    return(as.numeric(m[,1]))
  } else {
    ditxt <- m[,2] + 2*m[,1] - 1
    ditxt[ditxt == -1] <- NA
    if (any(!is.na(ditxt) & ditxt > 1)) {
      stop("treatment dichotoization overlaps")
    }
    return(ditxt)
  }

}
