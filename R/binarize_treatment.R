# (Internal) Uses a dichotomy formula to create a binary version of the
# treatment variable.
.binarize_treatment <- function(trt, dichot) {

  if (!is(dichot, "formula")) {
    stop("`dichotomy` must be formula")
  }

  if (!is.data.frame(trt)) {
    stop(paste("`trt` is expected to be a named `data.frame`",
               "(e.g. from `treatment(des)`)"))
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
    stop("At least one side for dichotomy formula must not be `.`")
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
      stop("treatment dichotomy overlaps")
    }
    return(ditxt)
  }

}
