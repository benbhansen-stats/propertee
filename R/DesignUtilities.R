##' These are special function used only in the definition of Design class
##' objects. They identify the clusters, blocks and forcing variables.
##'
##' @title Special terms in Design
##' @param ... any number of variables of the same length.
##' @return the variables with appropriate labels
##' @export
##' @rdname DesignSpecials
cluster <- function (...)
{
  #browser()
  allf <- list(...)
  do.call(cbind, allf)
}

##' @rdname DesignSpecials
##' @export
block <- cluster

##' @rdname DesignSpecials
##' @export
forcing <- cluster

# Perform checks on formula for creation of Design.
# Checks performed:
# - Ensure presence of cluster()
# - Disallow multiple cluster(), block(), or forcing() terms
# - Disallow forcing() unless in RDD
checkDesignFormula <- function(form, allowForcing = FALSE) {
  tt <- terms(form, c("cluster", "block", "forcing"))
  specials <- attr(tt, "specials")

  if (attr(tt, "response") == 0) {
    stop("Must specify a treatment variable as the left side of the formula.")
  }

  if (is.null(specials$cluster)) {
    stop("Must specify at least one clustering variable.")
  } else if (length(specials$cluster) > 1) {
    stop("Specify only one cluster() (cluster() can take multiple variables).")
  }

  if (!is.null(specials$block) && length(specials$block) > 1) {
    stop("Specify only one block() (block() can take multiple variables).")
  }

  if (!allowForcing && !is.null(attr(tt, "specials")$forcing)) {
    stop("forcing() only allowed in RD_Design")
  } else if (allowForcing && length(specials$forcing) > 1) {
    stop("Specify only one forcing() (forcing() can take multiple variables).")
  }

  TRUE
}
