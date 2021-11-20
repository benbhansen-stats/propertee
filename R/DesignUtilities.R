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
unitid <- cluster

##' @rdname DesignSpecials
##' @export
block <- cluster

##' @rdname DesignSpecials
##' @export
forcing <- cluster

# Internal Function
# Perform checks on formula for creation of Design.
# Checks performed:
# - Ensure presence of exactly one of cluster() or unitid()
# - Disallow multiple cluster(), unitid(), block(), or forcing() terms
# - Disallow forcing() unless in RDD
.check_design_formula <- function(form, allowForcing = FALSE) {
  tt <- terms(form, c("cluster", "unitid", "block", "forcing"))
  specials <- attr(tt, "specials")

  if (attr(tt, "response") == 0) {
    stop("Must specify a treatment variable as the left side of the formula.")
  }

  if (is.null(specials$cluster) & is.null(specials$unitid)) {
    stop("Must specify a clustering or unitid variable.")
  } else if (length(specials$cluster) + length(specials$unitid) > 1) {
    # there's 2+ entered; need to figure out what combination

    # at least one of each
    if (length(specials$cluster) >= 1 & length(specials$unitid) >= 1) {
      stop(paste("Only one of cluster() or unitid() can be entered.",
                 "(`cluster()` and `unitid()` can take multiple variables)"))
    }

    # multiple `cluster()`, no `unitid()`
    if (length(specials$cluster) >= 1 & length(specials$unitid) == 0) {
      stop(paste("Only one cluster can be entered.",
                 "(`cluster()` can take multiple variables)"))
      }

    # multiple `unitid()`, no `cluster()`
    if (length(specials$cluster) == 0 & length(specials$unitid) >= 1) {
      stop(paste("Only one unitid can be entered.",
                 "(`unitid()` can take multiple variables)"))
    }
    stop("This error should never be seen, but something went wrong!")
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
