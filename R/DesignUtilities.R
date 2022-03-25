##' These are special function used only in the definition of Design class
##' objects. They identify the units of assignment, blocks and forcing variables.
##'
##' `unit_of_assignment`, `uoa`, `cluster` and `unitid` are all aliases of each
##' other; you must include one and only one in each Design. The choice of which
##' to use will have no impact on any analysis.
##'
##' @title Special terms in Design
##' @param ... any number of variables of the same length.
##' @return the variables with appropriate labels
##' @export
##' @rdname DesignSpecials
unit_of_assignment <- function(...) {
  allf <- list(...)
  do.call(cbind, allf)
}

##' @rdname DesignSpecials
##' @export
unitid <- unit_of_assignment

##' @rdname DesignSpecials
##' @export
cluster <- unit_of_assignment

##' @rdname DesignSpecials
##' @export
uoa <- unit_of_assignment

##' @rdname DesignSpecials
##' @export
block <- unit_of_assignment

##' @rdname DesignSpecials
##' @export
forcing <- unit_of_assignment

# Internal Function
# Perform checks on formula for creation of Design.
# Checks performed:
# - Ensure presence of exactly one of unit_of_assignment(), cluster() or unitid()
# - Disallow multiple block(), or forcing() terms
# - Disallow forcing() unless in RDD
.check_design_formula <- function(form, allow_forcing = FALSE) {
  tt <- terms(form, c("unit_of_assignment", "uoa", "cluster", "unitid", "block", "forcing"))
  specials <- attr(tt, "specials")

  if (attr(tt, "response") == 0) {
    stop("Must specify a treatment variable as the left side of the formula.")
  }

  spec_uas <- specials$unit_of_assignment
  spec_uoa <- specials$uoa
  spec_clu <- specials$cluster
  spec_uni <- specials$unitid

  len_uas <- length(spec_uas)
  len_uoa <- length(spec_uoa)
  len_clu <- length(spec_clu)
  len_uni <- length(spec_uni)

  if (is.null(spec_uas) & is.null(spec_uoa) & is.null(spec_clu) & is.null(spec_uni)) {
    stop("Must specify a unit_of_assignment, cluster or unitid variable.")
  } else if (len_uas + len_uoa + len_clu + len_uni > 1) {
    # there's 2+ entered; need to figure out what combination

    if ((len_uas >= 1) + (len_uoa >= 1) + (len_clu >= 1) + (len_uni >= 1) > 1) {
      # There's more than one specified.
      stop(paste("Only one of `unit_of_assignment()`, `cluster()` or",
                 "`unitid()` can be entered."))
    } else {
      which <- switch(    (len_uoa >= 1) +
                          (len_uas >= 1) +
                      2 * (len_clu >= 1) +
                      3 * (len_uni >= 1),
                      "`unit_of_assignment()`",
                      "`cluster()`",
                      "`unitid()`")
      stop(paste0("Only one instance of ", which, " can be entered. (",
                 which, " can take multiple variables)"))
    }
  }

  if (!is.null(specials$block) && length(specials$block) > 1) {
    stop("Specify only one block() (block() can take multiple variables).")
  }

  if (!allow_forcing && !is.null(attr(tt, "specials")$forcing)) {
    stop("forcing() only allowed in rd_design")
  } else if (allow_forcing && length(specials$forcing) > 1) {
    stop("Specify only one forcing() (forcing() can take multiple variables).")
  }

  TRUE
}

# Internal function to rename cluster/unitid/uoa in a formula to
# unit_of_assignment
.update_form_to_unit_of_assignment <- function(form) {
    rename_list <- list("cluster" = as.name("unit_of_assignment"),
                        "uoa" = as.name("unit_of_assignment"),
                        "unitid" = as.name("unit_of_assignment"))
    as.formula(do.call("substitute", list(form, rename_list)))
}
