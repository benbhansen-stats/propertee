#' @include Design.R DesignAccessors.R
NULL

##' These are special functions used only in the definition of \code{Design}
##' objects. They identify the units of assignment, blocks and forcing
##' variables.
##'
##' These functions have no use outside of the formula in creating a
##' \code{Design}.
##'
##' \code{unit_of_assignment}, \code{uoa}, \code{cluster} and \code{unitid} are
##' all aliases of each other; you must include one and only one in each
##' \code{Design}. The choice of which to use will have no impact on any
##' analysis, only on some output display.
##'
##' See \code{rct_design}, \code{obs_design}, or \code{rd_design} for examples
##' of their usage.
##'
##' @title Special terms in \code{Design} creation formula
##' @param ... any number of variables of the same length.
##' @return the variables with appropriate labels
##' @export
##' @rdname DesignSpecials
unit_of_assignment <- function(...) {
  allf <- list(...)
  return(do.call(cbind, allf))
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

#
##' Checks performed:
##' * Ensure presence of exactly one of \code{unit_of_assignment()},
##'   \code{cluster()} or \code{unitid()}.
##' * Disallow multiple \code{block()} or multiple \code{forcing()} terms.
##' * Disallow \code{forcing()} unless in RDD.
##'
##' @title (Internal) Perform checks on formula for creation of Design.
##' @param form A formula passed to \code{.new_Design()}
##' @param allow_forcing Binary whether \code{forcing()} is allowed (\code{TRUE}
##'   for RDD, \code{FALSE} for RCT and Obs).
##' @return \code{TRUE} if all checks pass, otherwise errors.
##' @keywords internal
.check_design_formula <- function(form, allow_forcing = FALSE) {
  tt <- terms(form, c("unit_of_assignment", "uoa", "cluster",
                      "unitid", "block", "forcing"))
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

  if (is.null(spec_uas) &
      is.null(spec_uoa) &
      is.null(spec_clu) &
      is.null(spec_uni)) {
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

  invisible(TRUE)
}


##' Internally, we always refer to uoa/cluster/unitid as "unit_of_assignment"
##' @title (Internal) Rename cluster/unitid/uoa in a formula to
##'   unit_of_assignment for internal consistency
##' @param form A formula passed to \code{.new_Design()}
##' @return The formula with "cluster"/"unitid"/"uoa" replace with
##'   "unit_of_assignment"
##' @keywords internal
.update_form_to_unit_of_assignment <- function(form) {
  form <- deparse(form)
  form <- gsub("cluster\\(", "unit_of_assignment(", form)
  form <- gsub("uoa\\(", "unit_of_assignment(", form)
  form <- gsub("unitid\\(", "unit_of_assignment(", form)
  return(as.formula(form))
}

##' Check if \code{Design} has access to a binary treatment
##'
##' These functions determine if a \code{Design} is dichotomized
##' (\code{is_dichotomized()}) or has access to a binary treatment variable
##' (\code{has_binary_treatment()}) or has either
##' (\code{is_binary_or_dichotomized()}).
##'
##' \code{is_dichotomized()} checks for the the presence of a \code{@dichotomy}
##' slot in the \code{Design}, either passed directly into the various
##' \code{*_design()} \code{Design} creators, or added afterwards with
##' \code{dichotomy(my_design)<-}.
##'
##' \code{has_binary_treatment()} returns \code{TRUE} if the treatment is
##' numeric with only values of \code{0} and \code{1}, or is logical. It also
##' allows \code{NA} in the treatment.
##'
##' \code{is_binary_or_dichotomized} returns \code{TRUE} if either
##' \code{is_dichotomized()} and \code{has_binary_treatment()} return
##' \code{TRUE}.
##'
##' @param des a \code{Design} object
##' @return Logical; see details.
##' @export
##' @rdname design_treatment_status
##' @examples
##' data(simdata)
##' des1 <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
##' des2 <- rct_design(o ~ uoa(cid1, cid2), data = simdata)
##' des3 <- rct_design(o ~ uoa(cid1, cid2), data = simdata,
##'                    dichotomy = o >= 3 ~ .)
##'
##' if (is_binary_or_dichotomized(des1)) {
##'   wt <- ate(des1, data = simdata)
##' }
##' if (is_binary_or_dichotomized(des2)) {
##'   wt <- ate(des2, data = simdata) # would fail but won't be run
##' }
##' if (is_binary_or_dichotomized(des3)) {
##'   wt <- ate(des3, data = simdata)
##' }
##'
is_dichotomized <- function(des) {
  if (!inherits(des, "Design")) {
    stop("des must be a Design object.")
  }

  # legnth 3 formula implies a LHS, a ~, and a RHS.
  return(length(des@dichotomy) == 3)
}

##' @export
##' @rdname design_treatment_status
has_binary_treatment <- function(des) {
  if (!inherits(des, "Design")) {
    stop("des must be a Design object.")
  }

  if (is.factor(treatment(des)[, 1])) {
    # Short circuit if treatment is a factor
    return(FALSE)
  }

  return(all(treatment(des)[, 1] %in% c(0, 1, TRUE, FALSE, NA)))
}

##' @export
##' @rdname design_treatment_status
is_binary_or_dichotomized <- function(des) {
  if (!inherits(des, "Design")) {
    stop("des must be a Design object.")
  }

  return(is_dichotomized(des) || has_binary_treatment(des))
}

##' Check if \code{Design} objects are identical
##' @param x A \code{Design} object.
##' @param y A \code{Design} object.
##' @param dichotomy_force Logical, default \code{FALSE}. If \code{FALSE}, the
##'   \code{dichotomy()} of \code{x} and \code{y} is ignored in the comparison.
##'   (In other words, two \code{Design}s which differ only in their
##'   \code{@dichotomy} slot will be considered identical.) If \code{TRUE}, the
##'   \code{dichotomy} must also be in agreement between \code{x} and \code{y}.
##' @return Logical, are \code{x} and \code{y} identical?
##' @export
identical_Designs <- function(x, y, dichotomy_force = FALSE) {
  if (!dichotomy_force) {
        x@dichotomy <- stats::formula(env = globalenv())
        x@call$dichotomy <- NULL
        y@dichotomy <- stats::formula(env = globalenv())
        y@call$dichotomy <- NULL
  }
  return(identical(x, y))
}

##' @title Identify "small" blocks (one treated or one control unit of assignment)
##' @param des A \code{Design} object.
##' @return Logical vector with length given by the number of blocks in the
##' Design
##' @export
identify_small_blocks <- function(des) {
  blk_txt_cts <- design_table(des, "t", "b")
  is_small_blk <- apply(blk_txt_cts, 1, function(blk) any(blk == 1))
  return(is_small_blk)
}

##' @title Make a dataframe that links units of assignment with clusters
##' @param des A \code{Design} object.
##' @param cluster A character vector of column names to use as clusters. Columns
##' must exist in the dataframe used to create the \code{Design} object. Defaults
##' to NULL, in which case the column names specified in the \code{unitid()},
##' \code{unit_of_assignment()}, or \code{cluster()} function in the \code{Design}
##' formula will be used.
##' @return A dataframe where the number of rows coincides with the number of
##' distinct unit of assignment or cluster combinations (depending on whether
##' `cluster` is a more or less granular level than the assignment level) and
##' the columns correspond to the unit of assignment columns and a "cluster"
##' column
##' @export
make_uoa_cluster_df <- function(des, cluster = NULL) {
  if (!inherits(des, "Design")) stop("Must be provided a valid `Design` object")
  uoa_cols <- var_names(des, "u")
  q_df <- NULL
  des_cl <- des@call
  if (is.null(subset_cl <- des_cl$subset)) desdata_cl <- des_cl$data else {
    desdata_cl <- quote(subset(x = df, subset_arg))
    desdata_cl$x <- des_cl$data
    desdata_cl[[3]] <- subset_cl
  }
  for (f in seq_len(sys.nframe())) {
    q_df <- tryCatch({
      eval(desdata_cl, envir = parent.frame(f))
    }, error = function(e) return(NULL))
    if (!is.null(q_df) & inherits(q_df, "data.frame")) break
  }
  
  if (is.null(q_df)) {
    stop("Could not find design data in the call stack")
  }
  
  if (!is.null(cluster) & !all(cluster %in% colnames(q_df))) {
    stop(paste("Could not find", cluster, "column in the design data"))
  }
  
  q_df <- q_df[, c(uoa_cols, cluster), drop = FALSE]
  grab_uoas_fn <- switch(
    des@unit_of_assignment_type,
    "unitid" = unitids,
    "unit_of_assignment" = units_of_assignment,
    "cluster" = clusters
  )
  uoas <- grab_uoas_fn(des)
  
  out <- unique(merge(uoas, q_df, by = uoa_cols, all.y = TRUE))
  rownames(out) <- NULL
  if (nrow(out) < nrow(uoas)) warning(paste("Some units of assignment in the Design",
                                            "were not found in the data used to",
                                            "create the Design object. Ensure",
                                            "the original data has not been",
                                            "modified."))

  if (is.null(cluster_cols <- cluster)) cluster_cols <- uoa_cols
  out$cluster <- apply(
    out[, cluster_cols, drop = FALSE],
    1,
    function(...) paste(..., collapse = "_")
  )

  return(out[, c(uoa_cols, "cluster")])
}
