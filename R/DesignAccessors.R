############### Treatment

##' @export
##' @rdname Design_extractreplace
setGeneric("treatment", function(x, binary = FALSE, newdata = NULL, by = NULL, ...) {
  standardGeneric("treatment")
})

##' Allows access to the elements which define a \code{Design}, enabling their
##' extraction or replacement.
##'
##' For \code{treatment()}, when argument \code{binary} is \code{FALSE}, the
##' treatment variable passed into the \code{Design} is returned as a one-column
##' \code{data.frame}. If \code{binary = TRUE} is passed, and the \code{Design}
##' either has a binary treatment variable, or has a \code{dichotomy}, a binary
##' one-column \code{data.frame} will be returned. If the \code{Design} does not
##' have access to binary treatment (non-binary treatment and no
##' \code{dichotomy} specified), passing \code{binary = TRUE} will error.
##'
##' \code{binary = "ifany"} is the most permissible; returning the dichotomized
##' treatment variable if \code{@dichotomy} exists, otherwise returning the
##' original treatment without error.
##'
##' The one-column \code{data.frame} returned by \code{treatment()} is named as
##' entered in the \code{Design} creation, but if a \code{dichotomy} is in the
##' \code{Design}, the column name is \code{"__z"}.
##' @title Accessors and Replacers for \code{Design} objects
##' @param x \code{Design} object
##' @param binary If \code{FALSE} (default), return a \code{data.frame}
##'   containing the named treatment variable. If \code{TRUE} and \code{x} has a
##'   formula in \code{@dichotomy}, return a \code{data.frame} containing a
##'   binary treatment variable with the name \code{"z__"}. Errors on
##'   \code{TRUE} if treatment is non-binary \code{@dichotomy} is \code{NULL}.
##'   If \code{"ifany"}, returns a binary treatment if possible (if treatment is
##'   already binary, or there's a valid \code{@dichotomy}), otherwise return
##'   original treatment.
##' @param newdata Optionally an additional \code{data.frame}. If passed, and
##'   the unit of assignment variable is found in \code{newdata}, then the
##'   requested variable type for each unit of \code{newdata} is returned. See
##'   \code{by} argument if the name of the unit of assignment differs.
##' @param by optional; named vector or list connecting names of cluster/unit of
##'   assignment variables in \code{design} to cluster/unit of assignment
##'   variables in \code{data}. Names represent variables in the Design; values
##'   represent variables in the data. Only needed if variable names differ.
##' @param ... Ignored.
##' @return \code{data.frame} containing treatment variable
##' @export
##' @rdname Design_extractreplace
setMethod("treatment", "Design", function(x, binary = FALSE, newdata = NULL, by = NULL, ...) {
  binary <- as.character(binary)
  if (!binary %in% c("TRUE", "FALSE", "ifany")) {
    stop(paste("Valid input to `binary=` argument include only TRUE, ",
               "FALSE, and 'ifany'."))
  }
  .design_accessors_newdata_validate(newdata, by)

  # Case 1: binary = FALSE
  if (binary == FALSE) {
    # Return original treatment
    if (!is.null(newdata)) {
      return(.get_col_from_new_data(x, newdata, type = "t", by))
    }
    return(x@structure[x@column_index == "t"])
  }

  # Case 2: binary = TRUE
  if (binary == TRUE) {

    # Case 2a: binary = TRUE, treatment is stored as binary
    if (has_binary_treatment(x)) {
      # Treatment is binary, return original treatment
      if (!is.null(newdata)) {
        return(.get_col_from_new_data(x, newdata, type = "t", by))
      }
      return(x@structure[x@column_index == "t"])
    }

    # Case 2b: binary = TRUE, stored treatment is non-binary but has dichotomy
    if (is_dichotomized(x)) {
      # Has dichotomization, return that
      if (!is.null(newdata)) {
        treatment(x) <- .bin_txt(x)
        return(.get_col_from_new_data(x, newdata, type = "t", by))
      }
      return(data.frame(z__ = .bin_txt(x)))
    }

    # Case 2c: binary = TRUE, treatment is non-binary and no dichotomy
    stop(paste("No binary treatment can be produced. Treatment is",
               "non-binary and `x` does not contain a `@dichotomy`."))
  }

  # Case 3: binary = "ifany"
  if (binary == "ifany") {

    # Case 3a: binary = "ifany", no dichotomy
    if (!is_dichotomized(x)) {
      # Return original treatment
      if (!is.null(newdata)) {
        return(.get_col_from_new_data(x, newdata, type = "t", by))
      }
      return(x@structure[x@column_index == "t"])
    }

    # Case 3b: binary = "ifany", dichotomy
    if (is_dichotomized(x)) {
      # Has dichotomization, return that
      if (!is.null(newdata)) {
        treatment(x) <- .bin_txt(x)
        return(.get_col_from_new_data(x, newdata, type = "t", by))
      }
      return(data.frame(z__ = .bin_txt(x)))
    }
  }

})

##' @export
##' @rdname Design_extractreplace
setGeneric("treatment<-", function(x, value) standardGeneric("treatment<-"))

##' @param value Replacement. Either a \code{vector}/\code{matrix} of
##'   appropriate dimension, or a named \code{data.frame} if renaming variable
##'   as well.
##' @export
##' @rdname Design_extractreplace
setMethod("treatment<-", "Design", function(x, value) {
  value <- .convert_to_data.frame(value, x, "t")

  x@structure[x@column_index == "t"] <- value
  names(x@structure)[x@column_index == "t"] <- colnames(value)
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

##' If the \code{Design} has a \code{@dichotomy}, or has a treatment variable
##' consisting only of 0/1 or \code{NA}, then returns the binary treatment.
##' Otherwise (it has a non-binary treatment and lacks a dichotomy) it errors.
##' @title (Internal) Extracts treatment as binary \code{vector} if possible or
##'   else errors.
##' @param des A \code{Design}
##' @return A \code{vector} of binary treatments
##' @keywords internal
.bin_txt <- function(des) {
  if (!is_dichotomized(des)) {
    tt <- treatment(des, binary = FALSE)[, 1]
    if (!all(tt %in% c(0:1, NA))) {
      stop("binary treatment cannot be obtained")
    }
    return(tt)
  }
  return(.apply_dichotomy(treatment(des, binary = FALSE),
                          des@dichotomy))
}

##' Given a treatment variable (passed as a named \code{data.frame}) and a
##' dichotomy formula (see help on \code{rct_design()} for details on
##' specification), returns \code{vector} containing only \code{0}, \code{1}, or
##' \code{NA}.
##' @title (Internal) Applies dichotomy to treatment
##' @param txt A named \code{data.frame} containing a single column of the
##'   treatment, such as that produed by `treatment(mydesign)`.
##' @param dichotomy A dichotomization formula. See the details in the Details
##'   for the help of \code{rct_design()}.
##' @return A \code{vector} of binary treatments
##' @keywords internal
.apply_dichotomy <- function(txt, dichotomy) {

  if (!inherits(dichotomy, "formula")) {
    stop("`dichotomy` must be formula")
  }

  if (!is.data.frame(txt)) {
    stop(paste("`txt` is expected to be a named `data.frame`",
               "(e.g. from `treatment(des)`)"))
  }

  lhs_dot <- rhs_dot <- FALSE
  if (dichotomy[[3]] == ".") {
    # control group is .
    # control goes first since txt switches LHS and RHS
    dichotomy[[3]] <- 1
    rhs_dot <- TRUE
  }
  if (dichotomy[[2]] == ".") {
    # treatment group is .
    dichotomy[[2]] <- dichotomy[[3]]
    dichotomy[[3]] <- 1
    lhs_dot <- TRUE
  }
  if (lhs_dot & rhs_dot) {
    stop("At least one side for dichotomy formula must not be `.`")
  }

  m <- model.frame(dichotomy, txt, na.action = na.pass)

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

############### Units of Assignment

##' @export
##' @rdname Design_extractreplace
setGeneric("units_of_assignment", function(x, newdata = NULL, by = NULL) {
  standardGeneric("units_of_assignment")
})

##' @export
##' @rdname Design_extractreplace
setMethod("units_of_assignment", "Design", function(x, newdata = NULL, by = NULL) {
  if (x@unit_of_assignment_type == "unitid") {
    stop("Design specified with `unitid()`, not `unit_of_assignment()`")
  }
  if (x@unit_of_assignment_type == "cluster") {
    stop("Design specified with `cluster()`, not `unit_of_assignment()`")
  }
  .design_accessors_newdata_validate(newdata, by)
  if (!is.null(newdata)) {
    return(.get_col_from_new_data(x, newdata, type = "u", by))
  }
  return(x@structure[x@column_index == "u"])
})

##' @export
##' @rdname Design_extractreplace
setGeneric("units_of_assignment<-", function(x, value) {
  standardGeneric("units_of_assignment<-")
})

##' @export
##' @rdname Design_extractreplace
setMethod("units_of_assignment<-", "Design", function(x, value) {

  value <- .convert_to_data.frame(value, x, "u")

  x <- .update_structure(x, value, "u")

  dupuoa <- duplicated(units_of_assignment(x))
  dupall <- duplicated(x@structure[x@column_index != "f"])
  if (any(dupuoa)) {

    if (sum(dupuoa) != sum(dupall)) {
      stop(paste("Fewer new units of assignment then original, but new",
                 "collapsed units would have non-constant treatment and/or",
                 "block structure"))
    }
    warning("Fewer new units of assignment then original, collapsing")

    x@structure <- x@structure[-dupuoa, ]
  }
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

############### Cluster

##' @export
##' @rdname Design_extractreplace
setGeneric("clusters", function(x, newdata = NULL, by = NULL) standardGeneric("clusters"))

##' @export
##' @rdname Design_extractreplace
setMethod("clusters", "Design", function(x, newdata = NULL, by = NULL) {
  if (x@unit_of_assignment_type == "unitid") {
    stop("Design specified with `unitid()`, not `cluster()`")
  }
  if (x@unit_of_assignment_type == "unit_of_assignment") {
    stop("Design specified with `unit_of_assignment()`, not `cluster()`")
  }
  .design_accessors_newdata_validate(newdata, by)
  if (!is.null(newdata)) {
    return(.get_col_from_new_data(x, newdata, type = "u", by))
  }
  return(x@structure[x@column_index == "u"])
})

##' @export
##' @rdname Design_extractreplace
setGeneric("clusters<-", function(x, value) standardGeneric("clusters<-"))

##' @export
##' @rdname Design_extractreplace
setMethod("clusters<-", "Design", function(x, value) {

  value <- .convert_to_data.frame(value, x, "u")

  x <- .update_structure(x, value, "u")

  dupclust <- duplicated(clusters(x))
  dupall <- duplicated(x@structure[x@column_index != "f"])
  if (any(dupclust)) {

    if (sum(dupclust) != sum(dupall)) {
      stop(paste("Fewer new clusters then original, but new collapsed",
                 "clusters would have non-constant treatment and/or",
                 "block structure"))
    }
    warning("Fewer new clusters then original, collapsing")

    x@structure <- x@structure[-dupclust, ]
  }
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

############### Unitid

##' @export
##' @rdname Design_extractreplace
setGeneric("unitids", function(x) standardGeneric("unitids"))

##' @export
##' @rdname Design_extractreplace
setMethod("unitids", "Design", function(x) {
  if (x@unit_of_assignment_type == "cluster") {
    stop("Design specified with `cluster()`, not `unitid()`")
  }
  if (x@unit_of_assignment_type == "unit_of_assignment") {
    stop("Design specified with `unit_of_assignment()`, not `unitid()`")
  }
  return(x@structure[x@column_index == "u"])
})

##' @export
##' @rdname Design_extractreplace
setGeneric("unitids<-", function(x, value) standardGeneric("unitids<-"))

##' @export
##' @rdname Design_extractreplace
setMethod("unitids<-", "Design", function(x, value) {

  value <- .convert_to_data.frame(value, x, "u")

  x <- .update_structure(x, value, "u")

  dupids <- duplicated(unitids(x))
  dupall <- duplicated(x@structure[x@column_index != "f"])
  if (any(dupids)) {

    if (sum(dupids) != sum(dupall)) {
      stop(paste("Fewer new unitids then original, but new collapsed",
                 "units would have non-constant treatment and/or",
                 "block structure"))
    }
    warning("Fewer new unitids then original, collapsing")

    x@structure <- x@structure[-dupids, ]
  }
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

############### Blocks

##' @export
##' @rdname Design_extractreplace
setGeneric("blocks", function(x, newdata = NULL, by = NULL) standardGeneric("blocks"))

##' @export
##' @rdname Design_extractreplace
setMethod("blocks", "Design", function(x, newdata = NULL, by = NULL) {
  .design_accessors_newdata_validate(newdata, by)
  if (!is.null(newdata)) {
    return(.get_col_from_new_data(x, newdata, type = "b", by))
  }
  return(x@structure[x@column_index == "b"])
})

##' @export
##' @rdname Design_extractreplace
setGeneric("blocks<-", function(x, value) standardGeneric("blocks<-"))

##' @export
##' @rdname Design_extractreplace
setMethod("blocks<-", "Design", function(x, value) {
  value <- .convert_to_data.frame(value, x, "b")

  x <- .update_structure(x, value, "b")
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

############### Forcing

##' @export
##' @rdname Design_extractreplace
setGeneric("forcings", function(x, newdata = NULL, by = NULL) standardGeneric("forcings"))

##' @export
##' @rdname Design_extractreplace
setMethod("forcings", "Design", function(x, newdata = NULL, by = NULL) {
  if (x@type != "RD") {
    stop("Forcing variable only used in RD designs")
  }
  .design_accessors_newdata_validate(newdata, by)
  if (!is.null(newdata)) {
    return(.get_col_from_new_data(x, newdata, type = "f", by))
  }
  return(x@structure[x@column_index == "f"])
})

##' @export
##' @rdname Design_extractreplace
setGeneric("forcings<-", function(x, value) standardGeneric("forcings<-"))

##' @export
##' @rdname Design_extractreplace
setMethod("forcings<-", "Design", function(x, value) {
  if (x@type != "RD") {
    stop("Forcing variable only used in RD designs")
  }

  value <- .convert_to_data.frame(value, x, "f")

  x <- .update_structure(x, value, "f")
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

############### dichotomy

##' Extract or replace dichotomy
##' @param x \code{Design} object
##' @param value Replacement \code{dichotomy} formula, or \code{NULL} to remove
##' @return Dichomization formula
##' @export
##' @rdname Design_extract_dichotomy
setGeneric("dichotomy", function(x) standardGeneric("dichotomy"))

##' @export
##' @rdname Design_extract_dichotomy
setMethod("dichotomy", "Design", function(x) {
  return(x@dichotomy)
})

##' @export
##' @rdname Design_extract_dichotomy
setGeneric("dichotomy<-", function(x, value) standardGeneric("dichotomy<-"))

##' @export
##' @rdname Design_extract_dichotomy
setMethod("dichotomy<-", "Design", function(x, value) {
  if (is.null(value)) {
    value <- stats::formula()
  }
  x@dichotomy <- value
  validObject(x)
  return(x)
})


############### Helper Functions

##' Helper function for \code{Design} replacers to ensure replacement is a
##' properly named \code{data.frame}
##'
##' When given a replacement set of values (e.g \code{vector} or \code{matrix}),
##' this ensures that the replacement is a named \code{data.frame}.
##'
##' Input \code{vector}: Since it cannot be named, a vector can only be used to
##' replace an existing component. If the existing component has more than 1
##' column, uses the name of the first column.
##'
##' Input \code{matrix} or \code{data.frame}: If unnamed and replacing existing
##' component, must have no more columns than original component. (If less
##' columns, uses the name of the first few columns.) If named, can replace any
##' number of columns.
##'
##' @title (Internal) Ensures replacement column for \code{Design} is a
##'   \code{data.frame}.
##' @param value A \code{vector} or \code{data.frame} containing a replacement.
##' @param design A \code{Design}
##' @param type One of "t", "f", "u" or "b"
##' @return \code{data.frame} containing named column(s)
##' @keywords internal
.convert_to_data.frame <- function(value, design, type) {
  if (!type %in% c("t", "f", "u", "b")) {
    stop("Invalid type argument")
  }

  # colnames(value) will return NULL if passed a vector, or a matrix/df without
  # names
  if (is.null(colnames(value))) {
    null_name <- TRUE
  } else {
    null_name <- FALSE
  }

  if (!inherits(value, "data.frame")) {
    # Gives a nice error (cannot convert class to data.frame) so no need for
    # tryCatch here
    value <- as.data.frame(value)
  }

  if (nrow(design@structure) != nrow(value)) {
    stop("replacement entries do not have same number of rows as current")
  }

  if (null_name) {
    old_names <- var_names(design, type)
    if (length(old_names) > ncol(value)) {
      # If replacement has less columns, and we need to rename, take the first
      # few names from the original component.
      old_names <- old_names[seq_len(ncol(value))]
    } else if (length(old_names) < ncol(value)) {
      if (length(old_names) == 0) {
        # No original component in design
        stop("Adding a new component requires a named `data.frame`")
      } else {
        # Original component is shorter
        stop(paste("Replacement component has more columns than original,",
                   "replacement must be named"))
      }
    }
    colnames(value) <- old_names
  }

  return(value)
}

##' Assumes \code{.convert_to_data.frame()} has already been called on
##' \code{new}
##' @title (Internal) Replaces \code{type} columns in \code{design} with
##'   \code{new}
##' @param design A \code{Design}
##' @param new A named \code{data.frame} with the replacement, should be the
##'   output of \code{.convert_to_data.frame()}.
##' @param type One of "t", "f", "u" or "b\".
##' @return The updated \code{Design}
##' @keywords internal
.update_structure <- function(design, new, type) {
  design@structure <-
    cbind.data.frame(design@structure[design@column_index != type], new)

  design@column_index <- c(design@column_index[design@column_index != type],
                           rep(type, ncol(new)))
  names(design@column_index) <- colnames(design@structure)
  validObject(design)
  return(design)
  }

##' Helper function to update the \code{call} with the appropriate variable
##' names after they've been modified. Called within \code{Design} replacers.
##'
##' It's return should be stuck into the design via \code{des@call$formula <-
##' .update_call_formula(des)}
##' @title (Internal) Updates `des@call`'s formula with the currently defined
##'   variable names.
##' @param design A \code{Design}
##' @return An updated formula
##' @keywords internal
.update_call_formula <- function(design) {

  .collapse <- function(d, type) {
    paste(var_names(d, type), collapse = ",")
  }

  # Get treatment, ~, and uoa/cluster
  form <- paste0(var_names(design, "t"), "~", design@unit_of_assignment_type,
                 "(", .collapse(design, "u"), ")")

  # Get block if included
  if (length(var_names(design, "b")) > 0) {
    form <- paste0(form, "+ block(", .collapse(design, "b"), ")")
  }

  if (design@type == "RD") {
    form <- paste0(form, "+ forcing(", .collapse(design, "f"), ")")
  }
  return(formula(form))
}

##' @title (Internal) Extract specified \code{type} from new data set
##' @param design A \code{Design}
##' @param newdata A \code{data.frame}, which may or may not be the one which
##'   was used to create \code{design}. It must have the units of assignment
##'   variable(s) (though \code{by=} argument can be used if the name differ),
##'   and will appropriately merge with the \code{design} the blocks, treatment
##'   or forcings.
##' @param type One of "t", "f", or "b".
##' @param by optional; named vector or list connecting names of cluster/unit of
##'   assignment variables in \code{design} to cluster/unit of assignment
##'   variables in \code{data}. Names represent variables in the Design; values
##'   represent variables in the data. Only needed if variable names differ.
##' @return The column(s) belonging to the requested \code{type} in
##' @keywords internal
##' @importFrom stats model.matrix
.get_col_from_new_data <- function(design, newdata, type, by = NULL) {

  if (!is.null(by)) {
    design <- .update_by(design, newdata, by)
  }

  form_for_design <- as.formula(paste("~",
                                       paste(c(var_names(design, "u"),
                                               var_names(design, type)),
                                             collapse = "+"),
                                       " - 1"))
  form_for_newdata <- as.formula(paste("~",
                                       paste(var_names(design, "u"),
                                             collapse = "+"),
                                       " - 1"))

  design_data <- stats::model.matrix(form_for_design, design@structure)
  newdata_data <- stats::model.matrix(form_for_newdata, newdata)

  merged <- merge(design_data, newdata_data, by = var_names(design, "u"))

  return(merged[var_names(design, type)])


}

##' @title (Internal) Checks newdata/by argument for design accessors
##' @param newdata newdata argument from e.g. \code{treatment()},
##'   \code{blocks()}, etc
##' @param byargument from e.g. \code{treatment()},
##'   \code{blocks()}, etc. See \code{.check_by()}
##' @return Invisibly \code{TRUE}. Warns or errors as appropriate.
##' @keywords internal
.design_accessors_newdata_validate <- function(newdata, by) {

  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) {
      warning(paste("`newdata` is not a data.frame, errors or",
                    "unpredictable results may occur"))
    }
  }
  if (!is.null(by)) {
    .check_by(by)
  }

  invisible(TRUE)
}
