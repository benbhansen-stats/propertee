#' @include StudySpecification.R
NULL

############### Treatment

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("treatment",
           function(x, newdata = NULL, dichotomy = NULL, by = NULL, ...) {
  standardGeneric("treatment")
})

##' @title Accessors and Replacers for \code{StudySpecification} objects
##'
##' @description Allows access to the elements which define a
##'   \code{StudySpecification}, enabling their extraction or replacement.
##'
##' @details For [treatment()], when argument \code{binary} is \code{FALSE}, the
##'   treatment variable passed into the \code{StudySpecification} is returned
##'   as a one-column \code{data.frame} regardless of whether it is binary or
##'   \code{x} has a \code{dichotomy}
##'
##'   If a \code{dichotomy} is passed, a binary one-column \code{data.frame}
##'   will be returned. If not and \code{binary} is \code{TRUE}, unless the
##'   \code{StudySpecification} has a binary treatment, [treatment()] will
##'   error. If \code{binary} is \code{"ifany"}, it will return the original
##'   treatment in this case.
##'
##'   The one-column \code{data.frame} returned by [treatment()] is named as
##'   entered in the \code{StudySpecification} creation, but if a
##'   \code{dichotomy} is passed, the column name is \code{"__z"} to try and
##'   avoid any name conflicts.
##'
##'   For the \code{value} when using replacers, the replacement must have the
##'   same number of rows as the \code{StudySpecification} (the same number of
##'   units of assignment). The number of columns can differ (e.g. if the
##'   \code{StudySpecification} were defined with two variable uniquely
##'   identifying blocks, you can replace that with a single variable uniquely
##'   identifying blocks, as long as it respects other restrictions.)
##'
##'   If the replacement value is a \code{data.frame}, the name of the columns
##'   is used as the new variable names. If the replacement is a \code{matrix}
##'   or \code{vector}, the original names are retained. If reducing the number
##'   of variables (e.g., moving from two variables uniquely identifying to a
##'   single variable), the appropriate number of variable names are retained.
##'   If increasing the number of variables, a \code{data.frame} with names must
##'   be provided.
##'
##' @param x a \code{StudySpecification} object
##' @param newdata optional; an additional \code{data.frame}. If passed, and the
##'   unit of assignment variable is found in \code{newdata}, then the requested
##'   variable type for each unit of \code{newdata} is returned. See \code{by}
##'   argument if the name of the unit of assignment differs.
##' @param dichotomy optional; a formula specifying how to dichotomize a
##'   non-binary treatment variable. See the Details section of the \code{ett()}
##'   or \code{att()} help pages for information on specifying this formula
##' @param by optional; named vector or list connecting names of unit of
##'   assignment/unitid/cluster variables in \code{x} to unit of
##'   assignment/unitid/cluster variables in \code{data}. Names represent
##'   variables in \code{x}; values represent variables in \code{newdata}. Only
##'   needed if variable names differ.
##' @param implicit Should a block-less \code{StudySpecification} return a
##'   constant 1 when extracting \code{blocks}?
##' @param ... ignored.
##' @return \code{data.frame} containing requested variable, or an updated
##'   \code{StudySpecification}. [treatment()] works slightly differently, see
##'   \code{Details}.
##' @export
##' @rdname StudySpecification_extractreplace
##' @examples
##' data(simdata)
##' spec <- obs_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
##' blocks(spec) # empty
##' blocks(spec) <- data.frame(blks = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5))
##' blocks(spec)
##' blocks(spec) <- c(5, 5, 4, 4, 3, 3, 2, 2, 1, 1)
##' blocks(spec) # notice that variable is not renamed
setMethod("treatment", "StudySpecification",
          function(x, newdata = NULL, dichotomy = NULL, by = NULL, ...) {
  .specification_accessors_newdata_validate(newdata, by)

  if (is.null(dichotomy)) {
    if (!is.null(newdata)) {
      return(.get_col_from_new_data(x, newdata, type = "t", by))
    } else {
      return(x@structure[x@column_index == "t"])
    }
  }

  if (!is.null(dichotomy)) {
    return(data.frame(z__ = .bin_txt(x, newdata, dichotomy)))
  } else if (!has_binary_treatment(x)) {
    stop(paste("No binary treatment can be produced. Treatment is",
               "non-binary and no `dichotomy` provided."))
  }
})

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("treatment<-", function(x, value) standardGeneric("treatment<-"))

##' @param value replacement. Either a \code{vector}/\code{matrix} of
##'   appropriate dimension, or a named \code{data.frame} if renaming variable
##'   as well. See \code{Details}.
##' @export
##' @rdname StudySpecification_extractreplace
setMethod("treatment<-", "StudySpecification", function(x, value) {
  value <- .convert_to_data.frame(value, x, "t")

  x@structure[x@column_index == "t"] <- value
  names(x@structure)[x@column_index == "t"] <- colnames(value)
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

##' @title (Internal) Extracts treatment as binary \code{vector}
##' @details If a \code{dichotomy} is specified or the \code{StudySpecification}
##'   has a treatment variable consisting only of 0/1 or \code{NA}, then returns
##'   the binary treatment. Otherwise (it has a non-binary treatment and lacks a
##'   dichotomy) it errors.
##' @param spec A \code{StudySpecification}, used to get treatment assignment
##'   information
##' @param data A dataframe with unit of assignment information and, if a
##'   dichotomy is provided, columns specified therein
##' @param dichotomy Optional, a formula. See the Details section of the
##'   \code{ett()} or \code{att()} help pages for information on specifying the
##'   formula
##' @return A \code{vector} of binary treatments
##' @keywords internal
.bin_txt <- function(spec, data = NULL, dichotomy = NULL) {
  # get treatment from the specification
  tt <- treatment(spec)

  if (!is.null(data)) tt <- .expand_txt(tt, data, spec)

  if (!is.null(dichotomy)) {
    treatment <- .apply_dichotomy(tt, dichotomy)
  } else {
    treatment <- tt[,1]
  }

  if (!all(treatment %in% c(0:1, NA))) {
    stop("Must provide a dichotomy if the `StudySpecification` has a non-binary treatment")
  }

  return(treatment)
}

##' @title (Internal) Expand treatment variable from a \code{StudySpecification}
##'   to a dataframe with unit of assignment information
##' @param txt A dataframe with one column corresponding to the treatment. Can
##'   be dichotomized or as it's stored in \code{spec}
##' @param data A dataframe with unit of assignment information
##' @param spec A \code{StudySpecification}, used to align unit of assignment
##'   information with \code{txt}
##' @keywords internal
.expand_txt <- function(txt, data, spec) {
  # if there's a `data` argument, merge it to the treatment info using the
  # units of assignment
  if (!all(var_names(spec, "u") %in% colnames(data))) {
    stop("Not all unit of assignment variables can be found in `data`")
  }
  treatment_uoa <- cbind(txt, spec@structure[, var_names(spec, "u"), drop = FALSE])
  txt <- .merge_preserve_order(data, treatment_uoa, by = var_names(spec, "u"), all.x = TRUE)

  txtname <- var_names(spec, "t")
  txt <- tryCatch(txt[, txtname, drop = FALSE],
                  error = function(e) {
                    # if treatment variable already exists in data, there
                    # will be a .x and .y version; e.g. z.x and z.y, so
                    # we'll extract the ".y" version (the second one)
                    # since the merge above has the treatment from the
                    # StudySpecification second.
                    txt[, paste0(txtname, ".y"), drop = FALSE]
                  })
  colnames(txt) <- txtname

  return(txt)
}

##' @title (Internal) Applies dichotomy to treatment
##' @description Given a dichotomy formula and a \code{data.frame} with a
##'   treatment variable and any variables in the formula, returns a
##'   \code{vector} containing only \code{0}, \code{1}, or \code{NA}.
##' @param txt A named \code{data.frame} containing a column of the treatment,
##'   such as that produed by `treatment(myspecification)`, and any variables
##'   specified in \code{dichotomy}.
##' @param dichotomy A formula specifying how to dichotomize the non-binary
##'   treatment column in \code{txt} (or a call that evaluates to a formula).
##'   See the Details section of the \code{ett()} or \code{att()} help pages for
##'   information on specifying this formula
##' @return A \code{vector} of binary treatments
##' @keywords internal
.apply_dichotomy <- function(txt, dichotomy) {
  if (!is.data.frame(txt)) stop("`txt` is expected to be a named `data.frame`")
  if (!all(setdiff(all.vars(dichotomy), ".") %in% colnames(txt))) {
    stop(paste("Could not find variables specified in dichotomy. Provide a",
               "data argument with these columns, or ensure provided data argument",
               "has them."))
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
##' @rdname StudySpecification_extractreplace
setGeneric("units_of_assignment", function(x, newdata = NULL, by = NULL) {
  standardGeneric("units_of_assignment")
})

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("units_of_assignment", "StudySpecification", function(x, newdata = NULL, by = NULL) {
  if (x@unit_of_assignment_type == "none") {
    stop("StudySpecification specified without unit of assignment")
  }
  if (x@unit_of_assignment_type == "unitid") {
    stop("StudySpecification specified with `unitid()`, not `unit_of_assignment()`")
  }
  if (x@unit_of_assignment_type == "cluster") {
    stop("StudySpecification specified with `cluster()`, not `unit_of_assignment()`")
  }
  .specification_accessors_newdata_validate(newdata, by)
  if (!is.null(newdata)) {
    return(.get_col_from_new_data(x, newdata, type = "u", by))
  }
  return(x@structure[x@column_index == "u"])
})

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("units_of_assignment<-", function(x, value) {
  standardGeneric("units_of_assignment<-")
})

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("units_of_assignment<-", "StudySpecification", function(x, value) {

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
##' @rdname StudySpecification_extractreplace
setGeneric("clusters", function(x, newdata = NULL, by = NULL) standardGeneric("clusters"))

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("clusters", "StudySpecification", function(x, newdata = NULL, by = NULL) {
  if (x@unit_of_assignment_type == "none") {
    stop("StudySpecification specified without clusters")
  }
  if (x@unit_of_assignment_type == "unitid") {
    stop("StudySpecification specified with `unitid()`, not `cluster()`")
  }
  if (x@unit_of_assignment_type == "unit_of_assignment") {
    stop("StudySpecification specified with `unit_of_assignment()`, not `cluster()`")
  }
  .specification_accessors_newdata_validate(newdata, by)
  if (!is.null(newdata)) {
    return(.get_col_from_new_data(x, newdata, type = "u", by))
  }
  return(x@structure[x@column_index == "u"])
})

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("clusters<-", function(x, value) standardGeneric("clusters<-"))

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("clusters<-", "StudySpecification", function(x, value) {

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
##' @rdname StudySpecification_extractreplace
setGeneric("unitids", function(x) standardGeneric("unitids"))

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("unitids", "StudySpecification", function(x) {
  if (x@unit_of_assignment_type == "none") {
    stop("StudySpecification specified without unit IDs")
  }
  if (x@unit_of_assignment_type == "cluster") {
    stop("StudySpecification specified with `cluster()`, not `unitid()`")
  }
  if (x@unit_of_assignment_type == "unit_of_assignment") {
    stop("StudySpecification specified with `unit_of_assignment()`, not `unitid()`")
  }
  return(x@structure[x@column_index == "u"])
})

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("unitids<-", function(x, value) standardGeneric("unitids<-"))

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("unitids<-", "StudySpecification", function(x, value) {

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

############### ..uoa..
# This is when a spec is specified without uoa/cluster/unitid, e.g. z ~ 1. See
# #193 for details

##' @title Return ..uoa.. column
##' @param spec A studyspecification
##' @return The ..uoa.. column
##' @keywords internal
..uoa.. <- function(spec) {
  if (spec@unit_of_assignment_type != "none") {
    stop("Internal error: ..uoa..() only valid for `unit_of_assignment_type` == 'none'")
  }
  return(spec@structure["..uoa.."])
}

############### Blocks

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("blocks", function(x,
                              newdata = NULL,
                              by = NULL,
                              ...) standardGeneric("blocks"))

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("blocks", "StudySpecification", function(x,
                                       newdata = NULL,
                                       by = NULL,
                                       ...,
                                       implicit = FALSE) {
  .specification_accessors_newdata_validate(newdata, by)
  if (!is.null(newdata)) {
    return(.get_col_from_new_data(x,
                                  newdata,
                                  type = "b",
                                  by,
                                  implicitBlock = implicit,
                                  ...))
  }
  if (!any(x@column_index == "b") & implicit) {
    return(data.frame(".blocks_internal" = rep(1, nrow(x@structure))))
  }
  return(x@structure[x@column_index == "b"])
})

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("blocks<-", function(x, value) standardGeneric("blocks<-"))

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("blocks<-", "StudySpecification", function(x, value) {
  value <- .convert_to_data.frame(value, x, "b")

  x <- .update_structure(x, value, "b")
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

##' @export
##' @rdname StudySpecification_extractreplace
has_blocks <- function(x) {
  return("b" %in% x@column_index)
}

############### Forcing

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("forcings", function(x, newdata = NULL, by = NULL) standardGeneric("forcings"))

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("forcings", "StudySpecification", function(x, newdata = NULL, by = NULL) {
  if (x@type != "RD") {
    stop("Forcing variable only used in RD specifications")
  }
  .specification_accessors_newdata_validate(newdata, by)
  if (!is.null(newdata)) {
    return(.get_col_from_new_data(x, newdata, type = "f", by))
  }
  return(x@structure[x@column_index == "f"])
})

##' @export
##' @rdname StudySpecification_extractreplace
setGeneric("forcings<-", function(x, value) standardGeneric("forcings<-"))

##' @export
##' @rdname StudySpecification_extractreplace
setMethod("forcings<-", "StudySpecification", function(x, value) {
  if (x@type != "RD") {
    stop("Forcing variable only used in RD specifications")
  }

  value <- .convert_to_data.frame(value, x, "f")

  x <- .update_structure(x, value, "f")
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  return(x)
})

############### Helper Functions

##' Helper function for \code{StudySpecification} replacers to ensure
##' replacement is a properly named \code{data.frame}
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
##' @title (Internal) Ensures replacement column for \code{StudySpecification}
##'   is a \code{data.frame}.
##' @param value A \code{vector} or \code{data.frame} containing a replacement.
##' @param specification A \code{StudySpecification}
##' @param type One of "t", "f", "u" or "b"
##' @return \code{data.frame} containing named column(s)
##' @keywords internal
.convert_to_data.frame <- function(value, specification, type) {
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

  if (nrow(specification@structure) != nrow(value)) {
    stop("replacement entries do not have same number of rows as current")
  }

  if (null_name) {
    old_names <- var_names(specification, type)
    if (length(old_names) > ncol(value)) {
      # If replacement has less columns, and we need to rename, take the first
      # few names from the original component.
      old_names <- old_names[seq_len(ncol(value))]
    } else if (length(old_names) < ncol(value)) {
      if (length(old_names) == 0) {
        # No original component in specification
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
##' @title (Internal) Replaces \code{type} columns in \code{specification} with
##'   \code{new}
##' @param specification A \code{StudySpecification}
##' @param new A named \code{data.frame} with the replacement, should be the
##'   output of \code{.convert_to_data.frame()}.
##' @param type One of "t", "f", "u" or "b\".
##' @return The updated \code{StudySpecification}
##' @keywords internal
.update_structure <- function(specification, new, type) {
  specification@structure <-
    cbind.data.frame(specification@structure[specification@column_index != type], new)

  specification@column_index <- c(specification@column_index[specification@column_index != type],
                           rep(type, ncol(new)))
  names(specification@column_index) <- colnames(specification@structure)
  validObject(specification)
  return(specification)
  }

##' Helper function to update the \code{call} with the appropriate variable
##' names after they've been modified. Called within \code{StudySpecification}
##' replacers.
##'
##' It's return should be stuck into the specification via
##' \code{spec@call$formula <- .update_call_formula(spec)}
##' @title (Internal) Updates `spec@call`'s formula with the currently defined
##'   variable names.
##' @param specification A \code{StudySpecification}
##' @return An updated formula
##' @keywords internal
.update_call_formula <- function(specification) {

  .collapse <- function(d, type) {
    paste(var_names(d, type), collapse = ",")
  }

  # Get treatment, ~, and uoa/cluster
  form <- paste0(var_names(specification, "t"), "~", specification@unit_of_assignment_type,
                 "(", .collapse(specification, "u"), ")")

  # Get block if included
  if (length(var_names(specification, "b")) > 0) {
    form <- paste0(form, "+ block(", .collapse(specification, "b"), ")")
  }

  if (specification@type == "RD") {
    form <- paste0(form, "+ forcing(", .collapse(specification, "f"), ")")
  }
  return(formula(form))
}

##' @title (Internal) Extract specified \code{type} from new data set
##' @param specification A \code{StudySpecification}
##' @param newdata A \code{data.frame}, which may or may not be the one which
##'   was used to create \code{specification}. It must have the units of
##'   assignment variable(s) (though \code{by=} argument can be used if the name
##'   differ), and will appropriately merge with the \code{specification} the
##'   blocks, treatment or forcings.
##' @param type One of "t", "f", or "b".
##' @param by optional; named vector or list connecting names of unit of
##'   assignment/unitid/cluster variables in \code{specification} to unit of
##'   assignment/unitid/cluster variables in \code{data}. Names represent
##'   variables in the StudySpecification; values represent variables in the
##'   data. Only needed if variable names differ.
##' @param implicitBlock If the \code{StudySpecification} does not include a
##'   block, \code{TRUE} will return a constant 1 for the blocks if \code{type}
##'   requests it.
##' @param ... Additional arguments to \code{merge()}.
##' @return The column(s) belonging to the requested \code{type} in
##' @keywords internal
##' @importFrom stats model.matrix
.get_col_from_new_data <- function(specification,
                                   newdata,
                                   type,
                                   by = NULL,
                                   implicitBlock = FALSE,
                                   ...) {

  if (!is.null(by)) {
    specification <- .update_by(specification, newdata, by)
  }

  form_for_spec <- as.formula(paste("~",
                                      paste(c(var_names(specification, "u"),
                                              var_names(specification, type)),
                                            collapse = "+"),
                                      " - 1"))
  form_for_newdata <- as.formula(paste("~",
                                       paste(var_names(specification, "u"),
                                             collapse = "+"),
                                       " - 1"))

  specification_data <- stats::model.frame(form_for_spec, specification@structure)
  if (type == "b" &&
        var_names(specification, "b", implicitBlocks = TRUE)[1] == ".blocks_internal") {
    specification_data$.blocks_internal <- rep(1, nrow(specification_data))
  }

  if (specification@unit_of_assignment_type == "none") {
    newdata[["..uoa.."]] <- rownames(newdata)
  }



  newdata_data <- stats::model.frame(form_for_newdata, newdata, na.action = na.pass)

  merged <- .merge_preserve_order(newdata_data,
                                  specification_data,
                                  by = var_names(specification, "u"),
                                  sort = FALSE,
                                  ...)

  return(merged[var_names(specification, type, implicitBlocks = TRUE)])


}

##' @title (Internal) Checks newdata/by argument for specification accessors
##' @param newdata newdata argument from e.g. \code{treatment()},
##'   \code{blocks()}, etc
##' @param by from e.g. \code{treatment()}, \code{blocks()}, etc. See
##'   \code{.check_by()}
##' @return Invisibly \code{TRUE}. Warns or errors as appropriate.
##' @keywords internal
.specification_accessors_newdata_validate <- function(newdata, by) {

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
