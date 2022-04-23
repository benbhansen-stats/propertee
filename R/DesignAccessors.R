############### Treatment

##' @export
##' @rdname Design_extractreplace
setGeneric("treatment", function(x, binary = FALSE) standardGeneric("treatment"))

##' Extract and replace elements of the Design
##'
##' For \code{treatment()}, the value of \code{binary} also affects the
##' dimensions of what is return. \code{treatment(..., binary = FALSE)} (the
##' default) will return a \code{data.frame} with the treatment as a column, for
##' consistency with the other accessors. However, for \code{treatment(...,
##' binary = TRUE)}, the output is a vector.
##' @title Accessors for Design
##' @param x Design object
##' @param binary Logical; if FALSE (default), return a `data.frame` containing
##'   the named treatment variable. If TRUE and \code{x} has a formula in
##'   \code{@dichotomy}, return a `data.frame` containing a binary
##'   treatment variable with the name "z__". Has no effect if
##'   \code{@dichotomy} is \code{NULL}.
##' @return If \code{binary = FALSE} data.frame containing treatment variable.
##'   If \code{binary = TRUE}, a vector of 0/1 indicating treatment group
##'   membership.
##' @export
##' @rdname Design_extractreplace
setMethod("treatment", "Design", function(x, binary = FALSE) {
  if (binary & is_dichotomized(x)) {
    return(data.frame(z__ = .bin_txt(x)))
  }
  x@structure[x@column_index == "t"]
})

##' @export
##' @rdname Design_extractreplace
setGeneric("treatment<-", function(x, value) standardGeneric("treatment<-"))

##' @param value Replacement. Either a vector/matrix of appropriate dimension,
##'   or a named data.frame if renaming variable as well.
##' @export
##' @rdname Design_extractreplace
setMethod("treatment<-", "Design", function(x, value) {
  value <- .convert_to_data.frame(value, x, "t")

  x@structure[x@column_index == "t"] <- value
  names(x@structure)[x@column_index == "t"] <- colnames(value)
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  x
})

# (Internal) Extracts treatment as binary vector if possible (if design has
# dichotomy, or if treatment is 0/1), or else errors.
.bin_txt <- function(des) {
  if (!is_dichotomized(des)) {
    tt <- treatment(des, binary = FALSE)[, 1]
    if (!all(tt %in% 0:1)) {
      stop("binary treatment cannot be obtained")
    }
    return(tt)
  }
  .binarize_treatment(treatment(des, binary = FALSE),
                      des@dichotomy)
}

############### Units of Assignment

##' @export
##' @rdname Design_extractreplace
setGeneric("units_of_assignment", function(x) standardGeneric("units_of_assignment"))

##' @export
##' @rdname Design_extractreplace
setMethod("units_of_assignment", "Design", function(x) {
  if (x@unit_of_assignment_type == "unitid") {
    stop("Design specified with `unitid()`, not `unit_of_assignment()`")
  }
  if (x@unit_of_assignment_type == "cluster") {
    stop("Design specified with `cluster()`, not `unit_of_assignment()`")
  }
  x@structure[x@column_index == "u"]
})

##' @export
##' @rdname Design_extractreplace
setGeneric("units_of_assignment<-", function(x, value) standardGeneric("units_of_assignment<-"))

##' @export
##' @rdname Design_extractreplace
setMethod("units_of_assignment<-", "Design", function(x, value) {

  value <- .convert_to_data.frame(value, x, "u")

  x <- .update_structure(x, value, "u")

  dupuoa <- duplicated(units_of_assignment(x))
  dupall <- duplicated(x@structure[x@column_index != "f"])
  if (any(dupuoa)) {

    if (sum(dupuoa) != sum(dupall)) {
      stop(paste("Fewer new units of assignment then original, but new collapsed",
                 "units would have non-constant treatment and/or",
                 "block structure"))
    }
    warning("Fewer new units of assignment then original, collapsing")

    x@structure <- x@structure[-dupuoa, ]
  }
  x@call$formula <- .update_call_formula(x)
  validObject(x)
  x
})

############### Cluster

##' @export
##' @rdname Design_extractreplace
setGeneric("clusters", function(x) standardGeneric("clusters"))

##' @export
##' @rdname Design_extractreplace
setMethod("clusters", "Design", function(x) {
  if (x@unit_of_assignment_type == "unitid") {
    stop("Design specified with `unitid()`, not `cluster()`")
  }
  if (x@unit_of_assignment_type == "unit_of_assignment") {
    stop("Design specified with `unit_of_assignment()`, not `cluster()`")
  }
  x@structure[x@column_index == "u"]
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
  x
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
  x@structure[x@column_index == "u"]
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
  x
})

############### Blocks

##' @export
##' @rdname Design_extractreplace
setGeneric("blocks", function(x) standardGeneric("blocks"))

##' @export
##' @rdname Design_extractreplace
setMethod("blocks", "Design", function(x) {
  x@structure[x@column_index == "b"]
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
  x
})

############### Forcing

##' @export
##' @rdname Design_extractreplace
setGeneric("forcings", function(x) standardGeneric("forcings"))

##' @export
##' @rdname Design_extractreplace
setMethod("forcings", "Design", function(x) {
  if (x@type != "RD") {
    stop("Forcing variable only used in RD designs")
  }
  x@structure[x@column_index == "f"]
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
  x
})

############### dichotomy

##' Extract or replace dichotomy
##' @param x Design object
##' @param value Replacement dichotomy formula, or \code{NULL} to remove
##' @return Dichomization formula
##' @export
##' @rdname Design_extract_dichotomy
setGeneric("dichotomy", function(x) standardGeneric("dichotomy"))

##' @export
##' @rdname Design_extract_dichotomy
setMethod("dichotomy", "Design", function(x) {
  x@dichotomy
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
  x
})


############### Helper Functions

# Internal helper function
# Takes in a replacement item and a design,
# and if replacement isn't a data.frame already,
# converts it to a data.frame, extracting names from
# the design if original value isn't already named
.convert_to_data.frame <- function(value, design, type) {
  if (!is(value, "data.frame")) {
    if (is.null(colnames(value))) {
      null_name <- TRUE
    } else {
      null_name <- FALSE
    }
    value <- as.data.frame(value)
    if (nrow(design@structure) != nrow(value)) {
      stop("replacement entries do not have same number of rows as current")
    }
    if (null_name) {
      old_names <- var_names(design, type)
      if (length(old_names) > ncol(value)) {
        old_names <- old_names[seq_len(ncol(value))]
      } else if (length(old_names) < ncol(value)) {
        stop("additional variables must be named")
      }
      colnames(value) <- old_names
    }
  }
  if (nrow(design@structure) != nrow(value)) {
    stop("replacement entries do not have same number of rows as current")
  }
  value
}

# Internal helper function
# Replaces `type` columns in `design` with `new`. Assumes
# `.convert_to_data.frame` has already been called on `new`
.update_structure <- function(design, new, type) {
  design@structure <-
    cbind.data.frame(design@structure[design@column_index != type], new)

  design@column_index <- c(design@column_index[design@column_index != type],
                           rep(type, ncol(new)))
  names(design@column_index) <- colnames(design@structure)
  validObject(design)
  return(design)
  }

# (Internal) Updates `des@call`'s formula with the currently defined variable
# names. Returns the updated formula, should be stuck into the design via
# `des@call$formula <- .update_call_formula(des)`
.update_call_formula <- function(design) {

  .collapse <- function(d, type) {
    paste(var_names(d, type), collapse = ",")
  }

  # Get treatment, ~, and uoa/cluster
  form <- paste0(var_names(design, "t"), "~", design@unit_of_assignment_type, "(",
                 .collapse(design, "u"), ")")

  # Get block if included
  if (length(var_names(design, "b")) > 0) {
    form <- paste0(form, "+ block(", .collapse(design, "b"), ")")
  }

  if (design@type == "RD") {
    form <- paste0(form, "+ forcing(", .collapse(design, "f"), ")")
  }
  form <- formula(form)
}
# (Internal) Updates `des@call`'s formula with the currently defined variable
# names. Returns the updated formula, should be stuck into the design via
# `des@call$formula <- .update_call_formula(des)`
.update_call_formula <- function(design) {

  .collapse <- function(d, type) {
    paste(var_names(d, type), collapse = ",")
  }

  # Get treatment, ~, and uoa/cluster
  form <- paste0(var_names(design, "t"), "~", design@unit_of_assignment_type, "(",
                 .collapse(design, "u"), ")")

  # Get block if included
  if (length(var_names(design, "b")) > 0) {
    form <- paste0(form, "+ block(", .collapse(design, "b"), ")")
  }

  if (design@type == "RD") {
    form <- paste0(form, "+ forcing(", .collapse(design, "f"), ")")
  }
  form <- formula(form)
}
