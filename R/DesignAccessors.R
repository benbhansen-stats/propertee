
############### Treatment

##' @export
##' @rdname Design_extractreplace
setGeneric("treatment", function(x) standardGeneric("treatment"))

##' @title Extract and Replace elements of Design
##' @param x Design object
##' @return data.frame containing unit-leve information
##' @export
##' @rdname Design_extractreplace
setMethod("treatment", "Design", function(x) {
  x@structure[x@columnIndex == "t"]
})

##' @export
##' @rdname Design_extractreplace
setGeneric("treatment<-", function(x, value) standardGeneric("treatment<-"))

##' @param value Replacement. Either a vector/matrix of appropriate dimension,
##'   or a named data.frame if renaming variable as well.
##' @export
##' @rdname Design_extractreplace
setMethod("treatment<-", "Design", function(x, value) {
  value <- .convert_treatment_to_factor(value)

  value <- .convert_to_data.frame(value, x, "t")

  x@structure[x@columnIndex == "t"] <- value
  names(x@structure)[x@columnIndex == "t"] <- colnames(value)
  validObject(x)
  x
})

############### Units of Assignment

##' @export
##' @rdname Design_extractreplace
setGeneric("unitsOfAssignment", function(x) standardGeneric("unitsOfAssignment"))

##' @export
##' @rdname Design_extractreplace
setMethod("unitsOfAssignment", "Design", function(x) {
  if (x@unitOfAssignmentType == "unitid") {
    stop("Design specified with `unitid()`, not `unitOfAssignment()`")
  }
  if (x@unitOfAssignmentType == "cluster") {
    stop("Design specified with `cluster()`, not `unitOfAssignment()`")
  }
  x@structure[x@columnIndex == "u"]
})

##' @export
##' @rdname Design_extractreplace
setGeneric("unitsOfAssignment<-", function(x, value) standardGeneric("unitsOfAssignment<-"))

##' @export
##' @rdname Design_extractreplace
setMethod("unitsOfAssignment<-", "Design", function(x, value) {

  value <- .convert_to_data.frame(value, x, "u")

  x <- .update_structure(x, value, "u")

  dupclust <- duplicated(unitsOfAssignment(x))
  dupall <- duplicated(x@structure[x@columnIndex != "f"])
  if (any(dupclust)) {

    if (sum(dupclust) != sum(dupall)) {
      stop(paste("Fewer new units of assignment then original, but new collapsed",
                 "units would have non-constant treatment and/or",
                 "block structure"))
    }
    warning("Fewer new units of assignment then original, collapsing")

    x@structure <- x@structure[-dupclust,]
  }
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
  if (x@unitOfAssignmentType == "unitid") {
    stop("Design specified with `unitid()`, not `cluster()`")
  }
  if (x@unitOfAssignmentType == "unitOfAssignment") {
    stop("Design specified with `unitOfAssignment()`, not `cluster()`")
  }
  x@structure[x@columnIndex == "u"]
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
  dupall <- duplicated(x@structure[x@columnIndex != "f"])
  if (any(dupclust)) {

    if (sum(dupclust) != sum(dupall)) {
      stop(paste("Fewer new clusters then original, but new collapsed",
                 "clusters would have non-constant treatment and/or",
                 "block structure"))
    }
    warning("Fewer new clusters then original, collapsing")

    x@structure <- x@structure[-dupclust,]
  }
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
  if (x@unitOfAssignmentType == "cluster") {
    stop("Design specified with `cluster()`, not `unitid()`")
  }
  if (x@unitOfAssignmentType == "unitOfAssignment") {
    stop("Design specified with `unitOfAssignment()`, not `unitid()`")
  }
  x@structure[x@columnIndex == "u"]
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
  dupall <- duplicated(x@structure[x@columnIndex != "f"])
  if (any(dupids)) {

    if (sum(dupids) != sum(dupall)) {
      stop(paste("Fewer new unitids then original, but new collapsed",
                 "units would have non-constant treatment and/or",
                 "block structure"))
    }
    warning("Fewer new unitids then original, collapsing")

    x@structure <- x@structure[-dupids,]
  }
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
  x@structure[x@columnIndex == "b"]
})

##' @export
##' @rdname Design_extractreplace
setGeneric("blocks<-", function(x, value) standardGeneric("blocks<-"))

##' @export
##' @rdname Design_extractreplace
setMethod("blocks<-", "Design", function(x, value) {
  value <- .convert_to_data.frame(value, x, "b")

  .update_structure(x, value, "b")
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
  x@structure[x@columnIndex == "f"]
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

  .update_structure(x, value, "f")
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
      nullName <- TRUE
    } else {
      nullName <- FALSE
    }
    value <- as.data.frame(value)
    if (nrow(design@structure) != nrow(value)) {
      stop("replacement entries do not have same number of rows as current")
    }
    if (nullName) {
      oldNames <- varNames(design, type)
      if (length(oldNames) > ncol(value)) {
        oldNames <- oldNames[seq_len(ncol(value))]
      } else if (length(oldNames) < ncol(value)) {
        stop("additional variables must be named")
      }
      colnames(value) <- oldNames
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
    cbind.data.frame(design@structure[design@columnIndex != type], new)

  design@columnIndex<- c(design@columnIndex[design@columnIndex != type],
                         rep(type, ncol(new)))
  names(design@columnIndex) <- colnames(design@structure)
  validObject(design)
  return(design)
  }

# Internal helper function
# Converts treatment to factor
.convert_treatment_to_factor <- function(treatment) {
  if (!(is.data.frame(treatment) |
          (is.null(dim(treatment)) &
             (is.numeric(treatment) |
                is.factor(treatment) |
                is.logical(treatment))))) {
    stop("Treatment must be numeric/factor/logical vector or data.frame")
  }

  return_data_frame <- FALSE
  if (!is.null(dim(treatment))) {
    if (ncol(treatment) != 1) {
      stop("Only one treatment variable allowed")
    }
    # If the treatment is named (e.g. column in DF), we'll keep it later...
    return_data_frame <- TRUE
  }
  else {
    # ... if its not named, convert to DF for processing first
    treatment <- data.frame(treatment)
  }

  if (is.numeric(treatment[,1])) {
    if (any(!treatment[,1] %in% 0:1)) {
      stop("Numerical treatments must only contain values 0 and 1.")
    }
    treatment[,1] <- as.factor(treatment[,1])
  } else if (is.logical(treatment[,1])) {
    treatment[,1] <-  as.factor(treatment[,1])
  } else if (!is.factor(treatment[,1])) {
    stop("Treatment must be binary (0/1), logical, factor or ordered")
  }
  treatment[,1] <- droplevels(treatment[,1])
  if (!return_data_frame) {
    # if treatment wasn't passed as a DF, drop dimension
    treatment <- treatment[,1]
  }
  return(treatment)
}
