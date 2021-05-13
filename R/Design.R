setClass("Design",
         slots = c(structure = "data.frame",
                   columnIndex = "character",
                   type = "character"))

setValidity("Design", function(object) {
  if (any(dim(object@structure) == 0)) {
    return("@structure must have positive dimensions")
  }
  tr <- object@structure[object@columnIndex == "t"]
  if (ncol(tr) == 0) {
    return("Missing treatment index")
  }
  tr <- tr[,1]
  if (is.null(tr) || is.na(sd(tr)) || sd(tr) == 0 || any(!(tr %in% 0:1))) {
    return("Invalid treatment; must be binary and non-constant")
  }
  if (ncol(object@structure) != length(object@columnIndex)) {
    return("@columnIndex does not agree with number of columns in @structure")
  }
  if (!all(object@columnIndex %in% c("t", "c", "b", "f"))) {
    wrong <- object@columnIndex[!object@columnIndex %in% c("t", "c", "b", "f")]
    return(paste("@columnIndex design elements must be [t,c,b,f]. unknown elements:",
                 paste(wrong, collapse = ", ")))
  }
  if (!object@type %in% c("RCT", "RD", "Obs")) {
    return(paste("@type must be one of [RCT,RD,Obs]. unknown @type:", object@type))
  }
  TRUE
})


New_Design <- function(form, data, type, subset = NULL) {
  if (!is.null(subset)) {
    data <- subset(data, subset = subset)
  }

  m <- as.data.frame(as.matrix(model.frame(form, data)))

  index <- rep("t", ncol(m))

  # Handle clusters
  clusters <- grepl("^cluster", colnames(m))
  if (any(clusters)) {
    index[which(clusters)] <- "c"
    cvars <- colnames(m)[clusters][1]
    cvars <- sub("^cluster\\(", "", cvars)
    cvars <- sub("\\)[\\.0-9]*$", "", cvars)
    cvars <- gsub(" ", "", cvars)
    colnames(m)[clusters] <- strsplit(cvars, ",")[[1]]
  }

  # Handle blocks
  blocks <- grepl("^block", colnames(m))
  if (any(blocks)) {
    index[which(blocks)] <- "b"
    bvars <- colnames(m)[blocks][1]
    bvars <- sub("^block\\(", "", bvars)
    bvars <- sub("\\)[\\.0-9]*$", "", bvars)
    bvars <- gsub(" ", "", bvars)
    colnames(m)[blocks] <- strsplit(bvars, ",")[[1]]
  }

  # Handle forcing
  forcings <- grepl("^forcing", colnames(m))
  if (any(forcings)) {
    index[which(forcings)] <- "f"
    fvars <- colnames(m)[forcings][1]
    fvars <- sub("^forcing\\(", "", fvars)
    fvars <- sub("\\)[\\.0-9]*$", "", fvars)
    fvars <- gsub(" ", "", fvars)
    colnames(m)[forcings] <- strsplit(fvars, ",")[[1]]
  }

  m_collapse <- unique(m)

  differing <- duplicated(m_collapse[, index == "c"])
  if (any(differing)) {
    stop("Each of treatment assignment, block and forcing must be constant within cluster.")
    # TODO: Can we make this error more informative since we're doing the work?
    # At a minimum identify clusters where the issue arises; or even on what particular
    # variables?
  }

  rownames(m_collapse) <- NULL

  new("Design",
      structure = m_collapse,
      columnIndex = index,
      type = type)
}

##' Generates an RCT Design object with the given specifications.
##'
##' The `formula` must include `cluster()` to identify the units of assignment
##' (one or more variables), it may optionally contain `strata()` as well.
##' @title Specify RCT Design
##' @param formula defines the design components
##' @param data the data set.
##' @param subset optionally subset the data before creating the design object
##' @return a Design object of type "RCT" for use in further analysis
##' @export
RCT_Design <- function(formula, data, subset = NULL) {
  checkDesignFormula(formula)

  design <- New_Design(formula, data, type = "RCT", subset = subset)
  return(design)
}

##' Generates an RD Design object with the given specifications.
##'
##' The `formula` must include `cluster()` to identify the units of assignment
##' (one or more variables), it may optionally contain `strata()` and/or
##' `forcing()` as well.
##' @title Specify Regression Discontinuity Design
##' @param formula defines the design components
##' @param data the data set.
##' @param subset optionally subset the data before creating the design object
##' @return a Design object of type "RD" for use in further analysis
##' @export
RD_Design <- function(formula, data, subset = NULL) {
  checkDesignFormula(formula, allowForcing = TRUE)

  design <- New_Design(formula, data, type = "RD", subset = subset)
  return(design)
}

##' Generates an Observationl Data Design object with the given specifications.
##'
##' The `formula` must include `cluster()` to identify the units of assignment
##' (one or more variables), it may optionally contain `strata()` as well.
##' @title Specify Observational Data Design
##' @param formula defines the design components
##' @param data the data set.
##' @param subset optionally subset the data before creating the design object
##' @return a Design object of type "Obs" for use in further analysis
##' @export
Obs_Design <- function(formula, data, subset = NULL) {
  checkDesignFormula(formula)

  design <- New_Design(formula, data, type = "Obs", subset = subset)
  return(design)
}

##' @title Show a Design
##' @param object Design object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "Design", function(object) {
  destype <- switch(object@type,
                    "RCT" = "Randomized Control Trial",
                    "RD" = "Regression Discontinuity Design",
                    "Obs" = "Observational Study")
  cat(destype)
  cat("\n\n")
  cat(paste("Treatment:", varNames(object, "t")))
  cat("\n")
  cat(paste("Cluster  :", paste(varNames(object, "c"), collapse = ", ")))
  cat("\n")
  if (length(varNames(object, "b")) > 0) {
    cat(paste("Block    :", paste(varNames(object, "b"), collapse = ", ")))
    cat("\n")
  }
  if (length(varNames(object, "f")) > 0) {
    cat(paste("Forcing  :", paste(varNames(object, "f"), collapse = ", ")))
    cat("\n")
  }
  invisible(object)
})


##' @title Extract names of Design variables
##' @param x Design x
##' @param type one of "t", "c", "b", "f"; for "treatment", "cluster", "block",
##'   and "forcing"
##' @return character vector of variable names of the given type
##' @export
varNames <- function(x, type) {
  stopifnot(class(x) == "Design")
  stopifnot(length(type) == 1)
  stopifnot(type %in% c("t", "c", "b", "f"))
  names(x@structure)[x@columnIndex == type]
}


###########################################
##### Accessor and modifier functions #####
###########################################


############### Treatment

##' @export
##' @rdname Design_extractreplace
setGeneric("treatment", function(x) standardGeneric("treatment"))

##' @title Extract and Replace elements of Design
##' @param x Design object
##' @return data.frame containing cluster-level information
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

  value <- .convert_to_data.frame(value, x, "t")

  x@structure[x@columnIndex == "t"] <- value
  names(x@structure)[x@columnIndex == "t"] <- colnames(value)
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
  x@structure[x@columnIndex == "c"]
})

##' @export
##' @rdname Design_extractreplace
setGeneric("clusters<-", function(x, value) standardGeneric("clusters<-"))

##' @export
##' @rdname Design_extractreplace
setMethod("clusters<-", "Design", function(x, value) {

  value <- .convert_to_data.frame(value, x, "c")

  x@structure[x@columnIndex == "c"] <- value
  names(x@structure)[x@columnIndex == "c"] <- colnames(value)

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

  x@structure[x@columnIndex == "b"] <- value
  names(x@structure)[x@columnIndex == "b"] <- colnames(value)
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

  x@structure[x@columnIndex == "f"] <- value
  names(x@structure)[x@columnIndex == "f"] <- colnames(value)
  validObject(x)
  x
})


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
    if (any(dim(design@structure[design@columnIndex == type]) != dim(value))) {
      stop("dimensions of current entires and replacement disagree")
    }

    if (nullName) {
      colnames(value) <- varNames(design, type)
    }
  }
  if (any(dim(design@structure[design@columnIndex == type]) != dim(value))) {
    stop("dimensions of current entires and replacement disagree")
  }
  value
}
