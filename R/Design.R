setClass("Design",
         slots = c(structure = "data.frame",
                   columnIndex = "character",
                   type = "character"))

setValidity("Design", function(object) {
  if (any(dim(object@structure) == 0)) {
    return("@structure must have positive dimensions")
  }
  if (any(duplicated(colnames(object@structure)))) {
    return("variables cannot be used more than once")
  }
  tr <- object@structure[object@columnIndex == "t"]
  if (ncol(tr) == 0) {
    return("Missing treatment index")
  }
  if (ncol(tr) > 1) {
    return("Only one treatment variable allowed")
  }
  tr <- tr[,1]
  if (is.null(tr) || !is.factor(tr)) {
    return("Invalid treatment; must be factor")
  }
  if (ncol(object@structure) != length(object@columnIndex)) {
    return("@columnIndex does not agree with number of columns in @structure")
  }
  if (any(colnames(object@structure) != names(object@columnIndex))) {
    return("name disagree between @structure and @columnIndex")
  }
  if (!all(object@columnIndex %in% c("t", "c", "b", "f"))) {
    wrong <- object@columnIndex[!object@columnIndex %in% c("t", "c", "b", "f")]
    return(paste("@columnIndex design elements must be [t,c,b,f]. unknown elements:",
                 paste(wrong, collapse = ", ")))
  }
  if (!object@type %in% c("RCT", "RD", "Obs")) {
    return(paste("@type must be one of [RCT,RD,Obs]. unknown @type:", object@type))
  }
  if (object@type != "RD" && any(object@columnIndex == "f")) {
    return("Forcing variables only valid in RD")
  }
  TRUE
})


New_Design <- function(form, data, type, subset = NULL) {
  if (!is.null(subset)) {
    data <- subset(data, subset = subset)
  }

  m <- do.call(data.frame, c(model.frame(form, data), check.names = FALSE))

  index <- rep("t", ncol(m))

  rename_vars <- function(modelframe, index, type) {
    pos <- grepl(paste0("^", type), colnames(modelframe))
    if (any(pos)) {
      index[pos] <- substr(type, 1, 1)
      vars <- colnames(modelframe)[pos][1]
      vars <- sub(paste0("^", type, "\\("), "", vars)
      vars <- sub("\\)[\\.0-9]*$", "", vars)
      vars <- strsplit(gsub(" ", "", vars), ",")[[1]]
      colnames(modelframe)[pos] <- vars
    }
    return(list(modelframe, index))
  }

  o <- rename_vars(m, index, "cluster")
  m <- o[[1]]
  index <- o[[2]]
  o <- rename_vars(m, index, "block")
  m <- o[[1]]
  index <- o[[2]]
  o <- rename_vars(m, index, "forcing")
  m <- o[[1]]
  index <- o[[2]]

  m_collapse <- unique(m)

  rownames(m_collapse) <- NULL

  treatment <- m_collapse[, index == "t", drop = FALSE]
  m_collapse[, index == "t"] <- .convert_treatment_to_factor(treatment)

  differing <- duplicated(m_collapse[, index == "c"])
  if (any(differing)) {
    noncon <- m_collapse[differing, index == "c", drop = FALSE]
    cat("\nClusters with non-constant treatment, block or forcing:\n")
    if (nrow(noncon) >= 6) {
      print(noncon[1:5,,drop = FALSE])
      cat("...\n")
    } else {
      print(noncon)
    }

    stop("Each of treatment assignment, block and forcing must be constant within cluster.")
  }

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

  value <- .convert_treatment_to_factor(value)

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

  x <- .updateStructure(x, value, "c")

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

  .updateStructure(x, value, "b")
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

  .updateStructure(x, value, "f")
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
.updateStructure <- function(design, new, type) {
  design@structure <-
    cbind.data.frame(design@structure[design@columnIndex != type], new)

  design@columnIndex <- c(design@columnIndex[design@columnIndex != type],
                          rep(type, ncol(new)))
  names(design@columnIndex) <- colnames(design@structure)
  validObject(design)
  return(design)
  }

# Internal helper function
# Converts treatment to factor
.convert_treatment_to_factor <- function(treatment) {
  if (!is.null(dim(treatment))) {
    if (ncol(treatment) != 1) {
      stop("Only one treatment variable allowed")
    }
    treatment <- treatment[,,drop = TRUE]
  }
  if (is.numeric(treatment)) {
    if (any(!treatment %in% 0:1)) {
      stop("Numerical treatments must only contain values 0 and 1.")
    }
    treatment <- as.factor(treatment)
  } else if (is.logical(treatment)) {
    treatment <- as.factor(treatment)
  } else if (!is.factor(treatment)) {
    stop("Treatment must be binary (0/1), logical, factor or ordered")
  }
  return(treatment)
}
