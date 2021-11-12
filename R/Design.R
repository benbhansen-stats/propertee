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
  o <- rename_vars(m, index, "unitid")
  m <- o[[1]]
  index <- o[[2]]
  o <- rename_vars(m, index, "block")
  m <- o[[1]]
  index <- o[[2]]
  o <- rename_vars(m, index, "forcing")
  m <- o[[1]]
  index <- o[[2]]

  if (any(index == "u")) {
    attr(index, "clusterinput") <- "unitid"
    index[index == "u"] <- "c"
  } else {
    attr(index, "clusterinput") <- "cluster"
  }

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
  clusttype <- switch(attr(object@columnIndex, "clusterinput"),
                      "cluster" = "Cluster  :",
                      "unitid"  = "Unitid   :")

  cat(destype)
  cat("\n\n")
  cat(paste("Treatment:", varNames(object, "t")))
  cat("\n")
  cat(paste(clusttype, paste(varNames(object, "c"), collapse = ", ")))
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
