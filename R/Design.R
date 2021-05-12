setClass("Design",
         slots = c(structure = "data.frame",
                   columnIndex = "character",
                   type = "character"))

setValidity("Design", function(object) {
  if (any(dim(object@structure) == 0)) {
    return("@structure must have positive dimensions")
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

  new("Design",
      structure = m,
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
