setClass("Design",
         slots = c(structure = "data.frame",
                   columnIndex = "factor",
                   type = "character"))

setValidity("Design", function(object) {
  if (any(dim(object@structure) == 0)) {
    "@structure must have positive dimensions"
  } else if (ncol(object@structure) != length(object@columnIndex)) {
    "@columnIndex does not agree with number of columns in @structure"
  } else {
    #TODO: Add validity for `type`
    TRUE
  }
})

##' These are special function used only in the definition of Design class
##' objects. They identify the clusters, blocks and forcing variables.
##'
##' @title Special terms in Design
##' @param ... any number of variables of the same length.
##' @return the variables with appropriate labels
##' @export
##' @rdname DesignSpecials
forcing <- block <- cluster <- function (...)
{
  #browser()
  allf <- list(...)
  do.call(cbind, allf)
}

##' @rdname DesignSpecials
##' @export
cluster <- block

##' @rdname DesignSpecials
##' @export
forcing <- block

New_Design <- function(form, data, subset = NULL) {
  if (!is.null(subset)) {
    data <- subset(data, subset = subset)
  }

  m <- as.data.frame(as.matrix(model.frame(form, data)))

  index <- factor(rep("t", ncol(m)),
                  levels = c("t", "c", "b", "f"))

  # Handle clusters
  clusters <- grepl("^cluster", colnames(m))
  index[which(clusters)] <- "c"
  cvars <- colnames(m)[clusters][1]
  cvars <- sub("^cluster\\(", "", cvars)
  cvars <- sub("\\)[\\.0-9]*$", "", cvars)
  cvars <- gsub(" ", "", cvars)
  colnames(m)[clusters] <- strsplit(cvars, ",")[[1]]

  # Handle blocks
  blocks <- grepl("^block", colnames(m))
  index[which(blocks)] <- "b"
  bvars <- colnames(m)[blocks][1]
  bvars <- sub("^block\\(", "", bvars)
  bvars <- sub("\\)[\\.0-9]*$", "", bvars)
  bvars <- gsub(" ", "", bvars)
  colnames(m)[blocks] <- strsplit(bvars, ",")[[1]]

  # Handle forcing
  forcings <- grepl("^forcing", colnames(m))
  index[which(forcings)] <- "f"
  fvars <- colnames(m)[forcings][1]
  fvars <- sub("^forcing\\(", "", fvars)
  fvars <- sub("\\)[\\.0-9]*$", "", fvars)
  fvars <- gsub(" ", "", fvars)
  colnames(m)[forcings] <- strsplit(fvars, ",")[[1]]

  new("Design",
      structure = m,
      columnIndex = index,
      type = "None")
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

  design <- New_Design(formula, data, subset = subset)

  design@type <- "RCT"
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

  design <- New_Design(formula, data, subset = subset)

  design@type <- "RD"
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

  design <- New_Design(formula, data, subset = subset)

  design@type <- "Obs"
  return(design)
}

# Perform checks on formula for creation of Design.
# Checks performed:
# - Ensure presence of cluster()
# - Disallow multiple cluster(), block(), or forcing() terms
# - Disallow forcing() unless in RDD
checkDesignFormula <- function(form, allowForcing = FALSE) {
  tt <- terms(form, c("cluster", "block", "forcing"))
  specials <- attr(tt, "specials")

  if (attr(tt, "response") == 0) {
    stop("Must specify a treatment variable as the left side of the formula.")
  }

  if (is.null(specials$cluster)) {
    stop("Must specify at least one clustering variable.")
  } else if (length(specials$cluster) > 1) {
    stop("Specify only one cluster() (cluster() can take multiple variables).")
  }

  if (!is.null(specials$block) && length(specials$block) > 1) {
    stop("Specify only one block() (block() can take multiple variables).")
  }

  if (!allowForcing && !is.null(attr(tt, "specials")$forcing)) {
    stop("forcing() only allowed in RD_Design")
  } else if (allowForcing && length(specials$forcing) > 1) {
    stop("Specify only one forcing() (forcing() can take multiple variables).")
  }

  TRUE
}
