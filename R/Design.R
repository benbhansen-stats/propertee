setClass("Design",
         slots = c(structure = "data.frame",
                   columnIndex = "character",
                   type = "character",
                   unitOfAssignmentType = "character",
                   call = "call"))

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
  if (length(unique(tr)) < 2) {
    return("Invalid treatment; treatment can not be constant")
  }
  if (ncol(object@structure) != length(object@columnIndex)) {
    return("@columnIndex does not agree with number of columns in @structure")
  }
  if (any(colnames(object@structure) != names(object@columnIndex))) {
    return("name disagree between @structure and @columnIndex")
  }
  if (!all(object@columnIndex %in% c("t", "u", "b", "f"))) {
    wrong <- object@columnIndex[!object@columnIndex %in% c("t", "u", "b", "f")]
    return(paste("@columnIndex design elements must be [t,u,b,f]. unknown elements:",
                 paste(wrong, collapse = ", ")))
  }
  if (!object@type %in% c("RCT", "RD", "Obs")) {
    return(paste("@type must be one of [RCT,RD,Obs]. unknown @type:", object@type))
  }
  if (object@type != "RD" && any(object@columnIndex == "f")) {
    return("Forcing variables only valid in RD")
  }
  if (object@type == "RD" && !any(object@columnIndex == "f")) {
    return("RD designs must include at least one forcing variables")
  }
  if (!object@unitOfAssignmentType %in% c("cluster", "unitid", "unitOfAssignment")) {
    return('valid `unitOfAssignmentType`s are "unitOfAssignment", "cluster" or "unitid"')
  }
  TRUE
})


New_Design <- function(form, data, type, subset = NULL, call = NULL) {
  if (is.null(call) | !is.call(call)) {
    call <- match.call()
    warning("Invalid call passed to `New_Design`, using default. Please use RD_Design, RCT_Design, or Obs_Design instead of `New_Design` directly.")
  }

  if (!is.null(subset)) {
    data <- subset(data, subset = subset)
  }

  ### Track whether Design uses uoa/cluster/unitid for nicer output later

  if (grepl("unitOfAssignment\\([a-zA-Z]", deparse(form)) |
        grepl("uoa\\([a-zA-Z]", deparse(form))) {
    autype <- "unitOfAssignment"
  } else if (grepl("cluster\\([a-zA-Z]", deparse(form))) {
    autype <- "cluster"
  } else if (grepl("unitid\\([a-zA-Z]", deparse(form))) {
    autype <- "unitid"
  } else {
    stop("This error should never be hit!")
  }

  # Ensure whichever unit of assignment function is used, `unitOfAssignment` is
  # called
  form <- .updateFormToUnitOfAssignment(form)

  m <- do.call(data.frame, c(model.frame(form, data), check.names = FALSE))

  cd <- .rename_model_frame_columns(m)
  m <- cd[["renamedModelFrame"]]
  index <- cd[["index"]]

  # Ensure there are not variable transformations (e.g. as.factor(x)
  if (!all(names(m) %in% names(data))) {
    stop("Do not use variable transformations in formula.\nInstead modify all relevant data sets as appropriate.")
  }

  m_collapse <- unique(m)

  rownames(m_collapse) <- NULL

  treatment <- m_collapse[, index == "t", drop = FALSE]
  m_collapse[, index == "t"] <- .convert_treatment_to_factor(treatment)

  differing <- duplicated(m_collapse[, index == "u"])
  if (any(differing)) {
    noncon <- m_collapse[differing, index == "u", drop = FALSE]
    cat("\nUnits of assignment with non-constant treatment, block or forcing:\n")
    if (nrow(noncon) >= 6) {
      print(noncon[1:5,,drop = FALSE])
      cat("...\n")
    } else {
      print(noncon)
    }

    stop("Each of treatment assignment, block and forcing must be constant within unit of assignment.")
  }

  new("Design",
      structure = m_collapse,
      columnIndex = index,
      type = type,
      unitOfAssignmentType = autype,
      call = call)
}

##' Generates a Design object with the given specifications.
##'
##' Generates a randomized control treatment Design (`RCT_Design`), or an
##' observational Design (`Obs_Design`), or a regression discontinuity Design
##' (`RD_Design`).
##'
##' The formula must include exactly one of `unitOfAssignment()`, `cluster()`,
##' or `unitid()` to identify the units of assignment (one or more variables).
##' If defining an RD_Design, the formula must also include a `forcing()` entry.
##' The formula may optionally include a `block()` entry as well.
##'
##' Vvariable transformation inside the `*_Design` calls are not allowed, as
##' that can lead to mismatches between the Design objects and future data. In
##' other words, avoid calls such as `RCT_Design(as.factor(z) ~ ...`, and
##' instead transform the variable in ALL relevant data sets a prior.
##' @title Specify Design
##' @param formula defines the design components
##' @param data the data set.
##' @param subset optionally subset the data before creating the design object
##' @return a Design object of the requested type for use in further analysis
##' @export
##' @rdname Design_objects
RCT_Design <- function(formula, data, subset = NULL) {
  .check_design_formula(formula)

  design <- New_Design(form = formula,
                       data = data,
                       type = "RCT",
                       subset = subset,
                       call = match.call())
  return(design)
}

##' @export
##' @rdname Design_objects
RD_Design <- function(formula, data, subset = NULL) {
  .check_design_formula(formula, allowForcing = TRUE)


  design <- New_Design(form = formula,
                       data = data,
                       type = "RD",
                       subset = subset,
                       call = match.call())
  return(design)
}

##' @export
##' @rdname Design_objects
Obs_Design <- function(formula, data, subset = NULL) {
  .check_design_formula(formula)


  design <- New_Design(form = formula,
                       data = data,
                       type = "Obs",
                       subset = subset,
                       call = match.call())
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
  uoatype <- switch(object@unitOfAssignmentType,
                      "unitOfAssignment" = "Unit of Assignment  :",
                      "cluster" = "Cluster  :",
                      "unitid"  = "Unitid   :")

  cat(destype)
  cat("\n\n")
  cat(paste("Treatment:", varNames(object, "t")))
  cat("\n")
  cat(paste(uoatype, paste(varNames(object, "u"), collapse = ", ")))
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
##' @param type one of "t", "u", "b", "f"; for "treatment", "unitOfAssignment",
##'   "block", and "forcing"
##' @return character vector of variable names of the given type
##' @export
varNames <- function(x, type) {
  stopifnot(class(x) == "Design")
  stopifnot(length(type) == 1)
  stopifnot(type %in% c("t", "u", "b", "f"))
  names(x@structure)[x@columnIndex == type]
}

# Internal to properly rename columns to strip unitOfAssignment(), block(), etc
.rename_model_frame_columns <- function(modframe) {

  index <- rep("t", ncol(modframe))

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

  o <- rename_vars(modframe, index, "unitOfAssignment")
  modframe <- o[[1]]
  index <- o[[2]]
  o <- rename_vars(modframe, index, "block")
  modframe <- o[[1]]
  index <- o[[2]]
  o <- rename_vars(modframe, index, "forcing")
  modframe <- o[[1]]
  index <- o[[2]]

  return(list(renamedModelFrame = modframe,
              index = index))
}

##' Returns a table of number of units of assignment in each treatment group,
##' sorted by the size of the groups
##'
##' @title treatment group table
##' @param design A Design object
##' @param ... add'l optional arguments to `table`
##' @return a table of class `treatmentTable`
##' @export
treatmentTable <- function(design, ...) {
  tab <- table(design@structure[varNames(design, "t")], ...)
  tab <- sort(tab, decreasing = TRUE)
  return(tab)
}
