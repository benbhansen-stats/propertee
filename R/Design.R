Design <- setClass("Design",
                   slots = c(structure = "data.frame",
                             column_index = "character",
                             type = "character",
                             unit_of_assignment_type = "character",
                             call = "call",
                             treatment_binary = "list"))

## If user passes a dichotomize argument:
## treatment_binary[[1]] is formula
## treatment_binary[[2]] is original treatment

setValidity("Design", function(object) {
  if (any(dim(object@structure) == 0)) {
    return("@structure must have positive dimensions")
  }
  if (any(duplicated(colnames(object@structure)))) {
    return("variables cannot be used more than once")
  }
  tr <- object@structure[object@column_index == "t"]
  if (ncol(tr) == 0) {
    return("Missing treatment index")
  }
  if (ncol(tr) > 1) {
    return("Only one treatment variable allowed")
  }
  tr <- tr[, 1]
  if (is.null(tr) || (!is.numeric(tr) && !is.character(tr))) {
    return("Invalid treatment; must be numberic or character")
  }
  if (length(unique(tr)) < 2) {
    return("Invalid treatment; treatment can not be constant")
  }
  if (ncol(object@structure) != length(object@column_index)) {
    return("@column_index does not agree with number of columns in @structure")
  }
  if (any(colnames(object@structure) != names(object@column_index))) {
    return("name disagree between @structure and @column_index")
  }
  if (!all(object@column_index %in% c("t", "u", "b", "f"))) {
    wrong <- object@column_index[!object@column_index %in% c("t", "u", "b", "f")]
    return(paste("@column_index design elements must be [t,u,b,f]. unknown elements:",
                 paste(wrong, collapse = ", ")))
  }
  if (!object@type %in% c("RCT", "RD", "Obs")) {
    return(paste("@type must be one of [RCT,RD,Obs]. unknown @type:", object@type))
  }
  if (object@type != "RD" && any(object@column_index == "f")) {
    return("Forcing variables only valid in RD")
  }
  if (object@type == "RD" && !any(object@column_index == "f")) {
    return("RD designs must include at least one forcing variables")
  }
  if (!object@unit_of_assignment_type %in% c("cluster", "unitid", "unit_of_assignment")) {
    return('valid `unit_of_assignment_type`s are "unit_of_assignment", "cluster" or "unitid"')
  }
  if (!length(object@treatment_binary) %in% c(0, 2)) {
    return("@treatment_binary wrong length")
  }
  if (length(object@treatment_binary) == 2) {
    if (!is(object@treatment_binary[[1]], "formula")) {
      return("@treatment_binary[[1]] must be formula")
    }
    if (!is.data.frame(object@treatment_binary[[2]])) {
      return("@treatment_binary[[1]] must be vector")
    }
    if (nrow(object@treatment_binary[[2]]) != nrow(object@structure)) {
      return("@treatment_binary[[1]] and @structure must have same number of elements")
    }
  }
  TRUE
})


new_Design <- function(form, data, type, subset = NULL, call = NULL, dichotomize = NULL) {
  if (is.null(call) | !is.call(call)) {
    call <- match.call()
    warning("Invalid call passed to `new_Design`, using default. Please use rd_design, rct_design, or obs_design instead of `new_Design` directly.")
  }

  if (!is.null(subset)) {
    data <- subset(data, subset = subset)
  }

  ### Track whether Design uses uoa/cluster/unitid for nicer output later

  if (grepl("unit_of_assignment\\([a-zA-Z]", deparse(form)) |
        grepl("uoa\\([a-zA-Z]", deparse(form))) {
    autype <- "unit_of_assignment"
  } else if (grepl("cluster\\([a-zA-Z]", deparse(form))) {
    autype <- "cluster"
  } else if (grepl("unitid\\([a-zA-Z]", deparse(form))) {
    autype <- "unitid"
  } else {
    stop("This error should never be hit!")
  }

  # Ensure whichever unit of assignment function is used, `unit_of_assignment` is
  # called
  form <- .update_form_to_unit_of_assignment(form)

  m <- do.call(data.frame, c(model.frame(form, data), check.names = FALSE))

  cd <- .rename_model_frame_columns(m)
  m <- cd[["renamedModelFrame"]]
  index <- cd[["index"]]

  m_collapse <- unique(m)

  rownames(m_collapse) <- NULL

  # If there are any non-numeric, non-string treatment columns, convert to numeric
  # with a warning.
  treatment <- m_collapse[, index == "t", drop = FALSE]
  character_treatment <- vapply(treatment, is.character, TRUE)
  if (any(!character_treatment)) {
    # idenify non-character treatments which aren't already numeric
    non_numeric_treatment <- !character_treatment & !vapply(treatment, is.numeric, TRUE)
    if (any(non_numeric_treatment)) {
      warning(paste("Treatment is not numeric or character. Attempting to convert",
                    "to numeric.\n It is suggested to convert to numeric or",
                    "character prior to creating Design."))
      treatment[, non_numeric_treatment] <- lapply(treatment[, non_numeric_treatment, drop = FALSE],
                                                    as.numeric)
    }
    # TODO: HAndle factors/ordered more gracefully
  }


  if (!is.null(dichotomize)) {
    if ("z__" %in% names(data)) {
      # We name the dichotomized treatment variable "z__", this is just a
      # warning to users if they happen to have a column in the data of the same
      # name.
      warning(paste('Variable "z__" found in data. Due to dichotomization,',
                    '"z__" is used as the binary treatment variable.\n',
                    'This isn\' necessarily a problem, but be careful!'))
    }

    # First, store the dichotomization info
    treatment_binary <- list(dichotomize = dichotomize,
                             original_treatment = treatment)

    # Replace treatment with dichotomized version with name "z__"
    treatment <- data.frame(z__ = .binarize_treatment(treatment, dichotomize))

    # Replace treatment in `m_collapse`. Future-proofing for multiple-treatment
    # variables, update index as well.
    m_collapse <- cbind(treatment, m_collapse[index != "t"])
    index <- c("t", index[index != "t"])
  } else {
    # If no dichotomization, `@treatment_binary` is empty list.
    treatment_binary <- list()
    m_collapse[, index == "t"] <- treatment
  }


  differing <- duplicated(m_collapse[, index == "u"])
  if (any(differing)) {
    noncon <- m_collapse[differing, index == "u", drop = FALSE]
    cat("\nUnits of assignment with non-constant treatment, block or forcing:\n")
    if (nrow(noncon) >= 6) {
      print(noncon[1:5, , drop = FALSE])
      cat("...\n")
    } else {
      print(noncon)
    }

    stop("Each of treatment assignment, block and forcing must be constant within unit of assignment.")
  }

  new("Design",
      structure = m_collapse,
      column_index = index,
      type = type,
      unit_of_assignment_type = autype,
      call = call,
      treatment_binary = treatment_binary)
}

##' Generates a Design object with the given specifications.
##'
##' Generates a randomized control treatment Design (`rct_design`), or an
##' observational Design (`obs_design`), or a regression discontinuity Design
##' (`rd_design`).
##'
##' The formula must include exactly one of `unit_of_assignment()`, `uoa()`,
##' `cluster()`, or `unitid()` to identify the units of assignment (one or more
##' variables). If defining an rd_design, the formula must also include a
##' `forcing()` entry. The formula may optionally include a `block()` entry as
##' well.
##'
##' The treatment variable passsed into \code{formula} can be either a character
##' or an object that can be coerced into \code{numeric} (e.g. \code{factor} or
##' \code{logical}). If the treatment is not a \code{\numeric}
##' \code{0}/\code{1}, then in order to generate weights with \code{ate()} or
##' \code{ett()}, the \code{dichotomize} argument must be used to identify the
##' treatment and control groups. The \code{Design} creation functions
##' (\code{rct_design()}, \code{rd_design()}, \code{obs_design()}) all support
##' the \code{dichotomize} argument, or instead \code{dichotomize} can be passed
##' to \code{ett()} and \code{ate()} directly.
##'
##' The \code{dichotomize} argument should be a formula consisting of a
##' conditional statement on both the left-hand side (identifying treatment
##' levels associated with "treatment") and the right hand side (identifying
##' treatment levels associated with "control"). For example, if your treatment
##' variable was called \code{dose}, you might write:
##'
##' \code{dichotomize = dose > 250 ~ dose <= 250}
##'
##' The conditionals need not assign all values of treatment to control or
##' treatment, for example, \code{dose > 300 ~ dose < 200} does not handle
##' \code{200 <= dose <= 300}. Units of assignment with treatment values not
##' assigned to either treatment or control as assumed to not recieve either,
##' but will be maintained in the Design for proper standard error calculations.
##'
##' The period (\code{.}) can be used to assign all other units of assignment.
##' For example, we could have written the above example as
##'
##' \code{dichotomize = dose > 250 ~ .}
##'
##' or
##'
##' \code{dichotomize = . ~ dose <= 250}
##'
##' The \code{dichotomize} formula supports all Relational Operators, Logical
##' Operators, and \code{%in%}.
##'
##' @title Specify Design
##' @param formula defines the Design components
##' @param data the data set.
##' @param subset optionally subset the data before creating the Design object
##' @param dichotomize optionally, a formula defining the dichotomization of the
##'   treatment variable if it isn't already \code{0}/\code{1}. See details.
##' @return a Design object of the requested type for use in further analysis
##' @export
##' @rdname Design_objects
rct_design <- function(formula, data, subset = NULL, dichotomize = NULL) {
  .check_design_formula(formula)

  new_Design(form = formula,
             data = data,
             type = "RCT",
             subset = subset,
             call = match.call(),
             dichotomize = dichotomize)
}

##' @export
##' @rdname Design_objects
rd_design <- function(formula, data, subset = NULL, dichotomize = NULL) {
  .check_design_formula(formula, allow_forcing = TRUE)

  new_Design(form = formula,
             data = data,
             type = "RD",
             subset = subset,
             call = match.call(),
             dichotomize = dichotomize)
}

##' @export
##' @rdname Design_objects
obs_design <- function(formula, data, subset = NULL, dichotomize = NULL) {
  .check_design_formula(formula)

  new_Design(form = formula,
             data = data,
             type = "Obs",
             subset = subset,
             call = match.call(),
             dichotomize = dichotomize)
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
  uoatype <- switch(object@unit_of_assignment_type,
                    "unit_of_assignment" = "Unit of Assignment",
                    "cluster" = "Cluster",
                    "unitid"  = "Unitid")

  cat(destype)
  cat("\n\n")


  vartab <- var_table(object)
  # Add a nice separating line for printing
  vartab <- rbind(c("---------", "---------"),
                  vartab)
  print(data.frame(vartab), row.names = FALSE, right = FALSE)

  if (is_dichotomized(object)) {
    cat("\n")
    cat(paste("Treatment was dichotomized; original treatment variable:",
              names(object@treatment_binary[[2]])))
    cat("\n")
  }

  cat("\n")
  invisible(object)
})

##' @title Extract names of Design variables
##' @param x Design x
##' @param type one of "t", "u", "b", "f"; for "treatment", "unit_of_assignment",
##'   "block", and "forcing"
##' @return character vector of variable names of the given type
##' @export
var_names <- function(x, type) {
  stopifnot(class(x) == "Design")
  stopifnot(length(type) == 1)
  stopifnot(type %in% c("t", "u", "b", "f"))
  names(x@structure)[x@column_index == type]
}

# Internal to properly rename columns to strip unit_of_assignment(), block(), etc
.rename_model_frame_columns <- function(modframe) {

  index <- rep("t", ncol(modframe))

  rename_vars <- function(modelframe, index, type) {
    pos <- grepl(paste0("^", type), colnames(modelframe))
    if (any(pos)) {
      index[pos] <- substr(type, 1, 1)
      vars <- colnames(modelframe)[pos][1]
      vars <- sub(paste0("^", type, "\\("), "", vars)
      if (grepl("^cbind", vars)) {
        # If user called something like `uoa(cbind(a, b))`
        vars <- sub("^cbind\\(", "", vars)
        vars <- sub("\\)\\).*$", "", vars)
      } else {
        # If user just called `uoa(a, b)`
        vars <- sub("\\)[\\.0-9]*$", "", vars)
      }
      vars <- strsplit(gsub(" ", "", vars), ",")[[1]]
      colnames(modelframe)[pos] <- vars
    }
    return(list(modelframe, index))
  }

  o <- rename_vars(modframe, index, "unit_of_assignment")
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
##' @return a table of treatment by units
##' @export
treatment_table <- function(design, ...) {
  tab <- table(design@structure[var_names(design, "t")], ...)
  tab <- sort(tab, decreasing = TRUE)
  return(tab)
}

##' Returns a table containing the variables identified in each structure
##'
##' @title variable identification table
##' @param design A Design object
##' @param ... add'l optional arguments to `table`
##' @param compress Should multiple variables be compressed into a
##'   comma-separated string? Default TRUE.
##' @param report_all Should we report all possible structures even if they don't
##'   exist in the Design? Default FALSE.
##' @return a table of variables in the Design structure
##' @export
var_table <- function(design, ..., compress = TRUE, report_all = FALSE) {
  uoatype <- switch(design@unit_of_assignment_type,
                    "unit_of_assignment" = "Unit of Assignment",
                    "cluster" = "Cluster",
                    "unitid"  = "Unitid")

  rows <- list()
  rows[["t"]] <- c("Treatment", var_names(design, "t"))
  rows[["u"]] <- c(uoatype    , var_names(design, "u"))
  rows[["b"]] <- c("Block"    , var_names(design, "b"))
  rows[["f"]] <- c("Forcing"  , var_names(design, "f"))

  # Identify if we have more than one variable specified in a specific
  # structure, and if so how ,long
  maxvar <- max(vapply(rows, length, 1))
  if (maxvar < 2) {
    # Should never hit this error
    stop("Internal error: No variables identified!")
  }

  # If we have more than 1 variable in at least 1 structure, add NA padding to
  # each row
  if (maxvar > 2) {
    rows <- lapply(rows, function(x) {
      if (length(x) < maxvar) {
        x <- c(x, rep(NA, maxvar - length(x)))
      }
      return(x)
    })
  }
  out <- do.call(rbind, rows)

  if (!report_all) {
    # Drop any rows we don't need
    out <- out[!is.na(out[, 2]), ]
  }

  if (compress == TRUE) {
    # If compressing, collapse any repeats if necessary
    if (maxvar > 2) {
      out2 <- apply(out[, -1], 1, function(x) {
        x <- x[!is.na(x)]
        return(paste(x, collapse = ", "))
      })
      out <- cbind(out[, 1], out2)
    }
    colnames(out) <- c("Structure", "Variables")
  } else {
    if (maxvar == 2) {
      colnames(out) <- c("Structure", "Variables")
    } else {
      colnames(out) <- c("Structure", paste("Variable", seq_len(maxvar - 1)))
    }
  }
  rownames(out) <- NULL
  return(out)
}
