setClass("Design",
         slots = c(structure = "data.frame",
                   column_index = "character",
                   type = "character",
                   unit_of_assignment_type = "character",
                   call = "call",
                   dichotomy = "formula"))

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
  if (is.null(tr) ||
        (!is.numeric(tr) && !is.character(tr)) && !is.logical(tr)) {
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
    wrong <- object@column_index[!object@column_index %in%
                                   c("t", "u", "b", "f")]
    return(paste("@column_index design elements must be [t,u,b,f].",
                 "unknown elements:", paste(wrong, collapse = ", ")))
  }
  if (!object@type %in% c("RCT", "RD", "Obs")) {
    return(paste("@type must be one of [RCT,RD,Obs]. unknown @type:",
                 object@type))
  }
  if (object@type != "RD" && any(object@column_index == "f")) {
    return("Forcing variables only valid in RD")
  }
  if (object@type == "RD" && !any(object@column_index == "f")) {
    return("RD designs must include at least one forcing variables")
  }
  if (!object@unit_of_assignment_type %in%
        c("cluster", "unitid", "unit_of_assignment")) {
    return(paste('valid `unit_of_assignment_type`s are "unit_of_assignment",',
                 '"cluster" or "unitid"'))
  }
  if (!length(object@dichotomy) %in% c(0, 3)) {
    return("@dichotomy invalid")
  }
  TRUE
})


##' @title Create a new \code{Design} object.
##' @param form Formula to create Design, see help for \code{rcr_design()},
##'   \code{rd_design()} or \code{obs_design()} for details.
##' @param data The data set
##' @param type One of "RCT", "RD", or "Obs"
##' @param subset Any subset information
##' @param call The call generating the \code{Design}.
##' @param dichotomy If present, the dichotomization formula
##' @return A new Design object
##' @importFrom stats formula
##' @keywords internal
new_Design <- function(form,
                       data,
                       type,
                       subset = NULL,
                       call = NULL,
                       dichotomy = stats::formula()) {

  if (is.null(call) | !is.call(call)) {
    call <- match.call()
    warning(paste("Invalid call passed to `new_Design`, using default.",
                  "Please use rd_design, rct_design, or obs_design instead ",
                  "of `new_Design` directly."))
  }

  if (!is.null(subset)) {
    data <- subset(data, subset = subset)
  }

  # `formula()` is a close equivalent of a NULL formula
  if (is.null(dichotomy)) dichotomy <- stats::formula()
  # the fact that a formula has an environment is playing hell with testing. I
  # don't believe we'll ever need the environment in which the dichotomy
  # formula is created as we use it on wahtever data we need, so setting it
  # always to the generic interactive environment for simplicity.
  environment(dichotomy) <- globalenv()

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

  # Ensure whichever unit of assignment function is used, `unit_of_assignment`
  # is called
  form <- .update_form_to_unit_of_assignment(form)

  m <- do.call(data.frame,
               c(model.frame(form, data, na.action = na.pass),
                 check.names = FALSE))

  cd <- .rename_model_frame_columns(m)
  m <- cd[["renamedModelFrame"]]
  index <- cd[["index"]]

  m_collapse <- unique(m)

  rownames(m_collapse) <- NULL

  ########### Examine treatment variable.
  treatment <- m_collapse[, index == "t"]

  if (options()$flexida_warn_on_conditional_treatment &
                grepl("[<>=]", deparse(form[[2]]))) {
    # Case 1: LHS of form has a conditional (e.g. (dose > 50). Search for >, <,
    # = and if found, evaluate and convert to numeric, but warn users.
    warning(paste("It appears that you've identified the treatment group",
                  "using conditional logic\n(by including one of <, >, or =",
                  "in the left hand side of `form`).\n",
                  "This is supported, but it is recommended to instead",
                  "include the non-binary\ntreatment variable in the `form`",
                  "and use the `dichotomy` to define the groups.\n",
                  "Using `dichotomy` will make modifying the groups",
                  "easier in the future\nshould you need to adjust."))
    # If the user is using conditionals, we'll be converting logical to numeric
    # later but don't need to the add'l warning message.
    if (!is.logical(treatment)) {
      stop(paste("treatment has conditional but isn't logical.\n",
                 "This may happen if treatment variable name has\n",
                 "'<', '>', or '=' in it. Rename variable to proceed."))
    }
  }
  else if (is.factor(treatment)) {
    # Case 2: Input is a factor or ordinal. Warn user before converting to
    # numeric if levels are numeric, otherwise return as character.
    warning(paste("Factor treatment variables have their labels converted",
                  "to numeric/character as appropriate.\nIt is suggested",
                  "to do this conversion yourself to ensure it proceeds",
                  "as you expect."))
    newt <- levels(treatment)[treatment]
    # if levels were numeric, convert to numeric. Otherwise leave as is
    treatment <- tryCatch(as.numeric(newt),
                     warning = function(w) {
                       newt
                     }, error = function(e) {
                       newt # should never hit here (warning only)
                     })
  } else if (!is.numeric(treatment) &
             !is.character(treatment) &
             !is.logical(treatment)) {
    # Case 3: Treatment is not a factor, a conditional, a numeric, a logical or
    # a character. It is some other type. Warn user before converting to
    # numeric.
    warning(paste("Treatment variables which are not numeric or character are",
                  "converted into numeric.\nIt is STRONGLY suggested to do",
                  "this conversion yourself to ensure it proceeds",
                  "as you expect."))
    treatment <- as.numeric(treatment)
  }
  # Case 0: treatment is numeric, logical, or character. No additional steps
  # needed.

  m_collapse[, index == "t"] <- treatment

  differing <- duplicated(m_collapse[, index == "u"])
  if (any(differing)) {
    noncon <- m_collapse[differing, index == "u", drop = FALSE]
    cat(paste("\nUnits of assignment with non-constant treatment, block",
              "or forcing:\n"))
    if (nrow(noncon) >= 6) {
      print(noncon[1:5, , drop = FALSE])
      cat("...\n")
    } else {
      print(noncon)
    }

    stop(paste("Each of treatment assignment, block and forcing must be",
               "constant within unit of assignment."))
  }

  new("Design",
      structure = m_collapse,
      column_index = index,
      type = type,
      unit_of_assignment_type = autype,
      call = call,
      dichotomy = dichotomy)
}

##' Generates a randomized control treatment Design (\code{rct_design()}), or an
##' observational Design (\code{obs_design()}), or a regression discontinuity
##' Design (\code{rd_design()}).
##'
##' The formula must include exactly one of \code{unit_of_assignment()},
##' \code{uoa()}, \code{cluster()}, or \code{unitid()} to identify the units of
##' assignment (one or more variables). If defining an rd_design, the formula
##' must also include a \code{forcing()} entry. The formula may optionally
##' include a \code{block()} entry as well.
##'
##' The treatment variable passed into the left-hand side of \code{formula} can
##' either be \code{logical}, \code{numeric}, or \code{character}. If it is
##' anything else, it is attemped to be converted to one of those (for example,
##' \code{factor} and \code{ordered} are converted to \code{numeric} if the
##' levels are \code{numeric}, otherwise to \code{character}. If the treatment
##' is not \code{logical} or \code{numeric} with only values 0 and 1, in order
##' to generate weights with \code{ate()} or \code{ett()}, the \code{dichotomy}
##' argument must be used to identify the treatment and control groups. The
##' \code{Design} creation functions (\code{rct_design()}, \code{rd_design()},
##' \code{obs_design()}) all support the \code{dichotomy} argument, or instead
##' \code{dichotomy} can be passed to \code{ett()} and \code{ate()} directly.
##'
##' The \code{dichotomy} argument should be a formula consisting of a
##' conditional statement on both the left-hand side (identifying treatment
##' levels associated with "treatment") and the right hand side (identifying
##' treatment levels associated with "control"). For example, if your treatment
##' variable was called \code{dose}, you might write:
##'
##' \code{dichotomy = dose > 250 ~ dose <= 250}
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
##' \code{dichotomy = dose > 250 ~ .}
##'
##' or
##'
##' \code{dichotomy = . ~ dose <= 250}
##'
##' The \code{dichotomy} formula supports all Relational Operators, Logical
##' Operators, and \code{%in%}.
##'
##' Note that you can specify a conditional logic treatment in the formula (e.g.
##' \code{rct_design(dose > 250 ~ unitOfAssignment(...}) but we would suggest
##' instead passing the treatment variable directly and using \code{dichotomy}.
##' Otherwise changing the dichotomy will require re-creating the Design,
##' instead of simply using \code{dichotomy(design) <-} or passing
##' \code{dichotomy} to \code{ate()} or \code{ett()}.
##'
##' @title Generates a \code{Design} object with the given specifications.
##' @param formula defines the \code{Design} components
##' @param data the data set.
##' @param subset optionally subset the data before creating the \code{Design}
##'   object
##' @param dichotomy optionally, a formula defining the dichotomy of the
##'   treatment variable if it isn't already 0/1 or \code{logical}. See details.
##' @return a Design object of the requested type for use in further analysis
##' @export
##' @rdname Design_objects
rct_design <- function(formula, data, subset = NULL, dichotomy = NULL) {
  .check_design_formula(formula)

  new_Design(form = formula,
             data = data,
             type = "RCT",
             subset = subset,
             call = match.call(),
             dichotomy = dichotomy)
}

##' @export
##' @rdname Design_objects
rd_design <- function(formula, data, subset = NULL, dichotomy = NULL) {
  .check_design_formula(formula, allow_forcing = TRUE)

  new_Design(form = formula,
             data = data,
             type = "RD",
             subset = subset,
             call = match.call(),
             dichotomy = dichotomy)
}

##' @export
##' @rdname Design_objects
obs_design <- function(formula, data, subset = NULL, dichotomy = NULL) {
  .check_design_formula(formula)

  new_Design(form = formula,
             data = data,
             type = "Obs",
             subset = subset,
             call = match.call(),
             dichotomy = dichotomy)
}

##' @title Show a Design
##' @param object Design object
##' @return an invisible copy of \code{object}
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
    cat(paste("Dichotomy rule:", deparse(object@dichotomy)))
    cat("\n")
  }

  cat("\n")
  invisible(object)
})

##' @title Extract names of Design variables
##' @param x Design x
##' @param type one of "t", "u", "b", "f"; for "treatment",
##'   "unit_of_assignment", "block", and "forcing"
##' @return character vector of variable names of the given type
##' @export
var_names <- function(x, type) {
  stopifnot(class(x) == "Design")
  stopifnot(length(type) == 1)
  stopifnot(type %in% c("t", "u", "b", "f"))
  names(x@structure)[x@column_index == type]
}

# Internal to properly rename columns to strip unit_of_assignment(), block(),
# etc
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
##' @param ... add'l optional arguments to \code{table}
##' @return a table of treatment by units
##' @export
treatment_table <- function(design, ...) {
  tab <- table(design@structure[var_names(design, "t")], ...)
  tab <- sort(tab, decreasing = TRUE)
  return(tab)
}

##' Returns a table containing the variables identified in each structure
##'
##' When \code{compress} is \code{TRUE}, the result will always have two
##' columns. When \code{FALSE}, the result will have number of columns equal to
##' the largest number of variables in a particular role, plus one. E.g., a call
##' such as \code{rct_design(z ~ unitid(a, b, c, d) ...} will have 4+1=5 columns
##' in the output matrix.
##'
##' When \code{report_all} is \code{TRUE}, the matrix is guaranteed to have 3
##' rows (if the \code{design} is an RCT or Obs) or 4 rows (when the
##' \code{design} is a RD). When \code{FALSE}, the matrix will have minimum 2
##' rows (treatment and cluster/unitid/unif of assignment), with additional rows
##' for blocks and forcing if included in the \code{Design}.
##' @title variable identification table
##' @param design A Design object
##' @param compress Should multiple variables be compressed into a
##'   comma-separated string? Default \code{TRUE}.
##' @param report_all Should we report all possible structures even if they
##'   don't exist in the Design? Default \code{FALSE}.
##' @return a \code{matrix} of variables in the Design structure
##' @export
var_table <- function(design, compress = TRUE, report_all = FALSE) {
  uoatype <- switch(design@unit_of_assignment_type,
                    "unit_of_assignment" = "Unit of Assignment",
                    "cluster" = "Cluster",
                    "unitid"  = "Unitid")

  # Start with the "table" as a list; easier to handle non-equal number of
  # elements.
  rows <- list()
  rows[["t"]] <- c("Treatment", var_names(design, "t"))
  rows[["u"]] <- c(uoatype    , var_names(design, "u"))
  rows[["b"]] <- c("Block"    , var_names(design, "b"))
  if (design@type == "RD") {
    rows[["f"]] <- c("Forcing"  , var_names(design, "f"))
  }

  # Identify if we have more than one variable specified in a specific
  # structure, and if so how ,long
  maxvar <- max(vapply(rows, length, 1))
  if (maxvar < 2) {
    # Should never hit this error
    stop("Internal error: No variables identified!")
  }

  if (!report_all) {
    # Drop any rows we don't need
    rows <- rows[vapply(rows, length, 1) > 1]
  }

  # add NA padding to make all rows the same length so the `rbind` that follows
  # doesn't need to recycle length
  rows <- lapply(rows, function(x) {
    if (length(x) < maxvar) {
      x <- c(x, rep(NA, maxvar - length(x)))
    }
    return(x)
  })

  # list to data.frame
  out <- do.call(rbind, rows)

  if (compress == TRUE) {
    # If compressing, collapse any repeats if necessary
    rows <- lapply(rows, function(x) {
      x <- x[!is.na(x)]
      c(x[1], paste(x[-1], collapse = ", "))
    })
    out <- do.call(rbind, rows)
    colnames(out) <- c("Structure", "Variables")
    out[out == ""] <- NA
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

##' Useful for debugging purposes to ensure that there is concordance between
##' variables in the \code{Design} and data.
##'
##' Consider the following scenario: A \code{Design} is generated from some
##' dataset, "data1", which includes a block variable "b1". Within each unique
##' cluster/unit of assignment of "data1", it must be the case that "b1" is
##' constant. (Otherwise the creation of the \code{Design} will fail.)
##'
##' Next, a model is fit which includes weights generated from the
##' \code{Design}, but on dataset "data2". In "data2", the block variable "b1"
##' also exists, but due to some issue with data cleaning, does not agree with
##' "b1" in "data1".
##'
##' This could cause errors, either directly (via actual error messages) or
##' simply produce nonsense results. `design_data_concordance()` is designed to
##' help debug these scenarios but providing information on whether variables in
##' both "data1" (used in the creation of the \code{Design}) and "data2" (some
##' new dataset passed to \code{design_data_concordance()}) have any
##' inconsistencies.
##' @title Checks for variable agreement within units of assignment
##' @param design A Design object
##' @param data Data set, presumably not the same used to create \code{design}.
##' @param by optional; named vector or list connecting names of variables in
##'   \code{design} to variables in \code{data}. Names represent variables in
##'   \code{design}; values represent variables in \code{data}. Only needed if
##'   variable names differ.
##' @param warn_on_nonexistence Default \code{TRUE}. If a variable does not
##'   exist in \code{data}, should this be flagged? If \code{FALSE}, silently
##'   move on if a variable doesn't exist in \code{data}.
##' @return Invisibly \code{TRUE} if no warnings are produced, \code{FALSE} if
##'   any warnings are produced.
##' @export
design_data_concordance <- function(design,
                                    data,
                                    by = NULL,
                                    warn_on_nonexistence = TRUE) {
  if (!is.null(by)) {
    design <- .update_by(design, data, by)
  }

  if (!all(var_names(design, "u") %in% names(data))) {
    stop(paste("Missing cluster/unit of assignment variables in data set.\n",
               "Are you missing a `by=` argument?"))
  }

  merged <- merge(design@structure, data,
                  all = TRUE,
                  by = var_names(design, "u"))

  # Returns TRUE if any warnings were printed, FALSE otherwise
  .check <- function(type, design, data) {
    anywarnings <- FALSE

    vnames <- var_names(design, type)
    # If this particular type of variable exists in the design...
    if (length(vnames) > 0) {
      # Loop over all variables of `type`
      for (var in vnames) {
        # If variable name exists in the data ...
        xynames <- paste0(var, c(".x", ".y"))
        if (all(xynames %in% names(data))) {
          varcompare <- merged[, xynames]
          # Ensure differences are 0
          if (!all(apply(varcompare[, xynames], 1, diff) == 0)) {
            warning(paste0("Inconsistencies in variable `", var, "`"))
            anywarnings <- TRUE
          }
        } else {
          if (warn_on_nonexistence) {
            warning(paste0("Variable `", var, "` not found in data"))
            anywarnings <- TRUE
          }
        }
      }
    }
    if (anywarnings) {
      return(invisible(TRUE))
    } else {
      return(invisible(FALSE))
    }
  }

  treat_had_warnings <- .check("t", design, merged)
  block_had_warnings <- .check("b", design, merged)
  force_had_warnings <- .check("f", design, merged)

  if (any(treat_had_warnings,
          block_had_warnings,
          force_had_warnings)) {
    return(invisible(FALSE))
  } else {
    return(invisible(TRUE))
  }
}
