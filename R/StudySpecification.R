setClass("StudySpecification",
         slots = c(structure = "data.frame",
                   column_index = "character",
                   type = "character",
                   unit_of_assignment_type = "character",
                   call = "call"))

setValidity("StudySpecification", function(object) {
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
        (!is.factor(tr) && !is.numeric(tr) &&
         !is.character(tr)) && !is.logical(tr)) {
    return("Invalid treatment; must be factor, numeric or character")
  }
  if (length(table(tr)) < 2) {
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
    return(paste("@column_index specification elements must be [t,u,b,f].",
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
    return("RD specifications must include at least one forcing variables")
  }
  if (!object@unit_of_assignment_type %in%
        c("cluster", "unitid", "unit_of_assignment", "none")) {
    return(paste('valid `unit_of_assignment_type`s are "unit_of_assignment",',
                 '"cluster" or "unitid"'))
  }
  return(TRUE)
})


##' Helper function to create a new \code{StudySpecification}. Called internally
##' from \code{rct_spec()}, \code{rd_spec()} or \code{obs_spec()}.
##' @title (Internal) Create a new \code{StudySpecification} object.
##' @param form Formula to create StudySpecification, see help for
##'   \code{rcr_spec()}, \code{rd_spec()} or \code{obs_spec()} for details.
##' @param data The data set
##' @param type One of "RCT", "RD", or "Obs"
##' @param subset Any subset information
##' @param call The call generating the \code{StudySpecification}.
##' @param na.fail Should it error on NA's (\code{TRUE}) or remove them
##'   (\code{FALSE})?
##' @param called_from_lmitt Logical; was this called inside \code{lmitt()}, or
##'   was it called from \code{*_spec()} (default).
##' @return A new StudySpecification object
##' @importFrom stats formula complete.cases terms
##' @importFrom utils capture.output
##' @keywords internal
.new_StudySpecification <- function(form,
                       data,
                       type,
                       subset = NULL,
                       call = NULL,
                       na.fail = TRUE,
                       called_from_lmitt = FALSE) {

  if (is.null(call) | !is.call(call)) {
    call <- match.call()
    warning(paste("Invalid call passed to `.new_StudySpecification`, using default.",
                  "Please use rd_spec, rct_spec, or obs_spec instead ",
                  "of `.new_StudySpecification` directly."))
  }

  if (!is.null(subset)) {
    data <- subset(data, subset = subset)
  }

  ## #174 convert all data.frames
  datadf <- .as_data_frame(data)
  # Moved this prior to `environment` call to hopefully detect more
  # bad input.

  ## keep formula's environment
  env <- environment(terms(form, data = data))
  environment(form) <- env
  call$formula <- form

  data <- datadf


  ### Track whether StudySpecification uses uoa/cluster/unitid for nicer output
  ### later

  if (grepl("unit_of_assignment\\([a-zA-Z]", deparse1(form)) |
        grepl("uoa\\([a-zA-Z]", deparse1(form))) {
    autype <- "unit_of_assignment"
  } else if (grepl("cluster\\([a-zA-Z]", deparse1(form))) {
    autype <- "cluster"
  } else if (grepl("unitid\\([a-zA-Z]", deparse1(form))) {
    autype <- "unitid"
  } else {
    autype <- "none"
    if (options()$propertee_warn_on_no_unit_of_assignment &
                  !called_from_lmitt) {
      warning(paste("The StudySpecification was created without an explicit",
                    "unit of assignment/unit ID. Merges going forward will",
                    "be done by row. It is up to the user to ensure that",
                    "row order is not modified.\nTo prevent this warning,",
                    "provide an explicit unit of assignment/unit ID."))
    }

    data[["..uoa.."]] <- rownames(data)
    form <- update(form, . ~ . + unit_of_assignment(..uoa..))
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

  # #94 handling NA's in non-treatment columns
  completecases <- stats::complete.cases(m[, index != "t"])
  na_tx  <- is.na(m[, index == "t", drop=TRUE])
  if (!all(completecases)) {
    if ( na.fail & !all(completecases | na_tx) ) {
      stop(paste("Missing values cannot be found in unit of assignment,",
                 "block or cluster variables (unless treatment is also NA).",
                 "Use option `na.fail = FALSE` for automatic removal of",
                 "incomplete cases."))
    } else {
      m <- m[completecases, ]
    }
  }





  m_collapse <- unique(m)

  rownames(m_collapse) <- NULL

  ########### Examine treatment variable.
  treatment <- m_collapse[, index == "t"]

  if (options()$propertee_warn_on_conditional_treatment &
                grepl("[<>=]", deparse1(form[[2]]))) {
    # If the user is using conditionals, we'll be converting logical to numeric
    # later but don't need to the add'l warning message.
    if (!is.logical(treatment)) {
      stop(paste("treatment has conditional but isn't logical.\n",
                 "This may happen if treatment variable name has\n",
                 "'<', '>', or '=' in it. Rename variable to proceed."))
    }
  } else if (!is.factor(treatment) &
             !is.numeric(treatment) &
             !is.character(treatment) &
             !is.logical(treatment)) {
    # Case 3: Treatment is not a factor, a conditional, a numeric, a logical or
    # a character. It is some other type. Warn user before converting to
    # numeric.
    warning(paste("Treatment variables which are not factors, numeric or",
                  " character are converted into numeric.\nIt is STRONGLY",
                  "recommended to do this conversion yourself to ensure it",
                  "proceeds as you expect."))
    treatment <- as.numeric(treatment)
  }
  # Case 0: treatment is numeric, logical, or character. No additional steps
  # needed.

  m_collapse[, index == "t"] <- treatment

  differing <- duplicated(m_collapse[, index == "u"])
  if (any(differing)) {
    noncon <- m_collapse[differing, index == "u", drop = FALSE]

    # Format data using utils::str or base formatting
    if (nrow(noncon) >= 6) {
      data_text <- paste(capture.output(noncon[1:5, , drop = FALSE]), collapse = "\n")
      data_text <- paste0(data_text, "\n...")
    } else {
      data_text <- paste(capture.output(noncon), collapse = "\n")
    }

    stop(paste0("Units of assignment with non-constant treatment, block or forcing:\n",
                data_text, "\n",
                "Each of treatment assignment, block and forcing must be ",
                "constant within unit of assignment."))
  }

  return(new("StudySpecification",
             structure = m_collapse,
             column_index = index,
             type = type,
             unit_of_assignment_type = autype,
             call = call))
}

##' @title Generates a \code{StudySpecification} object with the given
##'   specifications.
##'
##' @description Generate a randomized control treatment StudySpecification
##'   ([rct_spec()]), or an observational StudySpecification ([obs_spec()]), or
##'   a regression discontinuity StudySpecification ([rd_spec()]).
##'
##' @details The formula should include exactly one [unit_of_assignment()] to
##'   identify the units of assignment (one or more variables). (\code{uoa},
##'   \code{cluster}, or \code{unitid} are synonyms for
##'   \code{unit_of_assignment}; the choice of which has no impact on the
##'   analysis. See below for a limited exception in which the
##'   \code{unit_of_assignment} specification may be omitted.) If defining an
##'   \code{rd_spec}, the formula must also include a [forcing()] entry. The
##'   formula may optionally include a [block()] as well. Each of these can take
##'   in multiple variables, e.g. to pass both a household ID and individual ID
##'   as unit of assignment, use \code{uoa(hhid, iid)} and not \code{uoa(hhid) +
##'   uoa(iid)}.
##'
##'   The treatment variable passed into the left-hand side of \code{formula}
##'   can either be \code{logical}, \code{numeric}, or \code{character}. If it
##'   is anything else, it attempts conversion to one of those types (for
##'   example, \code{factor} and \code{ordered} are converted to \code{numeric}
##'   if the levels are \code{numeric}, otherwise to \code{character}). If the
##'   treatment is not \code{logical} or \code{numeric} with only values 0 and
##'   1, in order to generate weights with [ate()] or [ett()], the
##'   \code{dichotomy} argument must be used in those functions to identify the
##'   treatment and control groups. See [ett()] for more details on specifying a
##'   \code{dichotomy}.
##'
##'   There are a few aliases for each version.
##'
##'   If the formula excludes a \code{unit_of_assignment()}, data merges are
##'   performed on row order. Such formulas can also be passed as the
##'   specification argument to lmitt(), and that is their primary intended use
##'   case. It is recommended that each formula argument passed to
##'   *_specification() include a \code{unit_of_assignment()}, \code{uoa()} or
##'   \code{cluster()} term identifying the key variable(s) with which
##'   \code{StudySpecification} data is to be merged with analysis data.
##'   Exceptions to this rule will be met with a warning. To disable the
##'   warning, run \code{options("propertee_warn_on_no_unit_of_assignment" =
##'   FALSE)}.
##'
##'   The units of assignment, blocks, and forcing variables must be
##'   \code{numeric} or \code{character}. If they are otherwise, an attempt is
##'   made to cast them into \code{character}.
##'
##' @param formula a \code{formula} defining the \code{StudySpecification}
##'   components. See `Details` for specification.
##' @param data the data set from which to build the StudySpecification. Note
##'   that this data need not be the same as used to estimate the treatment
##'   effect; rather the \code{data} passed should contain information about the
##'   units of treatment assignment (as opposed to the units of analysis).
##' @param subset optional, subset the data before creating the
##'   \code{StudySpecification} object
##' @param na.fail If \code{TRUE} (default), any missing data found in the
##'   variables specified in \code{formula} (excluding treatment) will trigger
##'   an error. If \code{FALSE}, non-complete cases will be dropped before the
##'   creation of the \code{StudySpecification}
##' @return a \code{StudySpecification} object of the requested type for use in
##'   further analysis.
##' @export
##' @rdname StudySpecification_objects
##' @examples
##' data(simdata)
##' spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2) + block(bid),
##'                   data = simdata)
##'
##' data(schooldata)
##' spec <- obs_spec(treatment ~ unit_of_assignment(schoolid) + block(state),
##'                   data = schooldata)
rct_spec <- function(formula,
                       data,
                       subset = NULL,
                       na.fail = TRUE) {
  .check_spec_formula(formula)
  if (!is.null(substitute(subset))) {
    # Moving this inside .new_StudySpecification causes some scoping issues; it
    # should probably be inside there, but this works fine and is easier to
    # handle #209
    subset <- eval(substitute(subset), data, parent.frame())
  }
  return(.new_StudySpecification(form = formula,
                     data = data,
                     type = "RCT",
                     subset = subset,
                     call = match.call(),
                     na.fail = na.fail))
}

##' @export
##' @rdname StudySpecification_objects
rd_spec <- function(formula,
                      data,
                      subset = NULL,
                      na.fail = TRUE) {
  .check_spec_formula(formula, allow_forcing = TRUE)
  if (!is.null(substitute(subset))) {
    # Moving this inside .new_StudySpecification causes some scoping issues; it
    # should probably be inside there, but this works fine and is easier to
    # handle #209
    subset <- eval(substitute(subset), data, parent.frame())
  }
  return(.new_StudySpecification(form = formula,
                     data = data,
                     type = "RD",
                     subset = subset,
                     call = match.call(),
                     na.fail = na.fail))
}

##' @export
##' @rdname StudySpecification_objects
obs_spec <- function(formula,
                       data,
                       subset = NULL,
                       na.fail = TRUE) {
  .check_spec_formula(formula)
  if (!is.null(substitute(subset))) {
    # Moving this inside .new_StudySpecification causes some scoping issues; it
    # should probably be inside there, but this works fine and is easier to
    # handle #209
    subset <- eval(substitute(subset), data, parent.frame())
  }
  return(.new_StudySpecification(form = formula,
                     data = data,
                     type = "Obs",
                     subset = subset,
                     call = match.call(),
                     na.fail = na.fail))
}


################### Aliases
##' @export
##' @rdname StudySpecification_objects
rct_specification <- rct_spec

##' @export
##' @rdname StudySpecification_objects
rd_specification <- rd_spec

##' @export
##' @rdname StudySpecification_objects
obs_specification <- obs_spec

##' @export
##' @rdname StudySpecification_objects
obsstudy_spec <- obs_spec

##' @export
##' @rdname StudySpecification_objects
obsstudy_specification <- obs_spec


##' @title Show a \code{StudySpecification}
##' @description Display information about a \code{StudySpecification} object
##' @param object \code{StudySpecification} object, usually a result of a call
##'   to [rct_spec()], [obs_spec()], or [rd_spec()].
##' @return \code{object}, invisibly.
##' @export
setMethod("show", "StudySpecification", function(object) {
  spectype <- switch(object@type,
                     "RCT" = "Randomized Control Trial",
                     "RD" = "Regression Discontinuity StudySpecification",
                     "Obs" = "Observational Study")
  uoatype <- switch(object@unit_of_assignment_type,
                    "unit_of_assignment" = "Unit of Assignment",
                    "cluster" = "Cluster",
                    "unitid"  = "Unitid")

  cat(spectype)
  cat("\n\n")


  vartab <- var_table(object)
  # Add a nice separating line for printing
  vartab <- rbind(c("---------", "---------"),
                  vartab)
  print(data.frame(vartab), row.names = FALSE, right = FALSE)

  cat("\n")
  invisible(object)
})

##' After calling \code{model.frame()} on the formula input to
##' \code{.new_StudySpecification()}, the names of the columns will include
##' function names, e.g. "block(blockvar)". This function strips all these.
##'
##' @title (Internal) Rename columns to strip function calls
##' @param modframe A \code{data.frame}.
##' @return The \code{data.frame} with function calls removed
##' @keywords internal
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

##' @title Extract Variable Names from \code{StudySpecification}
##'
##' @description Methods to extract the variable names to the elements of the
##'   structure of the \code{StudySpecification} (e.g. treatment, unit of
##'   analysis, etc)
##'
##' @details When \code{compress} is \code{TRUE}, the result will always have
##'   two columns. When \code{FALSE}, the result will have number of columns
##'   equal to the largest number of variables in a particular role, plus one.
##'   E.g., a call such as \code{rct_spec(z ~ unitid(a, b, c, d) ...} will have
##'   4+1=5 columns in the output matrix with \code{compress = FALSE}.
##'
##'   When \code{report_all} is \code{TRUE}, the matrix is guaranteed to have 3
##'   rows (when the \code{specification} is an RCT or Obs) or 4 rows (when the
##'   \code{specification} is a RD), with empty variable entries as appropriate.
##'   When \code{FALSE}, the matrix will have minimum 2 rows (treatment and unit
##'   of assignment/unitid/cluster), with additional rows for blocks and forcing
##'   if included in the \code{StudySpecification}.
##' @param specification a \code{StudySpecification} object
##' @param compress should multiple variables be compressed into a
##'   comma-separated string? Default \code{TRUE}. If \code{FALSE}, multiple
##'   columns can be created instead.
##' @param report_all should we report all possible structures even if they
##'   don't exist in the \code{StudySpecification}? Default \code{FALSE}.
##' @return \code{var_table} returns the requested table. \code{var_names}
##'   returns a vector of variable names.
##' @export
##' @rdname StudySpecification_var_names
##' @order 1
##' @examples
##' spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
##' var_table(spec)
##' var_table(spec, compress = FALSE)
##' var_names(spec, "t")
##' var_names(spec, "u")
##' var_names(spec, "b")
var_table <- function(specification, compress = TRUE, report_all = FALSE) {
  uoatype <- switch(specification@unit_of_assignment_type,
                    "unit_of_assignment" = "Unit of Assignment",
                    "cluster" = "Cluster",
                    "unitid"  = "Unitid")

  # Start with the "table" as a list; easier to handle non-equal number of
  # elements.
  rows <- list()
  rows[["t"]] <- c("Treatment", var_names(specification, "t"))
  rows[["u"]] <- c(uoatype    , var_names(specification, "u"))
  rows[["b"]] <- c("Block"    , var_names(specification, "b"))
  if (specification@type == "RD") {
    rows[["f"]] <- c("Forcing"  , var_names(specification, "f"))
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

##' @param type one of "t", "u", "b", "f"; for "treatment",
##'   "unit_of_assignment", "block", and "forcing" respectively
##' @param implicitBlocks If the \code{StudySpecification} is created without
##'   blocks, setting this to \code{TRUE} will return "\code{.blocks_internal}"
##'   as the variable name corresponding to the blocks.
##' @export
##' @rdname StudySpecification_var_names
##' @order 2
var_names <- function(specification, type, implicitBlocks = FALSE) {
  stopifnot(inherits(specification, "StudySpecification"))
  stopifnot(length(type) == 1)
  stopifnot(type %in% c("t", "u", "b", "f"))
  if (type == "b" && !("b" %in% specification@column_index) && implicitBlocks) {
    return(".blocks_internal")
  }
  return(names(specification@structure)[specification@column_index == type])
}

##' @title Check for variable agreement within units of assignment
##'
##' @description Useful for debugging purposes to ensure that there is
##'   concordance between variables in the \code{StudySpecification} and data.
##'
##' @details Consider the following scenario: A \code{StudySpecification} is
##'   generated from some dataset, "data1", which includes a block variable
##'   "b1". Within each unique unit of assignment/unitid/cluster of "data1", it
##'   must be the case that "b1" is constant. (Otherwise the creation of the
##'   \code{StudySpecification} will fail.)
##'
##'   Next, a model is fit which includes weights generated from the
##'   \code{StudySpecification}, but on dataset "data2". In "data2", the block
##'   variable "b1" also exists, but due to some issue with data cleaning, does
##'   not agree with "b1" in "data1".
##'
##'   This could cause errors, either directly (via actual error messages) or
##'   simply produce nonsense results. [specification_data_concordance()] is
##'   specificationed to help debug these scenarios by providing information on
##'   whether variables in both the data used in the creation of
##'   \code{specification} ("data1" in the above example) and some new dataset,
##'   \code{data}, ("data2" in the above example) have any inconsistencies.
##'
##' @param specification a \code{StudySpecification} object
##' @param data a new data set, presumably not the same used to create
##'   \code{specification}.
##' @param by optional; named vector or list connecting names of variables in
##'   \code{specification} to variables in \code{data}. Names represent
##'   variables in \code{specification}; values represent variables in
##'   \code{data}. Only needed if variable names differ.
##' @param warn_on_nonexistence default \code{TRUE}. If a variable does not
##'   exist in \code{data}, should this be flagged? If \code{FALSE}, silently
##'   move on if a variable doesn't exist in \code{data}.
##' @return invisibly \code{TRUE} if no warnings are produced, \code{FALSE} if
##'   any warnings are produced.
##' @export
specification_data_concordance <- function(specification,
                                    data,
                                    by = NULL,
                                    warn_on_nonexistence = TRUE) {
  if (!is.null(by)) {
    specification <- .update_by(specification, data, by)
  }

  if (!all(var_names(specification, "u" ) %in% names(data))) {
    stop(paste("Missing unit of assignment/unitid/cluster variables",
               "in data set.\nAre you missing a `by=` argument?"))
  }

  merged <- merge(specification@structure, data,
                  all = TRUE,
                  by = var_names(specification, "u"))

  # Returns TRUE if any warnings were printed, FALSE otherwise
  .check <- function(type, specification, data) {
    anywarnings <- FALSE

    vnames <- var_names(specification, type)
    # If this particular type of variable exists in the specification...
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

  treat_had_warnings <- .check("t", specification, merged)
  block_had_warnings <- .check("b", specification, merged)
  force_had_warnings <- .check("f", specification, merged)

  if (any(treat_had_warnings,
          block_had_warnings,
          force_had_warnings)) {
    return(invisible(FALSE))
  } else {
    return(invisible(TRUE))
  }
}
