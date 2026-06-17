#' @include SandwichLayer.R
NULL

##' @title Preparing Peters-Belson adjustment of \code{lmitt()} estimates
##' @description \code{cov_adj()} takes a fitted model object and predicts
##'   outcomes under the counterfactual of being in the control condition. It
##'   stores in the prediction vector additional information such that standard
##'   errors of adjusted \code{lmitt()} estimates reflect contributions to
##'   sampling variability from the adjustment model coefficient estimates.
##' @details If \code{specification} is provided or can be found in the call
##'   stack, then prior to generating adjustments, \code{cov_adj()} will
##'   identify the treatment variable and replace all values with a reference
##'   level. If the treatment has logical type, this reference level is
##'   \code{FALSE}, and if it has numeric type, this is the smallest
##'   non-negative value (which means 0 for 0/1 binary). Factor treatments are
##'   not currently supported for \code{StudySpecification} objects.\cr\cr More
##'   broadly, the \code{set_to_reference} argument can be used to set columns
##'   to values specified by the user prior to generating predictions. The
##'   argument takes a named list where the names specify the columns whose
##'   values should be replaced, and the entries specify the replacement values.
##'   If an entry is given as \code{default()}, internal routines will determine
##'   the reference level (the first level of a factor/factorized character
##'   variable or the minimum of a numeric/logical variable). \cr\cr If
##'   \code{newdata} is provided, the values of the output vector represent
##'   adjustments for the outcomes in \code{newdata}. If \code{cov_adj()} is
##'   called in the \code{offset} argument of \code{lmitt()}, predictions will
##'   correspond to rows in the \code{data} argument of \code{lmitt()}.
##'   Otherwise, adjustments are made based on the data to which \code{model} is
##'   fit. \cr\cr The \code{by} argument specifies columns for merging datasets.
##'   The names of the argument specify columns in the dataframe for building
##'   \code{specification}, while the entries specify columns in the dataframe
##'   for fitting \code{model}.
##' @param model any model that inherits from a \code{glm}, \code{lm}, or
##'   \code{robustbase::lmrob} object.
##' @param newdata a dataframe; optional. See Details for how \code{newdata}
##' is found in the call stack if none is provided.
##' @param specification a \code{StudySpecification} object; optional. If not
##' provided, it can be retrieved from the call stack if passed to suitable
##' calls.
##' @param by named vector; optional.
##' @param set_to_reference named list; optional.
##' @return A \code{SandwichLayer} if \code{specification} is not NULL or a
##'   \code{StudySpecification} object is found in the call stack, otherwise a
##'   \code{PreSandwichLayer} object
##' @export
##' @example inst/examples/cov_adj.R
cov_adj <- function(model,
                    newdata = NULL,
                    specification =  NULL,
                    by = NULL,
                    set_to_reference = NULL) {
  if (is.null(specification)) {
    specification <- .get_spec(NULL_on_error = TRUE)
  }

  if (is.null(newdata)) {
    form <- .update_ca_model_formula(model, by, specification)
    newdata <- tryCatch(
      .get_data_from_model("cov_adj", form),
      error = function(e) {
        warning(paste("Could not find direct adjustment data in the call",
                      "stack, or it did not contain the columns specified in",
                      "`by`. Using the covariance adjustment data to generate",
                      "the covariance adjustments"), call. = FALSE)
        tryCatch({
          data_call <- model$call$data
          if (is.null(data_call)) {
            stop("`model` must be fit using a `data` argument")
          }
          data <- eval(data_call, envir = environment(formula(model)))
          stats::model.frame(form, data, na.action = na.pass)
        }, error = function(e) {
          stop(paste(e$message, "in covariance adjustment data"), call. = FALSE)
        })
      })
  } else {
    newdata <- .as_data_frame(newdata)
  }


  if (!is.null(set_to_reference)) {
    if (is.null(names(set_to_reference))) stop(
      "`set_to_reference` must be a named list"
    )
    missing_cols <- setdiff(names(set_to_reference), names(newdata))
    if (length(missing_cols) > 0) {
      warning(
        warningCondition(paste("Columns",
                               paste(missing_cols, collapse = ", "),
                               "from `set_to_reference` not found in newdata"),
                         call = "cov_adj")
      )
    }
    for (c in seq_along(set_to_reference)) {
      col <- names(set_to_reference)[c]
      val <- set_to_reference[[c]]
      if (inherits(val, "set_to_reference_default")) {
        if (is.numeric(newdata[[col]])) {
          val <- min(newdata[[col]], na.rm = TRUE)
        } else if (is.logical(newdata[[col]])) {
          val <- as.logical(min(newdata[[col]], na.rm = TRUE))
        } else if (!inherits(newdata[[col]], "factor")) {
          val <- levels(factor(newdata[[col]]))[1L]
        } else {
          val <- levels(newdata[[col]])[1L]
        }
      }
      newdata[[col]] <- val
    }
  }
  if (!is.null(specification)) {
    trt_name <- var_names(specification,'t')
    if (trt_name %in% names(newdata))
      if (is.numeric(treatment(specification)[, 1])) {
        newdata[[trt_name]] <- min(abs(treatment(specification)[, 1]))
      } else if (is.logical(treatment(specification)[, 1])) {
        newdata[[trt_name]] <- FALSE
      } else if (is.factor(treatment(specification)[, 1])) {
        newdata[[trt_name]] <- levels(treatment(specification)[, 1])[1]
      } else {
        warning(
          paste("If treatment is an independent variable of the covariance",
                "model, predictions may (incorrectly) include treatment",
                "contributions. For now, setting these contributions to 0 is",
                "implemented only for logical or numeric treatments, but this",
                "treatment is a factor.")
        )
      }
    if (specification@unit_of_assignment_type == "none") {
      newdata$..uoa.. <- rownames(newdata)
    }
  } else if (is.null(set_to_reference)) {
    warning(
      paste("Without a specification, post-treatment variables in the",
            "covariance adjustment model cannot be ensured not to contribute",
            "to predictions that offset study outcomes. Use the",
            "set_to_reference argument or pass the StudySpecification to",
            "cov_adj() to avoid this warning.")
    )
  }

  psl <- .make_PreSandwichLayer(model, newdata)

  if (is.null(specification)) {
    return(psl)
  } else {
    sl <- as.SandwichLayer(psl, specification, by = by, Q_data = newdata)
    return(sl)
  }
}

##' Instruct \code{cov_adj()} to find a default reference value for columns in
##' a covariance adjustment model
##' @details
##' When set as the entry of a list passed to the \code{set_to_reference}
##' argument of \code{cov_adj()}, \code{cov_adj()} uses the following
##' replacement values depending on the column type:
##'  - factor/character: minimum level
##'  - numeric: minimum value
##'  - logical: FALSE
##' @return empty list that inherits from class \code{set_to_reference_default}
##' @export
default <- function() {
  structure(list(), class = "set_to_reference_default")
}

##' @title (Internal) Add columns for merging covariance adjustment and direct
##'   adjustment samples to model formula
##' @details This function is typically used prior to
##'   \code{.get_data_from_model()} and incorporates information provided in a
##'   \code{by} argument to ensure the necessary columns for merging the two
##'   samples are included in any \code{model.frame()} calls.
##' @inheritParams cov_adj
##' @return formula
##' @keywords internal
.update_ca_model_formula <- function(model, by = NULL, specification = NULL) {
  form <- deparse1(formula(model))
  if (!is.null(specification)) {
    form <- paste(
      c(form, paste(var_names(specification, "u"), collapse = " + ")),
      collapse = " + "
    )
  }
  if (!is.null(by) && is.null(names(by))) {
    form <- paste(c(form, paste(by, collapse = " + ")), collapse = " + ")
  } else if (!is.null(names(by))) {
    names(by)[names(by) == ""] <- by[names(by) == ""]
    form <- paste(c(form, paste(names(by), collapse = " + ")), collapse = " + ")
  }

  form <- as.formula(form)
  return(form)
}
