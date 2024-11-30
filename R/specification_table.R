##' @title Table of elements from a \code{StudySpecification}
##'
##' @description Produces a table (1-dimensional, or 2-dimensional if \code{y}
##'   is specified) of the elements of the \code{StudySpecification}.
##'
##' @param specification A \code{StudySpecification} object
##' @param x One of "treatment", "unit of assignment", (synonym "uoa"), "block".
##'   Abbreviations are accepted. "unit of assignment" can be replaced by
##'   "unitid" or "cluster" if the \code{StudySpecification} was created with
##'   that element.
##' @param y Optionally, another string similar to \code{x}. A 1-dimensional
##'   table is produced if \code{y} is left at its default, \code{NULL}.
##' @param sort Ignored if \code{y} is not \code{NULL}. If \code{FALSE}
##'   (default), one-way table is sorted according to "names" of levels. If set
##'   to \code{TRUE}, one-way table is sorted according to values.
##' @param decreasing If \code{sort} is \code{TRUE}, choose whether to sort
##'   descending (\code{TRUE}, default) or ascending (\code{FALSE}).
##' @param use_var_names If \code{TRUE}, name dimensions of table returned by
##'   variable names. If \code{FALSE} (default), name by their function (e.g.
##'   "treatment" or "blocks"). Passing the \code{dnn} argument in \code{...}
##'   (an argument of [table()]) overrides whatever is requested here.
##' @param ... additional arguments [table()]
##' @return A table of the requested variables.
##' @export
##' @rdname specification_table
##' @order 2
##' @examples
##' data(simdata)
##' spec <- obs_spec(z ~ unit_of_assignment(uoa1, uoa2) + block(bid),
##'                   data = simdata)
##' specification_table(spec, "treatment")
##' specification_table(spec, "treatment", "block", sort = TRUE, use_var_names = TRUE)
stable <- function(specification,
                   x,
                   y = NULL,
                   sort = FALSE,
                   decreasing = TRUE,
                   use_var_names = FALSE,
                   ...) {

  # Internal function to match partial names and standardize across uoa/unit of
  # assignment
  .fix_xy <- function(z, which) {
    if (!is.null(z)) {
      z <- tolower(z)
      if (grepl(paste0("^", z), "treatments")) {
        z <- "treatment"
      } else if (grepl(paste0("^", z), "blocks")) {
        z <- "blocks"
      } else if (grepl(paste0("^", z), "uoa") ||
                 grepl(paste0("^", z), "unit of assignment") ||
                 grepl(paste0("^", z), "units of assignment")) {
        z <- "units_of_assignment"
      } else if (grepl(paste0("^", z), "clusters")) {
        z <- "clusters"
      } else {
        stop(paste("Invalid input", which))
      }
    }
    return(z)
  }

  x <- .fix_xy(x, "x")
  y <- .fix_xy(y, "y")

  # Internal function to obtain proper data from the StudySpecification structure, and if
  # its multi-dimensional (e.g. 2 variables uniquely identify block) collapse
  # appropriately.
  .get_xydat <- function(z, specification) {
    if (!is.null(z)) {

      # Since .fix_xy ensures that `x` and `y` contain valid function names,
      # this executes them.
      zdat <- match.fun(z)(specification)

      if (ncol(zdat) == 0) {
        # This is hit if user request an element not in the specification (e.g.
        # requests block, but specification is created without blocks
        warning(paste(z, "not found in specification"))
        return(NULL)
      }

      # If we have multiple variables in the specification definition, collapse them
      if (ncol(zdat) > 1) {
        zdat <- apply(zdat, 1, paste, collapse = ",")
      } else {
        # If we have a single variable, drop second dimension
        zdat <- zdat[, 1]
      }
      return(zdat)
    } else {
      return(NULL)
    }
  }

  xdat <- .get_xydat(x, specification)
  ydat <- .get_xydat(y, specification)

  # If user passes in `dnn`, we need to not pass our own when creating table
  # below.
  dots_have_dnn <- "dnn" %in% names(list(...))

  # If x or y were asking for units of assignment, `x` or `y` are now
  # "units_of_assignment". Replace "_" with " ". Does nothing in other cases.
  x <- gsub("_", " ", x)
  y <- gsub("_", " ", y)


  if (is.null(ydat) | is.null(xdat)) {
    if (is.null(xdat)) {
      # We're here if the user requsted both an x and a y, but the x is invalid
      # (e.g. specification was created without block, but user requested blocks). In
      # this case, it should be equivalent to a user only requesting x, so move
      # y into x. (Note that warning to use arises earlier in `.get_xydat`.)
      x <- y
      xdat <- ydat
    }
    # 1-dimensional table

    if (dots_have_dnn) {
      outtable <- table(xdat, ...)
    } else {
      if (use_var_names) {
        x <- paste(var_names(specification, substr(x, 1, 1)), collapse = ", ")
      }
      outtable <- table(xdat, dnn = list(x), ...)
    }

    if (sort) {
      outtable <- sort(outtable, decreasing = decreasing)
    }
  } else {
    # 2-dimensional table

    if (dots_have_dnn) {
      outtable <- table(ydat, xdat, ...)
    } else {
      if (use_var_names) {
        x <- paste(var_names(specification, substr(x, 1, 1)))
        y <- paste(var_names(specification, substr(y, 1, 1)))
      }
      outtable <- table(ydat, xdat, dnn = list(y, x), ...)
    }

  }

  return(outtable)


}

##' @export
##' @rdname specification_table
##' @order 1
specification_table <- stable
