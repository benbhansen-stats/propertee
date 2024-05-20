##' @title Table of elements from a \code{Design}
##'
##' @description Produces a table (1-dimensional, or 2-dimensional if \code{y}
##'   is specified) of the elements of the \code{Design}.
##'
##' @param design A \code{Design} object
##' @param x One of "treatment", "unit of assignment", (synonym "uoa"), "block".
##'   Abbreviations are accepted. "unit of assignment" can be replaced by
##'   "unitid" or "cluster" if the \code{Design} was created with that element.
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
##' @param treatment_binary Should the treatment (if requested) be dichotomized
##'   of \code{design} contains a \code{dichotomy}? Ignored if \code{design}
##'   does not contain a \code{dichotomy}, or if neither \code{x} or \code{y} is
##'   "treatment".
##' @param ... additional arguments [table()]
##' @return A table of the requested variables.
##' @export
##' @rdname design_table
##' @order 2
##' @examples
##' data(simdata)
##' des <- obs_design(z ~ unit_of_assignment(uoa1, uoa2) + block(bid),
##'                   data = simdata)
##' design_table(des, "treatment")
##' design_table(des, "treatment", "block", sort = TRUE, use_var_names = TRUE)
dtable <- function(design,
                   x,
                   y = NULL,
                   sort = FALSE,
                   decreasing = TRUE,
                   use_var_names = FALSE,
                   treatment_binary = TRUE,
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

  # Internal function to obtain proper data from the Design structure, and if
  # its multi-dimensional (e.g. 2 variables uniquely identify block) collapse
  # appropriately.
  .get_xydat <- function(z, design) {
    if (!is.null(z)) {

      # Since .fix_xy ensures that `x` and `y` contain valid function names,
      # this executes them.
      if (z == "treatment" & is_binary_or_dichotomized(design)) {
        # Pass down `treatment_binary`
        zdat <- match.fun(z)(design, binary = treatment_binary)
      } else {
        zdat <- match.fun(z)(design)
      }
      if (ncol(zdat) == 0) {
        # This is hit if user request an element not in the design (e.g.
        # requests block, but design is created without blocks
        warning(paste(z, "not found in design"))
        return(NULL)
      }

      # If we have multiple variables in the design definition, collapse them
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

  xdat <- .get_xydat(x, design)
  ydat <- .get_xydat(y, design)

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
      # (e.g. design was created without block, but user requested blocks). In
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
        x <- paste(var_names(design, substr(x, 1, 1)), collapse = ", ")
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
        x <- paste(var_names(design, substr(x, 1, 1)))
        y <- paste(var_names(design, substr(y, 1, 1)))
      }
      outtable <- table(ydat, xdat, dnn = list(y, x), ...)
    }

  }

  return(outtable)


}

##' @export
##' @rdname design_table
##' @order 1
design_table <- dtable
