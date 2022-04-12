##' Converts between Design types
##' @title Design conversion
##' @param Design Design to convert
##' @param data Converting to an RD requires adding a `forcing` variable, which
##'   requires access to the original data.
##' @param ... No addiitonal options at present
##' @param loseforcing Converting from RD to another design will error to avoid
##'   losing the forcing variable. Setting `loseforcing = TRUE` allows the
##'   conversion to automatically drop the forcing variable.
##' @param forcing Converting to an RD requires adding a `forcing` variable.
##'   This should be entered as the update to a formula, e.g.
##'   "~ . + forcing(forcevar)".
##' @return The design of the updated type
##' @export
##' @importFrom stats as.formula update
##' @rdname designconversion
as_rct_design <- function(Design, ..., loseforcing = FALSE) {
  if (Design@type == "RD") {
    if (!loseforcing) {
      stop("Any `forcing` variables will be dropped. Pass option `loseforcing = TRUE` to proceed.")
    }
    Design <- .remove_forcing(Design)
  }
  Design@type <- "RCT"
  Design@call[[1]] <- as.name("rct_design")
  validObject(Design)
  return(Design)
}

##' @export
##' @rdname designconversion
as_obs_design <- function(Design, ..., loseforcing = FALSE) {
  if (Design@type == "RD") {
    if (!loseforcing) {
      stop("Any `forcing` variables will be dropped. Pass option `force = TRUE` to proceed.")
    }
    Design <- .remove_forcing(Design)
  }
  Design@type <- "Obs"
  Design@call[[1]] <- as.name("obs_design")
  validObject(Design)
  return(Design)
}

##' @export
##' @rdname designconversion
as_rd_design <- function(Design, data, ..., forcing) {
  if (!is(forcing, "formula")) {
    stop('`forcing` must be entered as a formula such as "~ . + forcing(a, b)"')
  }

  origcall <- Design@call

  Design <- rd_design(formula = stats::update(stats::as.formula(Design@call[[2]]), forcing),
                      data = data,
                      subset = eval(Design@call[["subset"]]),
                      dichotomize = Design@dichotomization)

  Design@call <- origcall
  Design@call[[1]] <- as.name("rd_design")
  Design@call[[2]] <- as.call(stats::update(stats::as.formula(origcall[[2]]), forcing))

  validObject(Design)
  return(Design)
}

# Internal function
# Removes the forcing column to prepare conversion from RD to another type.
.remove_forcing <- function(d) {
  d@structure <- d@structure[, d@column_index != "f"]
  d@column_index <- d@column_index[d@column_index != "f"]

  # Remove the "forcing" element from the formula
  d@call[[2]][[3]] <- d@call[[2]][[3]][!grepl("forcing", d@call[[2]][[3]])]
  return(d)
}
