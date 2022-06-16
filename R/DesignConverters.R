##' Converts between \code{Design} types
##' @title Design conversion
##' @param Design \code{Design} to convert
##' @param data Converting to an RD requires adding a \code{forcing} variable,
##'   which requires access to the original data.
##' @param ... No addiitonal options at present
##' @param loseforcing Converting from RD to another \code{Design} type will
##'   error to avoid losing the forcing variable. Setting \code{loseforcing =
##'   TRUE} allows the conversion to automatically drop the forcing variable.
##'   Default \code{FALSE}.
##' @param forcing Converting to an RD requires adding a \code{forcing}
##'   variable. This should be entered as the update to a formula, e.g.
##'   \code{"~ . + forcing(forcevar)"}.
##' @return \code{Design} of the updated type
##' @export
##' @importFrom stats as.formula update
##' @rdname designconversion
as_rct_design <- function(Design, ..., loseforcing = FALSE) {
  if (Design@type == "RD") {
    if (!loseforcing) {
      stop(paste("Any `forcing` variables will be dropped. Pass option",
                 "`loseforcing = TRUE` to proceed."))
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
      stop(paste("Any `forcing` variables will be dropped. Pass option",
                 " `force = TRUE` to proceed."))
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

  Design <- rd_design(formula =
                        stats::update(stats::as.formula(Design@call[[2]]),
                                      forcing),
                      data = data,
                      subset = eval(Design@call[["subset"]]),
                      dichotomy = Design@dichotomy)

  Design@call <- origcall
  Design@call[[1]] <- as.name("rd_design")
  Design@call[[2]] <- as.call(stats::update(stats::as.formula(origcall[[2]]),
                                            forcing))

  validObject(Design)
  return(Design)
}

##' In preparation for converting an RD \code{Design} to another \code{Design},
##' this will strip the forcing variable entirely. It is removed from the data
##' (both \code{@structure} and \code{@column_index}), as well as from the
##' formula stored in \code{@call}.
##'
##' Note that the output \code{Design} will fail a validity check (with
##' \code{validObject()}) due to an RD \code{Design} requiring a forcing
##' variable, so change the \code{@type} immediately.
##'
##' @title (Internal) Removes the forcing column entirely from a \code{Design}
##' @param des A \code{Design}
##' @return The \code{Design} without any forcing variable
##' @keywords internal
.remove_forcing <- function(des) {
  des@structure <- des@structure[, des@column_index != "f"]
  des@column_index <- des@column_index[des@column_index != "f"]

  # Remove the "forcing" element from the call formula
  des@call[[2]][[3]] <- des@call[[2]][[3]][!grepl("forcing",
                                                  des@call[[2]][[3]])]
  return(des)
}
