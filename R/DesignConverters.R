##' @title Convert \code{Design} between types
##'
##' @description Convert a \code{Design} between a observational study, a
##'   randomized control trial, and a regression discontinuity (created from
##'   [obs_design()], [rct_design()] and [rd_design()] respectively).
##'
##' @param Design a \code{Design} to convert
##' @param data converting to an RD requires adding a \code{forcing} variable,
##'   which requires access to the original data.
##' @param ... Ignored.
##' @param loseforcing converting from RD to another \code{Design} type will
##'   error to avoid losing the forcing variable. Setting \code{loseforcing =
##'   TRUE} allows the conversion to automatically drop the forcing variable.
##'   Default \code{FALSE}.
##' @param forcing converting to an RD requires adding a \code{forcing}
##'   variable. This should be entered as a formula which would be passed to
##'   [update()], e.g. \code{forcing = . ~ . + forcing(forcevar)}.
##' @return \code{Design} of the updated type
##' @export
##' @importFrom stats as.formula update
##' @rdname designconversion
##' @examples
##' des <- rct_design(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
##' des
##' as_obs_design(des)
##' as_rd_design(des, simdata, forcing = ~ . + forcing(force))
##' des2 <- rd_design(o ~ uoa(uoa1, uoa2) + forcing(force), data = simdata)
##' des2
##' # as_rct_design(des2) # this will produce an error
##' as_rct_design(des2, loseforcing = TRUE)
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
  if (!inherits(forcing, "formula")) {
    stop('`forcing` must be entered as a formula such as "~ . + forcing(a, b)"')
  }

  origcall <- Design@call

  Design <- rd_design(formula =
                        stats::update(stats::as.formula(Design@call[[2]]),
                                      forcing),
                      data = data,
                      subset = eval(Design@call[["subset"]]))

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
