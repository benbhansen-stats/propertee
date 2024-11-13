##' @title Convert \code{StudySpecification} between types
##'
##' @description Convert a \code{StudySpecification} between a observational
##'   study, a randomized control trial, and a regression discontinuity (created
##'   from \code{obs_spec}, \code{rct_spec} and \code{rd_spec} respectively).
##'
##' @param StudySpecification a \code{StudySpecification} to convert
##' @param data converting to an RD requires adding a \code{forcing} variable,
##'   which requires access to the original data.
##' @param ... Ignored.
##' @param loseforcing converting from RD to another \code{StudySpecification}
##'   type will error to avoid losing the forcing variable. Setting
##'   \code{loseforcing = TRUE} allows the conversion to automatically drop the
##'   forcing variable. Default \code{FALSE}.
##' @param forcing converting to an RD requires adding a \code{forcing}
##'   variable. This should be entered as a formula which would be passed to
##'   \code{update}, e.g. \code{forcing = . ~ . + forcing(forcevar)}.
##' @return \code{StudySpecification} of the updated type
##' @export
##' @importFrom stats as.formula update
##' @rdname specificationconversion
##' @examples
##' spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
##' spec
##' as_obs_spec(spec)
##' as_rd_spec(spec, simdata, forcing = ~ . + forcing(force))
##' spec2 <- rd_spec(o ~ uoa(uoa1, uoa2) + forcing(force), data = simdata)
##' spec2
##' # as_rct_spec(spec2) # this will produce an error
##' as_rct_spec(spec2, loseforcing = TRUE)
as_rct_spec <- function(StudySpecification, ..., loseforcing = FALSE) {
  if (StudySpecification@type == "RD") {
    if (!loseforcing) {
      stop(paste("Any `forcing` variables will be dropped. Pass option",
                 "`loseforcing = TRUE` to proceed."))
    }
    StudySpecification <- .remove_forcing(StudySpecification)
  }
  StudySpecification@type <- "RCT"
  StudySpecification@call[[1]] <- as.name("rct_spec")
  validObject(StudySpecification)
  return(StudySpecification)
}

##' @export
##' @rdname specificationconversion
as_obs_spec <- function(StudySpecification, ..., loseforcing = FALSE) {
  if (StudySpecification@type == "RD") {
    if (!loseforcing) {
      stop(paste("Any `forcing` variables will be dropped. Pass option",
                 " `force = TRUE` to proceed."))
    }
    StudySpecification <- .remove_forcing(StudySpecification)
  }
  StudySpecification@type <- "Obs"
  StudySpecification@call[[1]] <- as.name("obs_spec")
  validObject(StudySpecification)
  return(StudySpecification)
}

##' @export
##' @rdname specificationconversion
as_rd_spec <- function(StudySpecification, data, ..., forcing) {
  if (!inherits(forcing, "formula")) {
    stop('`forcing` must be entered as a formula such as "~ . + forcing(a, b)"')
  }

  origcall <- StudySpecification@call

  StudySpecification <- rd_spec(formula =
                        stats::update(stats::as.formula(StudySpecification@call[[2]]),
                                      forcing),
                      data = data,
                      subset = eval(StudySpecification@call[["subset"]]))

  StudySpecification@call <- origcall
  StudySpecification@call[[1]] <- as.name("rd_spec")
  StudySpecification@call[[2]] <- as.call(stats::update(stats::as.formula(origcall[[2]]),
                                            forcing))

  validObject(StudySpecification)
  return(StudySpecification)
}

##' In preparation for converting an RD \code{StudySpecification} to another
##' \code{StudySpecification}, this will strip the forcing variable entirely. It
##' is removed from the data (both \code{@structure} and \code{@column_index}),
##' as well as from the formula stored in \code{@call}.
##'
##' Note that the output \code{StudySpecification} will fail a validity check
##' (with \code{validObject()}) due to an RD \code{StudySpecification} requiring
##' a forcing variable, so change the \code{@type} immediately.
##'
##' @title (Internal) Removes the forcing column entirely from a
##'   \code{StudySpecification}
##' @param spec A \code{StudySpecification}
##' @return The \code{StudySpecification} without any forcing variable
##' @keywords internal
.remove_forcing <- function(spec) {
  spec@structure <- spec@structure[, spec@column_index != "f"]
  spec@column_index <- spec@column_index[spec@column_index != "f"]

  # Remove the "forcing" element from the call formula
  spec@call[[2]][[3]] <- spec@call[[2]][[3]][!grepl("forcing",
                                                  spec@call[[2]][[3]])]
  return(spec)
}
