#' @include lmitt.formula.R
NULL

setGeneric("model.frame")

#' Build a model frame from a \code{lmitt.terms} object
#' 
#' The model frame will reflect the formula stored in the \code{model_form}
#' attribute of the \code{lmitt.terms} object passed to the \code{formula}
#' argument. If the \code{absorb} attribute of that object is TRUE, the model
#' frame will have an additional column for block identifiers. The \code{terms}
#' attribute of the model frame is the original \code{lmitt.terms} object passed
#' to the \code{formula} argument.
#' @param formula \code{lmitt.terms}.
#' @param data dataframe.
#' @param ... Additional arguments to inform the model frame generation, such
#' as weights, subset, or offset calls.
#' @returns dataframe with a \code{terms} attribute storing the original
#' \code{lmitt.terms} object
#' @importFrom stats terms
#' @exportS3Method stats::model.frame
#' @examples
#' data(simdata)
#' spec <- rct_spec(z ~ 1, simdata)
#' form <- build_lmitt_formula("y", "1", rhstype = "intercept", absorb = FALSE,
#'                             specification = quote(spec),
#'                             data = quote(simdata))
#' tt <- terms(form, spec)
#' stats::model.frame(tt, simdata)
#' 
#' blkd_spec <- rct_spec(z ~ unitid(uoa1, uoa2) + block(bid), simdata)
#' form <- build_lmitt_formula("y", "1", rhstype = "intercept", absorb = TRUE,
#'                             specification = quote(blkd_spec),
#'                             data = quote(simdata))
#' tt <- terms(form, spec)
#' stats::model.frame(tt, simdata)
model.frame.lmitt.terms <- function(formula, data, ...) {
  mf <- match.call()
  mf[[1L]] <- quote(stats::model.frame)

  # create the assignment variable column with the name expected by the
  # lmitt.terms object
  rhs <- if (length(formula) == 2) formula[[2L]] else formula[[3L]]
  a.call <- match.call(get(rhs[[2L]][[1L]]), rhs[[2L]])
  a.call$data <- data
  rhs[[2L]] <- a.call
  data[[attr(formula, "term.labels")[1L]]] <- eval(rhs[[2L]],
                                                   environment(formula))
  
  # use model.frame.default to build the model frame
  mf$data <- data
  mf[[2L]] <- attr(formula, "model_form")
  mm <- mf
  mf <- eval(mf, parent.frame())
  
  # add a ..block.. column if model.matrix will eventually need it for block
  # centering
  if (attr(formula, "absorb")) {
    specification <- tryCatch(
      .get_spec(),
      error = function(e) {
        # if .get_spec fails, just look in the parent frame for `specification`
        if (is.null(spec_call <- a.call$specification)) {
          spec_call <- quote(specification)
        }
        eval(spec_call, environment(formula))
      }
    )
    blocks <- blocks(specification,
                     data,
                     all.x = TRUE,
                     implicit = TRUE)[,1]
    blocks <- blocks[setdiff(seq_along(blocks), stats::na.action(mf))]
    mf[["..block.."]] <- blocks
  }
  
  attr(mf, "terms") <- formula
  return(mf)
}
