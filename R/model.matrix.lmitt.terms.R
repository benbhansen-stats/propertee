setGeneric("model.matrix")

#' Build a model matrix from a \code{lmitt.terms} object and its associated
#' model frame
#' 
#' The model matrix reflects the formula in the \code{model_form} attribute of
#' the \code{lmitt.terms} object passed to the \code{object} argument, hence the
#' need for the \code{data} argument to have columns from the model frame rather
#' than the dataframe from which the model frame may have been generated.
#' Aside from the intercept column, columns in the output matrix are
#' block-centered if the \code{absorb} attribute of \code{object} is TRUE; the
#' block means are stored as a \code{block.means} attribute of the model matrix
#' in this case.
#' @returns matrix, possibly with a \code{block.means} attribute storing the
#' column means by block
#' @param object \code{lmitt.terms}.
#' @param data model frame generated using \code{object}.
#' @param ... Additional arguments for forming the model matrix such as
#' \code{xlev}. See \code{stats::model.matrix.default()}.
#' @exportS3Method stats::model.matrix
#' @examples
#' data(simdata)
#' spec <- rct_spec(z ~ 1, simdata)
#' form <- build_lmitt_formula("y", "1", rhstype = "intercept", absorb = FALSE,
#'                             specification = quote(spec),
#'                             data = quote(simdata))
#' tt <- stats::terms(form, spec)
#' mf <- stats::model.frame(tt, simdata)
#' model.matrix(tt, mf)
#' 
#' blkd_spec <- rct_spec(z ~ unitid(uoa1, uoa2) + block(bid), simdata)
#' form <- build_lmitt_formula("y", "1", rhstype = "intercept", absorb = TRUE,
#'                             specification = quote(blkd_spec),
#'                             data = quote(simdata))
#' tt <- stats::terms(form, spec)
#' mf <- stats::model.frame(tt, simdata)
#' stats::model.matrix(tt, mf)
model.matrix.lmitt.terms <- function(object, data, ...) {
  mc <- match.call()
  mc[[1L]] <- quote(stats::model.matrix.default)
  mc[[2L]] <- attr(object, "model_form")
  mm <- eval(mc, parent.frame(1L))
  colnames(mm) <- gsub("\\:", "_", colnames(mm))
  if (attr(object, "absorb")) {
    block.centered <- areg.center(
      mm[,seq(2,ncol(mm)),drop=FALSE],
      as.factor(eval(mc$data, parent.frame(1L))[["..block.."]]),
      mget("(weights)",
           envir = list2env(eval(mc$data, parent.frame(1L))),
           ifnotfound = list(NULL))$`(weights)`
    )
    mm[,seq(2,ncol(mm))] <- block.centered
    attr(mm, "block.means") <- attr(block.centered, "block.means")
  }
  return(mm)
}

#' Return a model matrix associated with a model fit by \code{lmitt()}
#' @param object \code{teeMod}, or a model fit by \code{lmitt()}.
#' @param ... Additional arguments for forming the model matrix such as
#' \code{xlev}. See \code{stats::model.matrix.default()}.
#' @exportS3Method stats::model.matrix
#' @examples
#' data(simdata)
#' spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)
#' tm <- lmitt(y ~ 1, spec, simdata)
#' stats::model.matrix(tm)
#' 
#' blkd_spec <- rct_spec(z ~ unitid(uoa1, uoa2) + block(bid), simdata)
#' tm <- lmitt(y ~ 1, blkd_spec, simdata, absorb = TRUE)
#' stats::model.matrix(tm)
model.matrix.teeMod <- function(object, ...) {
  mf <- model.frame(object, xlev = object$xlevels, ...)
  model.matrix(terms(object), mf, ...)
}