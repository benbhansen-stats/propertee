setGeneric("terms")

#' Build a terms object from a \code{lmitt.formula}
#' 
#' This method parses the terms in a \code{lmitt.formula} and creates a
#' terms-like object that is used to generate the model frame for an
#' \code{lmitt()} fit.
#' 
#' @returns A \code{lmitt.terms} object.
#' @details
#' The formula associated with the output--what one obtains calling 
#' \code{formula()} on the object--is the \code{lmitt.formula} passed to
#' \code{x} (see \code{formula.lmitt.terms()}), but the attributes of the output
#' reflect the formula stored in the \code{model_form} attribute. This allows
#' \code{model.frame()} to build a model frame based on \code{model_form}, which
#' has the exact terms to be included in the model fit, while future
#' \code{terms.lmitt.terms()} calls have the information from the
#' \code{lmitt.formula} to build the model frame from new data. The
#' \code{model_form} and \code{absorb} attributes are specific to the
#' \code{lmitt.terms} class, where the latter informs future
#' \code{model.matrix.lmitt.terms()} calls to block-center its columns. Since
#' \code{model.frame()} and \code{model.matrix()} are S3 generics, 
#' \code{lmitt.terms} is not a defined S4 class, but rather manually set within
#' \code{terms.lmitt.formula()}.
#' @param x \code{lmitt.formula}.
#' @param specification \code{StudySpecification}.
#' @param ... Additional arguments to inform attributes stored in the
#' \code{lmitt.terms} object.
#'  
#' @importFrom stats reformulate
#' @exportS3Method stats::terms
#' @examples
#' data(simdata)
#' spec <- rct_spec(z ~ 1, simdata)
#' form <- build_lmitt_formula("y", "1", rhstype = "intercept", absorb = FALSE)
#' terms(form, spec)
terms.lmitt.formula <- function(x, specification, ...) {
  mc <- match.call()
  # extract the components of the formula
  lhs <- x[[2L]]
  rhs <- x[[3L]]
  mod_special <- rhs[[3L]]
  modvar <- deparse1(mod_special[[2L]])
  
  # get the default terms object
  mt <- mc
  mt[[2L]] <- as.formula(paste(c(lhs, deparse(x[[1L]]), rhs), collapse = ""))
  mt[[1L]] <- quote(stats::terms)
  mt <- eval(mt, parent.frame(1L))
  attr(mt, ".Environment") <- environment(x)
  
  # make the formula that will be called during model fitting and store it in
  # the terms object
  tx_col <- paste0(var_names(specification, "t"), ".")
  if (mod_special$rhstype == "intercept") {
    mf_form <- as.formula(paste0(lhs, "~", tx_col),
                          env = environment(x))
  } else if (mod_special$rhstype == "continuous") {
    mf_form <- as.formula(paste0(lhs, "~", tx_col, "+", tx_col, ":",
                                 modvar, "+",
                                 modvar),
                          env = environment(x))
  } else {
    mf_form <- as.formula(paste0(lhs, "~", tx_col, ":",
                                 modvar, "+",
                                 modvar),
                          env = environment(x))
  }
  attr(mt, "model_form") <- mf_form
  
  # replace attributes of the lmitt.terms object that use names from the rhs of
  # x with attributes that use the names in mf_form
  terms.deparsed <- c(deparse1(lhs), tx_col, if (modvar != 1) modvar)
  rhs.deparsed <- terms.deparsed[-1]
  attr(mt, "term.labels") <- rhs.deparsed
  attr(mt, "factors") <- attr(mt, "factors")[seq_along(terms.deparsed),
                                             seq_along(rhs.deparsed),
                                             drop=FALSE]
  dimnames(attr(mt, "factors")) <- list(terms.deparsed, rhs.deparsed)
  vars <- do.call(
    "call",
    c(list("list"), lapply(terms.deparsed, str2lang)),
    quote = TRUE
  )
  attr(mt, "variables") <- vars
  
  # store the indicator to absorb block means in the terms object as well
  attr(mt, "absorb") <- mod_special$absorb

  # set the class
  class(mt) <- c("lmitt.terms", "terms", "formula")
  return(mt)
}

#' Return an \code{lmitt.terms} object with \code{terms()}
#' @exportS3Method stats::terms
#' @param x \code{lmitt.terms}.
#' @param ... Arguments passed to \code{...} do not modify the
#' \code{lmitt.terms} object given.
#' @returns The \code{lmitt.terms} object passed to \code{x} unchanged.
#' @examples
#' data(simdata)
#' spec <- rct_spec(z ~ 1, simdata)
#' form <- build_lmitt_formula("y", "1", rhstype = "intercept", absorb = FALSE)
#' tt <- terms(form, spec)
#' stats::terms(tt)
terms.lmitt.terms <- function(x, ...) {
  return(x)
}

setGeneric("formula")
#' Return the \code{lmitt.formula} associated with a \code{lmitt.terms} object
#' @exportS3Method stats::formula
#' @param x \code{lmitt.terms}.
#' @param ... Arguments passed to \code{...} do not modify the
#' \code{lmitt.terms} object given.
#' @returns The formula associated with the \code{lmitt.terms} object; in other
#' words, the \code{lmitt.formula} it was created from (as a standard
#' \code{formula} class), rather than the \code{model_form} attribute.
#' @examples
#' data(simdata)
#' spec <- rct_spec(z ~ 1, simdata)
#' form <- build_lmitt_formula("y", "1", rhstype = "intercept", absorb = FALSE)
#' tt <- terms(form, spec)
#' stats::formula(tt)
formula.lmitt.terms <- function(x, ...) {
  mc <- match.call()
  mc[[1L]] <- quote(stats:::formula.terms)
  eval(mc, parent.frame(1L))
}