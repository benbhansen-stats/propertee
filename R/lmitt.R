##' @title Linear Model for Intention To Treat
##' @param formula See \code{lm()}
##' @param design Optional, explicitly specify the \code{Design} to be used. If
##'   the \code{Design} is specified elsewhere in the model (e.g. passed as an
##'   argument to any of \code{ate()}, \code{ett()}, \code{cov_adj()} or
##'   \code{adopters()}) it will be found automatically and does not need to be
##'   passed here as well. (If different \code{Design} objects are passed
##'   (either through the \code{lm} in weights or covariance adjustment, or
##'   through this argument), an error will be produced.)
##' @param ... Additional arguments passed to \code{lm()}.
##' @return \code{DirectAdjusted} model.
##' @export
##' @importFrom stats lm predict weights
lmitt <- function(formula,
                 design = NULL,
                 ...) {

  mf <- match.call()

  no_adopters <- is.null(attr(terms(formula, specials = "adopters"),
                                "specials")$adopters)
  if (no_adopters) {
    formula <- update(formula, . ~ . : adopters())
  }

  if (inherits(design, "formula")) {
    # If there's a `forcing()`, user wants RDD. If not, force Obs. To do RCT,
    # must create Design manually.
    if (!is.null(attr(terms(design, specials = "forcing"),
                      "specials")$forcing)) {
      des_call <- "rd_design"
    } else {
      des_call <- "obs_design"
    }

    # Build new call. All calls must include formula and data
    new_d_call <- paste0(des_call, "(",
                         "formula = ", deparse(design),
                         ", data = ", deparse(mf$data))
    # If user passed subset or dichotomy, include those. We do this so the
    # `design@call` will be in agreement.
#    if (!is.null(mf$subset)) {
#      new_d_call <- paste0(new_d_call, ", subset = ", deparse(mf$subset))
#    }
    if (!is.null(mf$dichotomy)) {
      new_d_call <- paste0(new_d_call, ", dichotomy = ", deparse(mf$dichotomy))
    }
    new_d_call <- paste0(new_d_call, ")")
    # str2lang converts character into call
    design <- eval(str2lang(new_d_call))
  }


  mf <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "method", "model", "x", "y", "qr", "singular.ok",
               "contrasts", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::lm)
  mf[[2L]] <- formula
  model <- eval(mf, parent.frame())

  model$call[[1]] <- as.name("lmitt")

  return(as.lmitt(model, design))

}
