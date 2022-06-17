##' @title Linear Model for Intention To Treat
##' @param formula See \code{lm()}
##' @param design Optional, explicitly specify the \code{Design} to be used. If
##'   the \code{Design} is specified elsewhere in the model (e.g. passed as an
##'   argument to any of \code{ate()}, \code{ett()}, \code{cov_adj()} or
##'   \code{adopters()}) it will be found automatically and does not need to be
##'   passed here as well. (If the \code{Design} is found in the model, this
##'   argument is ignored.)
##' @param target Optional, explicitly specify the estimand. If the
##'   \code{weights} are generated using either \code{ate()} or \code{ett()},
##'   the \code{target]} will be found automatically. Otherwise, specify whether
##'   the goal is estimating ATE ("ate") or ETT ("ett"). (If weights are
##'   specified, this argument is ignored.)
##' @param ... Additional arguments passed to \code{lm()}.
##' @return \code{DirectAdjusted} model.
##' @export
##' @importFrom stats lm predict weights
lmitt <- function(formula,
                 design = NULL,
                 target = NULL,
                 ...) {

  mf <- match.call()

  no_adopters <- is.null(attr(terms(formula, specials = "adopters"),
                                "specials")$adopters)
  if (no_adopters) {
    formula <- update(formula, . ~ . : adopters())
  }

  if (is(design, "formula")) {
    # If there's a `forcing()`, user wants RDD. If not, force Obs. To do RCT,
    # must create Design manually.
    if (!is.null(attr(terms(design, specials = "forcing"),
                      "specials")$forcing)) {
      des_call <- "rd_design"
    } else {
      des_call <- "obs_design"
    }

    new_d_call <- call(des_call,
                       form = design,
                       data = mf$data,
                       subset = mf$subset,
                       dichotomy = mf$dichotomy)
    design <- eval(new_d_call)
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

  return(as.DirectAdjusted(model, design, target))

}
