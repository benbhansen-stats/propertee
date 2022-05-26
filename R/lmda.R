##' @title Direct Adjusted Linear Regression Model
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
##' @param data See \code{lm()}
##' @param subset See \code{lm()}
##' @param weights See \code{lm()}
##' @param na.action See \code{lm()}
##' @param method See \code{lm()}
##' @param model See \code{lm()}
##' @param x See \code{lm()}
##' @param y See \code{lm()}
##' @param qr See \code{lm()}
##' @param singular.ok See \code{lm()}
##' @param contrasts See \code{lm()}
##' @param offset See \code{lm()}
##' @param ... See \code{lm()}
##' @return \code{DirectAdjusted} model.
##' @export
##' @importFrom stats lm predict weights
lmda <- function(formula,
                 design = NULL,
                 target = NULL,
                 data,
                 subset,
                 weights,
                 na.action,
                 method = "qr",
                 model = TRUE,
                 x = FALSE,
                 y = FALSE,
                 qr = TRUE,
                 singular.ok = TRUE,
                 contrasts = NULL,
                 offset,
                 ...) {

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "method", "model", "x", "y", "qr", "singular.ok",
               "contrasts", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::lm)
  model <- eval(mf, parent.frame())

  model$call[[1]] <- as.name("lmda")

  return(as.DirectAdjusted(model, design, target))

}
