#' @include DirectAdjusted.R
NULL
# The above ensures that `DirectAdjusted` is defined prior to
# `DirectAdjustedPBPH`.

setClass("DirectAdjustedPBPH",
         contains = "DirectAdjusted")


##' Peters-Belson with Prognostic Heterogeneity
##'
##' Performs Peters-Belson with Prognostic Heterogeneity on the data.
##' @param mod1 First stage fitted model.
##' @param design \code{Design} object
##' @param data Data where variables in \code{form} live.
##' @return A \code{pbph} object which extends \code{lm}. Can be
##'   passed to \code{summary} or \code{confint}.
##'
##'   The return contains an additional object, \code{pbph}, which \code{lm}
##'   doesn't include. \code{pbph} contains \code{mod1}, a copy of the first
##'   stage model, and \code{data}, which includes the \code{data} augmented
##'   with \code{treatment} and \code{predicted}.
##' @export
lmitt_pbph <- function(mod1, design, data) {

  warning(paste("Note: If first stage model was not fit exclusively on the",
                "control group, results currently may not be valid.\n"))

  # Get predicted values.
  pred <- predict(mod1, type = "response", newdata = data)
  progscore <- pred - mean(pred[treatment(design, newdata = data)[, 1] == 1])

  newform <- update(formula(mod1), . ~ assigned() + progscore:assigned() + 0)
  # Forcing the environment here breaks a lot of the things we've worked on more
  # generally, but is the fastest/easiest way to proceed in this very limited
  # situation.
  environment(newform) <- environment()

  mod2 <- lm(newform,
             data = data,
             offset = cov_adj(mod1, design = design))

  mod2 <- as.lmitt(mod2)



  return(new("DirectAdjustedPBPH", mod2))
}


##' @title Show an DirectAdjustedPBPH
##' @param object DirectAdjustedPBPH object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "DirectAdjustedPBPH", function(object) {
  print(object$coefficients)
  invisible(object)
})

##' @export
summary.DirectAdjustedPBPH <- function(object,
                                   vcov.type = "PBPH",
                                   ...) {
  out <- summary(as(object, "lm"))

  if (object$rank > 0) {
    covmat <- vcovDA(object, type = vcov.type, ...)
    out$coefficients[, 2L] <- sqrt(diag(covmat))
    out$coefficients[, 3L] <- out$coefficients[, 1L] / out$coefficients[, 2L]
    out$coefficients[, 4L] <- 2*stats::pt(abs(out$coefficients[, 3L]),
                                          object$df.residual,
                                          lower.tail = FALSE)
    out$vcov.type <- attr(covmat, "type")
  }

  class(out) <- "summary.DirectAdjustedPBPH"
  out$DirectAdjusted <- object
  return(out)
}


##' @title Print summary of \code{DirectAdjustedPBPH} object
##' @param x \code{summary.DirectAdjustedPBPH} object
##' @param digits the number of significant digits to use when printing.
##' @param signif.stars logical. If ‘TRUE’, ‘significance stars’ are printed for
##'   each coefficient.
##' @param ... Other args
##' @return object, invisibly
##' @importFrom stats pt printCoefmat
##' @export
print.summary.DirectAdjustedPBPH <- function(x,
                                         digits =
                                           max(3L, getOption("digits") - 3L),
                                         signif.stars =
                                           getOption("show.signif.stars"),
                                         ...) {

  df <- x$df

  if (x$DirectAdjusted@lmitt_fitted) {
    coefmatname <- "Treatment Effects"
  } else {
    coefmatname <- "Coefficients"
  }

  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
    return(invisible(x))
  }
  else {
    if (nsingular <- df[3L] - df[1L])
      cat("\n", coefmatname, ": (", nsingular,
          " not defined because of singularities)\n",
          sep = "")
    else cat("\n", coefmatname, ":\n")
    coefs <- x$coefficients
    if (any(aliased <- x$aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4,
                      dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    stats::printCoefmat(coefs,
                        digits = digits,
                        signif.stars = signif.stars,
                        na.print = "NA", ...)
  }

  cat(paste0("Std. Error calculated via type \"", x$vcov.type, "\"\n\n"))
  invisible(x)
}

##' @exportS3Method
vcov.DirectAdjustedPBPH <- function(object, ...) {
  cl <- match.call()

  if (is.null(cl[["type"]])) {
    cl$type <- "PBPH"
  }

  cl$x <- cl$object
  argmatch <- match(c("x", "type", "cluster"), names(cl), nomatch = 0L)
  new_cl <- cl[c(1L, argmatch)]
  new_cl[[1L]] <-  quote(vcovDA)
  vmat <- eval(new_cl, parent.frame())
  return(vmat)
}


.vcov_PBPH <- function(x, ...) {
  camod <- x$model$"(offset)"@fitted_covariance_model

  clstrs_ctrl <- x$model$"(offset)"@keys[,1]
  clstrs <- .make_uoa_ids(x)
  txt <- .bin_txt(x@Design)


  bnm <- .get_bread_and_meat_PBPH(camod, x, txt,
                                  clstrs_ctrl,
                                  clstrs)


  covmat <- (bnm$b22 %*% (bnm$m22 + bnm$b21 %*% bnm$b11 %*% bnm$m11 %*% bnm$b11 %*% t(bnm$b21)) %*% bnm$b22)
  return(covmat)
}

##' @importFrom sandwich bread meatCL
.get_bread_and_meat_PBPH <- function(mod1,
                                     mod2,
                                     treatment,
                                     clusters_mod1,
                                     clusters_mod2) {

  b11 <- sandwich::bread(mod1) / length(residuals(mod1))
  b22 <- sandwich::bread(mod2) / length(residuals(mod2))
  b21 <- bread21(mod2, eta = mod2$coef[2])

  m11 <-  sandwich::meatCL(mod1, clusters_mod1) * length(residuals(mod1))
  m22 <-  sandwich::meatCL(mod2, clusters_mod2) * length(residuals(mod2))

  return(list(b11 = b11,
              b21 = b21,
              b22 = b22,
              m11 = m11,
              m22 = m22))

}

##' @importFrom stats model.response model.matrix
bread21 <- function(model, eta ) {
  camod <- model$model$"(offset)"@fitted_covariance_model
  des <- model@Design

  resp <- stats::model.response(model$model)[treatment(des) == 1]
  covs <- stats::model.matrix(formula(camod), data = model$call$data)
  covs <- covs[treatment(des) == 1, , drop = FALSE]
  pred <- model.frame(model)$progscore[treatment(des) == 1]

  # Replace tauhat with tauhat_eta0
  # When eta = eta_0, the RHS is constant, so we just need the mean.
  tau <- mean(resp - (1 + eta) * pred)
  b21.1 <- apply(-(1 + eta) * covs, 2, sum)
  b21.2 <- apply( (resp - tau - 2 * (1 + eta) * pred) * covs, 2, sum)

  rbind(b21.1, b21.2)
}
