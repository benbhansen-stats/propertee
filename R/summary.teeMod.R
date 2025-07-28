##' @title Summarizing \code{teeMod} objects
##'
##' @description [summary()] method for class \code{teeMod}
##'
##' @details If a \code{teeMod} object is fit with a \code{SandwichLayer}
##'   offset, then the usual \code{stats::summary.lm()} output is enhanced by
##'   the use of covariance-adjusted sandwich standard errors, with t-test
##'   values recalculated to reflect the new standard errors.
##'
##' @param object \code{teeMod} object
##' @param vcov.type A string indicating the desired variance estimator. See
##'   [vcov_tee()] for details on accepted types.
##' @param ... Additional arguments to [vcov_tee()], such as the desired finite
##'   sample heteroskedasticity-robust standard error adjustment.
##' @return object of class \code{summary.teeMod}
##' @export
##' @method summary teeMod
##' @rdname teeMod_summary
summary.teeMod <- function(object,
                                   vcov.type = "HC0",
                                   ...) {
  orig.coefficients <- object$coefficients
  object$coefficients <- replace(orig.coefficients, is.na(orig.coefficients), 0)
  out <- summary(as(object, "lm"))

  if (object$rank > 0) {
    dots <- list(...)
    toprint <- !grepl("^offset:", names(orig.coefficients))
    out$coefficients <- rbind(out$coef, matrix(NA, nrow = sum(toprint) - nrow(out$coef), ncol = 4))
    rownames(out$coefficients) <- names(orig.coefficients)[toprint]

    out$coefficients[, 1L] <- orig.coefficients[toprint]
    covmat <- vcov_tee(object, type = vcov.type, ...)
    dof <- vapply(seq_len(object$rank),
                  .get_dof,
                  numeric(1L),
                  x = object, vcov_type = vcov.type,
                  cls = .make_uoa_ids(object, substr(vcov.type, 1, 2), dots$cluster), ...)
    out$coefficients[!is.na(orig.coefficients), 2L] <- sqrt(diag(covmat))[
      names(orig.coefficients)[toprint & !is.na(orig.coefficients)]]
    out$coefficients[, 3L] <- out$coefficients[, 1L] / out$coefficients[, 2L]
    out$coefficients[, 4L] <- 2*stats::pt(abs(out$coefficients[, 3L]),
                                          c(dof, rep(.get_dof(object, "CR0", 0),
                                                     sum(toprint) - nrow(out$coef))),
                                          lower.tail = FALSE)
    out$vcov.type <- attr(covmat, "type")
    out$vcov.cov_adj_bias_correction <- attr(covmat, "cov_adj_correction")
  }
  class(out) <- "summary.teeMod"
  out$teeMod <- object
  return(out)
}


##' @param x \code{summary.teeMod} object
##' @param digits the number of significant digits to use when printing.
##' @param signif.stars logical. If ‘TRUE’, ‘significance stars’ are printed for
##'   each coefficient.
##' @importFrom stats pt printCoefmat
##' @export
##' @rdname teeMod_summary
print.summary.teeMod <- function(x,
                                         digits =
                                           max(3L, getOption("digits") - 3L),
                                         signif.stars =
                                           getOption("show.signif.stars"),
                                         ...) {

  cat("\nCall:\n", paste(deparse1(x$teeMod@lmitt_call), sep = "\n",
                         collapse = "\n"), "\n", sep = "")

  df <- x$df

  if (x$teeMod@lmitt_fitted) {
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
    if (x$teeMod@lmitt_fitted) {
      toprint <- (
        grepl(paste0("^\\`?", var_names(x$teeMod@StudySpecification, "t"), "\\."), rownames(coefs)) |
          grepl(paste0("^", formula(x$teeMod)[[2L]], ":"), rownames(coefs))
      )
    } else {
      toprint <- !grepl("^offset:", rownames(coefs))
    }
    stats::printCoefmat(coefs[toprint, , drop = FALSE],
                        digits = digits,
                        signif.stars = signif.stars,
                        na.print = "NA", ...)
  }

  if (sum(toprint) > 0 & any(!is.na(coefs[toprint, 1]))) {
    # Only print if we estimate at least one treatment effect
    cat(paste0("Std. Error calculated via type \"", x$vcov.type, "\"\n\n"))
    if (!is.null(bc <- x$vcov.cov_adj_bias_correction)) {
      cat(paste0("Residuals from covariance adjustment model adjusted by an \"",
                 bc, "\" bias correction\n\n"))
    }
  } else {
    cat("\n")
  }
  invisible(x)
}
