#' @importFrom stats summary.lm residuals
bread.lm <- function(x, ...) {
  if (!is.null(x$na.action))
    class(x$na.action) <- "exclude"
  sx <- stats::summary.lm(x)
  n <- length(residuals(x))
  return(sx$cov.unscaled * n)
}


#' @importFrom stats residuals weights
bread.glm <- function(x, ...) {
  if (!is.null(x$na.action))
    class(x$na.action) <- "exclude"
  sx <- summary(x)
  wres <- as.vector(stats::residuals(x, "working")) * stats::weights(x, "working")
  dispersion <- if (substr(x$family$family, 1L, 17L) %in%
                      c("poisson", "binomial", "Negative Binomial"))
                  1
  else sum(wres^2)/sum(weights(x, "working"))
  return(sx$cov.unscaled * length(sx$deviance.resid) * dispersion)
}

#' @importFrom stats summary.lm residuals
bread.mlm <- function(x, ...) {
  if (!is.null(x$na.action))
    class(x$na.action) <- "exclude"
  K <- x$rank
  piv <- x$qr$pivot[seq_len(K)]
  cf <- x$coef[piv,,drop=FALSE]
  rval <- summary.lm(x)
  n <- nrow(residuals(x))

  rval <- kronecker(
    structure(diag(ncol(cf)), .Dimnames = rep.int(list(colnames(cf)),  2L)),
    structure(rval$cov.unscaled, .Dimnames = rep.int(list(rownames(cf)), 2L)) * n,
    make.dimnames = TRUE)

  return(rval)
}
