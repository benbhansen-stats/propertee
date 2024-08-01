#' @importFrom stats summary.lm
bread.lm <- function(x, ...) {
  if (!is.null(x$na.action))
    class(x$na.action) <- "omit"
  sx <- stats::summary.lm(x)
  return(sx$cov.unscaled * length(sx$residuals))
}


#' @importFrom stats residuals weights
bread.glm <- function(x, ...) {
  if (!is.null(x$na.action))
    class(x$na.action) <- "omit"
  sx <- summary(x)
  wres <- as.vector(stats::residuals(x, "working")) * stats::weights(x, "working")
  dispersion <- if (substr(x$family$family, 1L, 17L) %in%
                      c("poisson", "binomial", "Negative Binomial"))
                  1
  else sum(wres^2)/sum(weights(x, "working"))
  return(sx$cov.unscaled * length(sx$deviance.resid) * dispersion)
}
