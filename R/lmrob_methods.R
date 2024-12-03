##' @title Generate matrix of estimating equations for \code{lmrob()} fit
##' @details This is part of a workaround for an issue in the robustbase code
##'   affecting sandwich covariance estimation. The issue in question is issue
##'   #6471, robustbase project on R-Forge. This function contributes to
##'   providing sandwich estimates of covariance-adjusted standard errors for
##'   robust linear covariance adjustment models.
##' @param x An \code{lmrob} object produced using an MM/SM estimator chain
##' @param ... Additional arguments to be passed to \code{estfun}
##' @return A \eqn{n\times }(p+1) matrix where the first column corresponds to
##'   the scale estimate and the remaining \eqn{p} colums correspond to the
##'   coefficients
##' @author Ben B. Hansen
##' @rdname lmrob_methods
##' @exportS3Method
estfun.lmrob <- function(x, ...) {
  ctrl <- x$control
  if (!is(ctrl, "lmrobCtrl") & !inherits(ctrl, "list")) {
    stop("Model object must have a `control` element of type `lmrobCtrl` (robustbase 0.99-0 and newer) or `list`")
  }
  if (!(ctrl$method %in% c("SM", "MM"))) {
    stop("estfun.lmrob() supports only SM or MM estimates")
  }
  if (is.null(ctrl$psi)) {
    stop("parameter psi is not defined")
  }

  xmat <- stats::model.matrix(x)
  xmat <- stats::naresid(x$na.action, xmat)
  psi <- chi <- ctrl$psi
  stopifnot(is.numeric(c.chi <- ctrl$tuning.chi),
            is.numeric(c.psi <- ctrl$tuning.psi))
  r0 <- x$init$resid
  r <- x$residuals
  scale <- x$scale
  n <- length(r)
  stopifnot(n == length(r0), is.matrix(xmat), n == nrow(xmat))
  p <- ncol(xmat)
  r0.s <- r0 / scale
  w0 <- robustbase::Mchi(r0.s, cc = c.chi, psi = chi)
  Usigma <- scale(w0, center=TRUE, scale=FALSE)
  colnames(Usigma) <- "sigma"
  r.s <- r / scale
  w <- robustbase::Mpsi(r.s, cc = c.psi, psi = psi)
  Ubeta <- w * xmat
  rval <- cbind(Usigma, Ubeta)
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL

  return(rval)
}

##' @title Extract bread matrix from an \code{lmrob()} fit
##' @details This is part of a workaround for an issue in the robustbase code
##'   affecting sandwich covariance estimation. The issue in question is issue
##'   #6471, robustbase project on R-Forge. This function contributes to
##'   providing sandwich estimates of covariance-adjusted standard errors for
##'   robust linear covariance adjustment models.
##'
##' @param x An \code{lmrob} object produced using an MM/SM estimator chain
##' @param ... Additional arguments to be passed to \code{bread}
##' @return A \eqn{p\times }(p+1) matrix where the first column corresponds to
##'   the scale estimate and the remaining \eqn{p} colums correspond to the
##'   coefficients
##' @author Ben B. Hansen
##' @rdname lmrob_methods
##' @exportS3Method
bread.lmrob <- function(x, ...) {
  ctrl <- x$control
  if (!is(ctrl, "lmrobCtrl") & !inherits(ctrl, "list")) {
    stop("Model object must have a `control` element of type `lmrobCtrl` (robustbase 0.99-0 and newer) or `list`")
  }
  if (!(ctrl$method %in% c("SM", "MM"))) {
    stop("estfun.lmrob() supports only SM or MM estimates")
  }
  if (is.null(ctrl$psi)) {
    stop("parameter psi is not defined")
  }

  psi <- chi <- ctrl$psi
  stopifnot(is.numeric(c.chi <- ctrl$tuning.chi),
            is.numeric(c.psi <- ctrl$tuning.psi))
  r0 <- x$init$resid
  r <- x$residuals
  scale <- x$scale
  xmat <- stats::model.matrix(x)
  bb <- 1 / 2
  n <- length(r)
  stopifnot(n == length(r0), is.matrix(xmat), n == nrow(xmat))
  p <- ncol(xmat)
  r.s <- r / scale
  r0.s <- r0 / scale
  w <- robustbase::Mpsi(r.s, cc = c.psi, psi = psi, deriv = 1)
  w0 <- robustbase::Mchi(r0.s, cc = c.chi, psi = chi, deriv = 1)
  x.wx <- crossprod(xmat, xmat * w)
  A <- tryCatch(solve(x.wx) * scale, error = function(e) {
    tryCatch({
      out <- solve(x.wx, tol = 0) * scale
      warning("X'WX is almost singular", call. = FALSE)
      out
    }, error = function(e) {
      stop("X'WX is singular", call. = FALSE)
    })
  })

  # At this point A has no sample size scaling, as in robustbase:::.vcov.avar1
  # The lack of scaling there precisely compensates for the lack of scaling of
  # the crossproduct
  a <- A %*% (crossprod(xmat, w * r.s) / mean(w0 * r0.s))
  colnames(a) <- "sigma"

  # Now we restore sample size scaling to A
  A <- n * A
  rval <- cbind(a, A)

  return(rval)
}
