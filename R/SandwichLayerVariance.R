#' @include Design.R SandwichLayer.R
NULL

#' @title Compute directly adjusted cluster-robust sandwich variance estimates
#' @param x A \code{DirectAdjusted} model
#' @param type A string indicating the desired variance estimator. Currently
#' accepts "CR1"
#' @param ... Arguments to be passed to the internal variance estimation function
#' @export
#' @rdname var_estimators
vcovDA <- function(x, type = c("CR1"), ...) {
  type <- match.arg(type)

  var_func <- switch(
    type,
    "CR1" = .vcovMB_CR1
  )
  
  est <- var_func(x, ...)
  return(est)
}

#' @title Compute directly-adjusted CR1 cluster-robust sandwich variance estimates
#' under model-based assumptions
#' @param x A \code{DirectAdjusted} model
#' @param ... Arguments to be passed to sandwich::meatCL
#' @return \code{.vcovMB_CR1()}: A \eqn{2\times 2} matrix where the dimensions are
#' given by the intercept and treatment variable terms in the direct adjustment
#' model
#' @keywords internal
#' @rdname var_estimators
.vcovMB_CR1 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  # compute blocks
  a21 <- .get_a21(x)
  a11inv <- .get_a11_inverse(x)
  b12 <- .get_b12(x)

  a22inv <- .get_a22_inverse(x)
  b22 <- .get_b22(x, ...)
  b11 <- .get_b11(x, ...)

  meat <- (
    b22 -
      a21 %*% a11inv %*% b12 -
      t(b12) %*% t(a11inv) %*% t(a21) +
      a21 %*% a11inv %*% b11 %*% t(a11inv) %*% t(a21)
  )
  vmat <- a22inv %*% meat %*% a22inv

  return(vmat)
}

#' @title (Internal) Compute variance blocks
#' @details The \bold{B12 block} is the covariance matrix of the cluster-level
#'   estimating equations for the covariance and direct adjustment models. It
#'   has a row for each term in the covariance model and a column for each term
#'   in the direct adjustment model. For any row that does not appear in both
#'   the experimental design and the covariance model data, its contribution to
#'   this matrix will be 0. Thus, if there is no overlap between the two
#'   datasets, this will return a matrix of 0's.
#' @param x A \code{DirectAdjusted} model
#' @return \code{.get_b12()}: A \eqn{p\times 2} matrix where the number of rows
#'   are given by the number of terms in the covariance model and the number of
#'   columns correspond to intercept and treatment variable terms in the direct
#'   adjustment model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_b12 <- function(x) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  sl <- x$model$`(offset)`
  if (!is(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  uoanames <- var_names(x@Design, "u")
  zname <- var_names(x@Design, "t")

  # Check number of overlapping clusters n_QC; if n_QC <= 1, return 0
  if (ncol(sl@keys) == 1) {
    keys <- sl@keys[, 1]
  } else {
    keys <- Reduce(function(...) paste(..., sep = "_"), sl@keys)
    keys[grepl("NA", keys)] <- NA_integer_
  }
  
  no_uoas_overlap <- all(is.na(unique(keys)))
  if (no_uoas_overlap) {
    return(matrix(0,
                  nrow = dim(stats::model.matrix(sl@fitted_covariance_model))[2],
                  ncol = dim(stats::model.matrix(x))[2]))
  }
  
  # Sum est eqns to cluster level; since non-overlapping rows are NA in `keys`,
  # `by` call excludes them from being summed
  cmod_estfun <- sandwich::estfun(sl@fitted_covariance_model)
  cmod_aggfun <- ifelse(dim(cmod_estfun)[2] > 1, colSums, sum)
  cmod_eqns <- Reduce(rbind, by(cmod_estfun, keys, cmod_aggfun))

  # get rows from overlapping clusters in experimental data
  Q_uoas <- stats::expand.model.frame(x, uoanames, na.expand = TRUE)[uoanames]
  if (ncol(Q_uoas) == 1) {
    Q_uoas <- Q_uoas[, 1]
  } else {
    Q_uoas <- Reduce(function(...) paste(..., sep = "_"), Q_uoas)
    Q_uoas[grepl("NA", Q_uoas)] <- NA_integer_
  }

  msk <- Q_uoas %in% unique(keys[!is.na(keys)])
  damod_estfun <- sandwich::estfun(x)[msk, , drop = FALSE]
  damod_aggfun <- ifelse(dim(damod_estfun)[2] > 1, colSums, sum)
  damod_eqns <- Reduce(rbind, by(damod_estfun, Q_uoas[msk], damod_aggfun))

  matmul_func <- if (length(unique(keys[!is.na(keys)])) == 1) tcrossprod else crossprod
  return(matmul_func(cmod_eqns, damod_eqns))
}

#' @details The \bold{A22 block} is the diagonal element of the inverse expected
#'   Fisher Information matrix corresponding to the treatment estimate. As shown
#'   in the Details of \code{.get_b22()}, the estimating equations for a
#'   generalized linear model with a canonical link function can be written as
#'   \deqn{\psi_i = (r_i / Var(y_i)) * (d\mu_i/d\eta_i) * x_i} The expected
#'   information matrix A22 is then the negative Jacobian of \eqn{\psi_i}, which
#'   by likelihood theory is the variance-covariance matrix of \eqn{\psi_i}:
#'   \deqn{\psi_{i}\psi_{i}^{T} = E[(r_i / Var(y_i) * (d\mu_i/d\eta_i))^{2}] *
#'   x_ix_i' = (d\mu_i/d\eta_i)^{2} / Var(y_i) * x_ix_i'} Considering the whole
#'   sample, this can be expressed in matrix form as \eqn{X'WX}. The output of
#'   this function is the inverse of the diagonal element corresponding to the
#'   treatment estimate.
#' @return \code{.get_a22_inverse()}: A \eqn{2\times 2} matrix where the
#' dimensions are given by the intercept and treatment variable terms in the
#' direct adjustment model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a22_inverse <- function(x) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  # Get expected information per sandwich_infrastructure vignette
  w <- if (is.null(x$weights)) 1 else x$weights
  out <- solve(crossprod(stats::model.matrix(x) * sqrt(w)))

  return(out)
}

#' @param ... Arguments to be passed to sandwich::meatCL
#' @details The \bold{B22 block} refers to a clustered variance estimate of the
#'   treatment effect estimate. The \code{stats} package offers family objects
#'   with canonical link functions, so the log-likelihood for a generalized
#'   linear model can be written in terms of the linear predictor as
#'   \deqn{L(y_i, \beta, \phi, w_i) = w_i * (y_i * \beta'x_i - b(\beta'x_i)) /
#'   \phi + h(y_i; \phi)} The estimating equations \eqn{\psi_i} given by the
#'   score function can then be expressed as \deqn{\psi_i = E[w_i * (y_i -
#'   \mu(\beta'x_i)) * x_i / \phi]} In section 4.4 of the second edition of
#'   Categorical Data Analysis, Agresti shows the derivative of the mean
#'   function with respect to the linear predictor is equivalent to the weighted
#'   variance for an observation divided by the estimate of the dispersion
#'   parameter. Thus, the above estimating equations can also be written as
#'   \deqn{\psi_i = (r_i / Var(y_i)) * (d\mu_i/d\eta_i) * x_i}
#'
#'   A matrix \eqn{C} of dimension \eqn{n\times J} is formed to indicate which
#'   unit each subject in the design belongs to, where \eqn{J} is the number of
#'   units at the level of treatment assignment. The treatment assignment-level
#'   estimating equations are then obtained via \eqn{C'\psi}, where \eqn{\psi}
#'   is the matrix of estimating equations at the subject level.
#' @references Agresti, Alan. Categorical Data Analysis. 2003. Open WorldCat,
#'   https://nbn-resolving.org/urn:nbn:de:101:1-201502241089.
#' @return \code{.get_b22()}: A \eqn{2\times 2} matrix where the
#'   dimensions are given by the intercept and treatment variable terms in the
#'   direct adjustment model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_b22 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  nq <- sum(summary(x)$df[1L:2L])

  # Get units of assignment for clustering
  dots <- list(...)
  if (is.null(dots$cluster)) {
    uoanames <- var_names(x@Design, "u")
    uoas <- stats::expand.model.frame(x, uoanames)[, uoanames, drop = FALSE]
    if (ncol(uoas) == 1) {
      uoas <- factor(uoas[,1])
    } else {
      uoas <- factor(Reduce(function(...) paste(..., sep = "_"), uoas))
    }
    dots$cluster <- uoas
  }
  dots$x <- x

  out <- do.call(sandwich::meatCL, dots) * nq

  return(out)
}

#' @details The \bold{A11 block} is the \eqn{p\times p} matrix corresponding to
#'   the unscaled inverse of the observed Fisher information of the covariance
#'   model. The observed information is given by the estimate of the negative
#'   Jacobian of the model's estimating equations. The unscaled version provided
#'   here divides by the number of observations used to fit the covariance
#'   model.
#' @return \code{.get_a11_inverse()}: A \eqn{p\times p} matrix where the
#'   dimensions are given by the number of terms in the covariance model
#'   including an intercept
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a11_inverse <- function(x) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  cmod <- sl@fitted_covariance_model
  nc <- sum(summary(cmod)$df[1L:2L])

  out <- sandwich::bread(cmod) / nc
  return(out)
}

#' @details The \bold{B11 block} is the block of the sandwich variance estimator
#'   corresponding to the variance-covariance matrix of the covariance model
#'   coefficient estimates. The estimates returned here are potentially
#'   clustered (by the clustering in the experimental design) if rows in the
#'   covariance model data also exist in the design. If there is no overlap
#'   between the two datasets, the variance-covariance matrix is estimated
#'   assuming the observations are independent.
#' @return \code{.get_b11()}: A \eqn{p\times p} matrix where the dimensions are
#'   given by the number of terms in the covariance model including an
#'   intercept
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_b11 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  cmod <- sl@fitted_covariance_model
  nc <- sum(summary(cmod)$df[1L:2L])

  # Get units of assignment for clustering
  dots <- list(...)
  if (is.null(dots$cluster)) {
    if (ncol(sl@keys) == 1) {
      uoas <- sl@keys[, 1]
    } else {
      uoas <- Reduce(function(...) paste(..., sep = "_"), sl@keys)
      uoas[grepl("NA", uoas)] <- NA_integer_
    }

    # Replace NA's for rows not in the experimental design with a unique cluster ID
    nuoas <- length(unique(uoas))
    nas <- is.na(uoas)
    uoas[nas] <- paste0(nuoas - 1 + seq_len(sum(nas)), "*")
    uoas <- factor(uoas)
    nuoas <- length(unique(uoas))
    
    # if the covariance model only uses one cluster, treat units as independent
    if (nuoas == 1) {
      uoas <- seq_len(length(uoas))
    }
    dots$cluster <- uoas
  }
  dots$x <- cmod

  out <- do.call(sandwich::meatCL, dots) * nc
  return(out)
}


#' @details The \bold{A21 block} is the block of the sandwich variance estimator
#'   corresponding to the gradient of the direct adjustment model with respect
#'   to the covariates. Some of the information needed for this calculation is
#'   stored in the \code{DirectAdjusted} object's \code{SandwichLayer} offset. This
#'   block is the crossproduct of the prediction gradient and the gradient of
#'   the conditional mean vector for the direct adjustment model summed to the
#'   cluster level. In other words, we take this matrix to be \deqn{\sum(d\psi_i
#'   / d\alpha) = -\sum(w_i/\phi) * (d\mu(\eta_i) / d\eta_i) *
#'   (d\upsilon(\zeta_i) / d\zeta_i) * (x_i c_i)x_i'} where \eqn{\mu} and
#'   \eqn{\eta_i} are the conditional mean function and linear predictor for the
#'   ith cluster in the direct adjustment model, and \eqn{\upsilon} and
#'   \eqn{\zeta_i} are the conditional mean function and linear predictor for
#'   the ith cluster in the covariance model.
#' @return \code{.get_a12()}: A \eqn{2\times p} matrix where the number of
#'   rows are given by intercept and treatment variable terms in the direct
#'   adjusted model, and the number of columns are given by the number of terms
#'   in the covariance model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a21 <- function(x) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  # Get contribution to the estimating equation from the direct adjustment model
  w <- if (is.null(x$weights)) 1 else x$weights

  damod_mm <- stats::model.matrix(formula(x),
                                  stats::model.frame(x, na.action = na.pass))
  msk <- (apply(!is.na(sl@prediction_gradient), 1, all) &
            apply(!is.na(damod_mm), 1, all))
  
  out <- crossprod(damod_mm[msk, , drop = FALSE] * w,
                   sl@prediction_gradient[msk, , drop = FALSE])

  return(out)
}

##' @title Generate matrix of estimating equations for \code{lmrob()} fit
##' @details This is part of a workaround for an issue in the robustbase code
##' affecting sandwich covariance estimation. The issue in question is issue
##' #6471, robustbase project on R-Forge. This function contributes to providing
##' sandwich estimates of covariance-adjusted standard errors for robust linear
##' covariance adjustment models.
##' @param x An \code{lmrob} object produced using an MM/SM estimator chain
##' @param ... Additional arguments to be passed to \code{estfun}
##' @return A \eqn{n\times }(p+1) matrix where the first column corresponds to
##' the scale estimate and the remaining \eqn{p} colums correspond to the
##' coefficients
##' @author lrd author 2
##' @rdname lmrob_methods
##' @exportS3Method
estfun.lmrob <- function(x, ...) {
  ctrl <- x$control
  if (!inherits(ctrl, "list")) {
    stop("Model object must have a `control` element of type `list`")
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
##' affecting sandwich covariance estimation. The issue in question is issue
##' #6471, robustbase project on R-Forge. This function contributes to providing
##' sandwich estimates of covariance-adjusted standard errors for robust linear
##' covariance adjustment models.
##'
##' @param x An \code{lmrob} object produced using an MM/SM estimator chain
##' @param ... Additional arguments to be passed to \code{bread}
##' @return A \eqn{p\times }(p+1) matrix where the first column corresponds to
##' the scale estimate and the remaining \eqn{p} colums correspond to the
##' coefficients
##' @author lrd author 2
##' @rdname lmrob_methods
##' @exportS3Method
bread.lmrob <- function(x, ...) {
  ctrl <- x$control
  if (!inherits(ctrl, "list")) {
    stop("Model object must have a `control` element of type `list`")
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
