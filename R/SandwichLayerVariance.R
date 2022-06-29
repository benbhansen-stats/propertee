#' @include Design.R SandwichLayer.R
NULL

#' @title (Internal) Compute variance blocks
#' @details The \bold{B12 block} is the covariance matrix of the cluster-level
#'   estimating equations for the covariance and direct adjustment models. It
#'   has a row for each term in the covariance model and a column for each term
#'   in the direct adjustment model. For any row that does not appear in both
#'   the experimental design and the covariance model data, its contribution to
#'   this matrix will be 0. Thus, if there is no overlap between the two
#'   datasets, this will return a matrix of 0's.
#' @param x A \code{Lmitted} model
#' @return \code{.get_b12()}: A \eqn{p\times k} matrix where \eqn{p} is the
#'   number of terms in the covariance model and \eqn{k} is the number of terms
#'   in the \code{Lmitted} model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_b12 <- function(x) {
  if (!inherits(x, "Lmitted")) {
    stop("x must be a Lmitted model")
  }

  sl <- x$model$`(offset)`
  if (class(sl) != "SandwichLayer") {
    stop(paste("Lmitted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  uoanames <- var_names(x@Design, "u")
  zname <- var_names(x@Design, "t")

  # Sum est eqns to cluster level; since non-overlapping rows are NA they are
  # excluded automatically in `by` call
  cmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(sl@fitted_covariance_model),
       lapply(uoanames, function(col) sl@keys[, col]),
       colSums))

  # cmod_eqns will be NULL if there's no overlap, so no need to continue
  if (is.null(cmod_eqns)) {
    return(matrix(0,
                  nrow = dim(stats::model.matrix(sl@fitted_covariance_model))[2],
                  ncol = dim(stats::model.matrix(x))[2]))
  }

  # get overlapping rows from experimental data joining with `keys`
  lmitt_uoas <- stats::expand.model.frame(x, uoanames)
  lmitt_uoas <- .merge_preserve_order(
    lmitt_uoas,
    merge(unique(sl@keys), x@Design@structure), # merge here to use txt col for finding NA's
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)

  msk <- !is.na(lmitt_uoas[, paste0(zname, ".y")])
  dmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(x)[msk, ],
       lapply(uoanames, function(col) lmitt_uoas[msk, col]),
       colSums))

  return(crossprod(cmod_eqns, dmod_eqns))
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
#' @return \code{.get_a22_inverse()}: A \eqn{1\times1} matrix
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a22_inverse <- function(x) {
  if (!inherits(x, "Lmitted")) {
    stop("x must be a Lmitted model")
  }

  # Get expected information per sandwich_infrastructure vignette
  if ("glm" %in% x@.S3Class) {
    if (x$family$family %in% c("binomial", "poisson")) {
      dispersion <- 1
    } else {
      dispersion <- sum((x$weights * x$residuals)^2) / sum(x$weights)
    }
    W <- x$family$mu.eta(x$linear.predictors) * x$weights / dispersion
  } else {
    W <- if (is.null(x$weights)) rep(1, length(x$fitted.values)) else x$weights
  }

  fim <- crossprod(stats::model.matrix(x) * W, stats::model.matrix(x))
  zname <- var_names(x@Design, "t")
  out <- solve(fim)[zname, zname, drop = FALSE]

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
#' @return \code{.get_b22()}: A \eqn{(p+1)\times(p+1)} matrix where the
#'   dimensions are given by the number of terms in the \code{Lmitted} model
#'   (\eqn{p}) and an Intercept term. This should be a \eqn{2\times2} matrix
#'   given the dichotomous handling of treatment variables in this package and
#'   the use of the covariance model to offer the covariance adjustment.
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_b22 <- function(x, ...) {
  if (!inherits(x, "Lmitted")) {
    stop("x must be a Lmitted model")
  }

  nq <- sum(summary(x)$df[1L:2L])

  # Get units of assignment for clustering
  uoanames <- var_names(x@Design, "u")
  uoas <- stats::expand.model.frame(x, uoanames)[, uoanames, drop = FALSE]
  if (ncol(uoas) == 1) {
    uoas <- factor(uoas[,1])
  } else {
    uoas <- factor(Reduce(function(...) paste(..., sep = "_"), uoas))
  }

  zname <- var_names(x@Design, "t")
  vmat <- sandwich::meatCL(x, cluster = uoas, ...) * nq
  out <- vmat[zname, zname, drop = FALSE]

  return(out)
}

#' @details The \bold{A11 block} is the \eqn{p\times p} matrix corresponding to
#'   the unscaled inverse of the observed Fisher information of the covariance
#'   model. The observed information is given by the estimate of the negative
#'   Jacobian of the model's estimating equations. The unscaled version provided
#'   here divides by the number of observations used to fit the covariance
#'   model.
#' @return \code{.get_a11_inverse()}: A \eqn{(p+1)\times(p+1)} matrix where the
#'   dimensions are given by the number of terms in the covariance model
#'   (\eqn{p}) and an Intercept term.
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a11_inverse <- function(x) {
  if (!inherits(x, "Lmitted")) {
    stop("x must be a Lmitted model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("Lmitted model must have an offset of class `SandwichLayer`",
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
#' @return \code{.get_b11()}: A \eqn{(p+1)\times(p+1)} matrix the dimensions are
#'   given by the number of terms in the covariance model (\eqn{p}) and an
#'   Intercept term
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_b11 <- function(x, ...) {
  if (!inherits(x, "Lmitted")) {
    stop("x must be a Lmitted model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("Lmitted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  cmod <- sl@fitted_covariance_model
  nc <- sum(summary(cmod)$df[1L:2L])

  # Get units of assignment for clustering
  if (ncol(sl@keys) == 1) {
    uoas <- sl@keys[, 1]
  } else {
    uoas <- Reduce(function(...) paste(..., sep = "_"), sl@keys)
  }

  # Replace NA's for rows not in the experimental design with a unique cluster ID
  nuoas <- length(unique(uoas))
  nas <- grepl("NA", uoas)
  uoas[nas] <- paste0(nuoas - 1 + seq_len(sum(nas)), "*")
  uoas <- factor(uoas)

  out <- sandwich::meatCL(cmod, cluster = uoas, ...) * nc
  return(out)
}
