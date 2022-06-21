#' @include Design.R SandwichLayer.R
NULL

#' @title (Internal) Get the B12 block of the sandwich variance estimator
#' @details This block is the covariance matrix of the cluster-level estimating
#' equations for the covariance and direct adjustment models. It has a row for
#' each term in the covariance model and a column for each term in the direct
#' adjustment model. For any row that does not appear in both the experimental
#' design and the covariance model data, its contribution to this matrix will be
#' 0. Thus, if there is no overlap between the two datasets, this will return a
#' matrix of 0's.
#' @param x A DirectAdjusted model object
#' @return A pxk matrix where p is the column count in the covariance model
#' matrix and k is the column count in the treatment model matrix
#' @keywords internal
.get_b12 <- function(x) {
  sl <- x$model$`(offset)`
  if (class(sl) != "SandwichLayer") {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
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
  .data <- stats::expand.model.frame(x, uoanames)
  dmod_data <- .merge_preserve_order(
    .data,
    merge(unique(sl@keys), x@Design@structure), # merge here to use txt col for finding NA's
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)
  
  msk <- !is.na(dmod_data[, paste0(zname, ".y")])
  dmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(x)[msk, ],
       lapply(uoanames, function(col) dmod_data[msk, col]),
       colSums))

  return(t(cmod_eqns) %*% dmod_eqns)
}

#' @title (Internal) Get the inverse of the A22 block of the sandwich variance estimator
#' @param x A DirectAdjusted model object
#' @details This block is the diagonal element of the inverse expected Fisher
#' Information matrix corresponding to the treatment estimate. As shown in the
#' Details of \code{\link{.get_b22}}, the estimating equations for a generalized
#' linear model can be written as \deqn{\psi_{i} = w_{i}(y_{i} - \mu(\beta'x_{i}))x_{i}
#' / \phi} The information matrix A22 is then the negative
#' Jacobian of \eqn{\psi_{i}}: \deqn{A_{22} = \partial\psi_{i}/\partial\beta =
#' -w_{i}(d\mu(x_{i}\beta)/dx_{i}\beta)x_{i}x_{i}' / \phi} In
#' matrix form, this can be expressed as \eqn{-X'WX}, where W is a diagonal
#' matrix of the vector of elementwise products of the terms not
#' found in the model matrix.
#' @return A 1x1 matrix
#' @keywords internal
.get_a22_inverse <- function(x) {
  if (!is(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  # Get Fisher information matrix
  if ("glm" %in% x@.S3Class) {
    dispersion <- stats::summary.glm(x)$dispersion
    mu.eta <- x$family$mu.eta(x$linear.predictors)
  } else {
    dispersion <- 1
    mu.eta <- rep(1, length(x$fitted.values))
  }

  W <- diag(as.numeric(-x$weights * mu.eta / dispersion))

  fim <- -(t(stats::model.matrix(x)) %*% W %*% stats::model.matrix(x))
  zname <- var_names(x@Design, "t")
  return(solve(fim)[zname, zname, drop = FALSE])
}

#' @title (Internal) Get the B22 block of the sandwich variance estimator
#' @param x A DirectAdjusted model object
#' @details This block refers to the treatment assignment-level variance estimate
#' of the treatment effect. The \code{stats} package offers family objects with
#' canonical link functions, so the log-likelihood for a generalized linear model
#' can be written in terms of the linear predictor as \deqn{L(y_{i}, \beta, \phi,
#' w_{i}) = w_{i}(y_{i}\beta'x_{i} - b(\beta'x_{i})) / \phi + h(y_{i}; \phi)}
#' The estimating equations \eqn{\psi_{i}} given by the score function can then be
#' expressed as \deqn{\psi_{i} = w_{i}(y_{i} - \mu(\beta'x_{i}))x_{i} / \phi}
#' 
#' Where J is the number of units at the level where the treatment was assigned
#' in the experimental design, a matrix C of dimension nxJ is formed to indicate
#' which unit each subject in the design belongs to. The treatment assignment-
#' level estimating equations are then obtained via \eqn{C'\psi}, where \eqn{\psi}
#' is the matrix of estimating equations at the subject level.
#' @return A (p+1)x(p+1) matrix where the dimensions are given by the number of
#' terms in the treatment model (p) and an Intercept term. Typically, this should
#' be a 2x2 matrix given the dichotomous handling of treatment variables in this
#' package and the use of the covariance model to offer the covariance adjustment. 
#' @keywords internal
.get_b22 <- function(x) {
  if (!is(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  # Get wide table of unit of assignment indicators
  uoanames <- var_names(x@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(x, uoanames)[, uoanames])
  
  # Get unit of assignment-level est. eqns
  # The sandwich package uses a different dispersion calculation from the stats
  # package, so to keep in line with the stats package, we use mostly the same
  # code as sandwich::estfun but swap the dispersion calculation
  if ("glm" %in% x@.S3Class) {
    eqns <- x$weights * x$residuals * stats::model.matrix(x) / stats::summary.glm(x)$dispersion
  } else {
    eqns <- x$weights * x$residuals * stats::model.matrix(x)
  }
  
  cluster_eqns <- crossprod(uoa_matrix, eqns)
  vcov <- crossprod(cluster_eqns)
  zname <- var_names(x@Design, "t")
  return(vcov[zname, zname, drop = FALSE])
}
