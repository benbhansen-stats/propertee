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
#' linear model with a canonical link function can be written as
#' \deqn{\psii = (ri / Var(yi)) * (d\mui/d\etai) * xi}
#' The expected information matrix A22 is then the negative Jacobian of \eqn{\psii},
#' which by likelihood theory is the variance-covariance matrix of \eqn{\psii}:
#' \deqn{\psi_{i}\psi_{i}^{T} = E[(ri / Var(yi) * (d\mui/d\etai))^{2}] * xixi' =
#' (d\mui/d\etai)^{2} / Var(yi) * xixi'} Considering the whole sample, this can be
#' expressed in matrix form as \eqn{X'WX}. The output of this function is the
#' inverse of the diagonal element corresponding to the treatment estimate.
#' @return A 1x1 matrix
#' @keywords internal
.get_a22_inverse <- function(x) {
  if (!is(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  # Get expected information per sandwich_infrastructure vignette
  if ("glm" %in% x@.S3Class) {
    varys <- x$family$variance(x$fitted.values)
    if (substr(x$family$family, 1, 5) == "quasi") {
      varys <- stats::summary.glm(x)$dispersion * varys
    }
    W <- x$family$mu.eta(x$linear.predictors)^2 / varys
  } else {
    W <- rep(1, length(x$fitted.values))
  }
  
  fim <- crossprod(stats::model.matrix(x) * W, stats::model.matrix(x))
  zname <- var_names(x@Design, "t")
  
  return(solve(fim)[zname, zname, drop = FALSE])
}

#' @title (Internal) Get the B22 block of the sandwich variance estimator
#' @param x A DirectAdjusted model object
#' @param ... Arguments to be passed to sandwich::meatCL
#' @details This block refers to a clustered variance estimate of the treatment
#' effect estimate. The \code{stats} package offers family objects with
#' canonical link functions, so the log-likelihood for a generalized linear model
#' can be written in terms of the linear predictor as \deqn{L(yi, \beta, \phi, wi)
#' = wi * (yi * \beta'xi - b(\beta'xi)) / \phi + h(yi; \phi)} The estimating equations
#' \eqn{\psii} given by the score function can then be expressed as
#' \deqn{\psii = E[wi * (yi - \mu(\beta'xi)) * xi / \phi]}
#' In section 4.4 of the second edition of Categorical Data Analysis, Agresti
#' shows the derivative of the mean function with respect to the linear predictor
#' is equivalent to the weighted variance for an observation divided by the estimate
#' of the dispersion parameter. Thus, the above estimating equations can also be
#' written as \deqn{\psii = (ri / Var(yi)) * (d\mui/d\etai) * xi}
#' 
#' A matrix C of dimension nxJ is formed to indicate which unit each subject in
#' the design belongs to, where J is the number of units at the level of treatment
#' assignment. The treatment assignment-level estimating equations are then obtained
#' via \eqn{C'\psi}, where \eqn{\psi} is the matrix of estimating equations at the
#' subject level.
#' @references Agresti, Alan. Categorical Data Analysis. 2003. Open WorldCat,
#' https://nbn-resolving.org/urn:nbn:de:101:1-201502241089.
#' @return A (p+1)x(p+1) matrix where the dimensions are given by the number of
#' terms in the treatment model (p) and an Intercept term. Typically, this should
#' be a 2x2 matrix given the dichotomous handling of treatment variables in this
#' package and the use of the covariance model to offer the covariance adjustment. 
#' @keywords internal
.get_b22 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
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
  out <- sandwich::meatCL(x, cluster = uoas, ...)[zname, zname, drop = FALSE] * nq

  return(out)
}

#' @title (Internal) Get the inverse of the A11 block of the sandwich variance estimator
#' @param x A DirectAdjusted model object
#' @details This block is the pxp matrix corresponding to the unscaled
#' inverse of the observed Fisher information of the covariance model. The
#' observed information is given by the estimate of the negative Jacobian of the
#' model's estimating equations. The unscaled version provided here divides by
#' the number of observations used to fit the covariance model.
#' @return a (p+1)x(p+1) matrix where the dimensions are given by the number of
#' terms in the covariance model (p) and an Intercept term.
#' @keywords internal
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

  return(sandwich::bread(cmod) / nc)
}

#' (Internal) Get the B11 block of the sandwich variance estimator
#' @param x A DirectAdjusted model object
#' @param ... Arguments to be passed to sandwich::meatCL
#' @details This is the block of the sandwich variance estimator corresponding to
#' the variance-covariance matrix of the covariance model coefficient estimates.
#' The estimates returned here are potentially clustered (by the clustering in the
#' experimental design) if rows in the covariance model data also exist in the
#' design. If there is no overlap between the two datasets, the variance-covariance
#' matrix is estimated assuming the observations are independent.
#' @return A (p+1)x(p+1) matrix the dimensions are given by the number of
#' terms in the covariance model (p) and an Intercept term
#' @keywords internal
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
  if (ncol(sl@keys) == 1) {
    uoas <- sl@keys[,1]
  } else {
    uoas <- Reduce(function(...) paste(..., sep = "_"), sl@keys)
  }

  nuoas <- length(unique(uoas))
  nas <- grepl("NA", uoas)
  uoas[nas] <- paste0(nuoas - 1 + seq_len(sum(nas)), "*")
  uoas <- factor(uoas)
  nuoas <- length(levels(uoas))

  return(sandwich::meatCL(cmod, cluster = uoas, ...) * nc)
}
