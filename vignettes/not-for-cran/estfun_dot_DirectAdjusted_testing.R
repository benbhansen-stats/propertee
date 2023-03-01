#' This script serves as validation for the adjustments made to `DirectAdjusted`'s
#' `sandwich::estfun` method in February 2023. The goal of this change was to 
#' write a new `estfun` method such that when `sandwich::vcovCL` was called on a
#' `DirectAdjusted` object (and thus invoking `estfun`) it would return a meat
#' matrix corresponding to the covariance matrix of our stacked estimating
#' equations. Thus, this script tests 1) whether the new method produces the
#' output matrix we expect and 2) whether the resulting sandwich covariance
#' estimate matches what we previously calculated: \eqn{B_{22} - A_{21}A_{11}^{-1}
#' B_{12} - B_{12}^{T}A_{11}^{-1}A_{21}^{T} + A_{21}A_{11}^{-1}B_{11}A_{11}^{-1}
#' A_{21}^{T}, with these matrices defend in the appendix of Carroll et al. (2006).
#' - Josh Wasserman, February 2023

library(sandwich)
library(testthat)
library(flexida)

# DEFINE CLASSES AND FUNCTIONS FOR TESTING
# Create placeholder classes for this script so these tests will hold regardless
# of package development
# `newPsiTildeDAClass` will be use new `estfun.DirectAdjusted` method, the
# `oldPsiTildeDAClass` will use the base class's `DirectAdjusted` method, and
# the rest of the functions we re-write here because they will no longer be
# exported as part of the `flexida` package
setClass("newPsiTildeDAClass",
         contains = "DirectAdjusted")

estfun.newPsiTildeDAClass <- function(object) {
  ## this vector indicates the hierarchy of `sandwich::estfun` methods to use
  ## to extract estimating equations for ITT model
  valid_classes <- c("glm", "lmrob", "svyglm", "lm")
  base_class <- match(object@.S3Class, valid_classes)
  if (all(is.na(base_class))) {
    stop(paste("ITT effect model must have been fitted using a function from the",
               "`flexida`, `stats`, `robustbase`, or `survey` package"))
  }
  psi <- getS3method("estfun", valid_classes[min(base_class, na.rm = TRUE)])(object)
  
  ## if ITT model offset doesn't contain info about covariance model, psi should
  ## be the matrix of estimating equations returned
  ca <- object$model$`(offset)`
  if (is.null(ca) | !inherits(ca, "SandwichLayer")) {
    return(psi)
  }
  
  ## otherwise, extract/compute the rest of the relevant matrices/quantities
  cmod <- ca@fitted_covariance_model
  C_uoas <- ca@keys
  phi <- estfun(cmod)
  nc <- nrow(phi)
  uoa_cols <- colnames(C_uoas)
  
  ## figure out if rows need to be added to the matrix of estimating equations
  Q_uoas <- stats::expand.model.frame(object, uoa_cols, na.expand = TRUE)[, uoa_cols, drop = FALSE]
  Q_uoas <- apply(Q_uoas, 1, function(...) paste(..., collapse = "_"))
  
  C_uoas <- apply(C_uoas, 1, function(...) paste(..., collapse = "_"))
  C_uoas[vapply(strsplit(C_uoas, "_"), function(x) all(x == "NA"), logical(1))] <- NA_character_
  
  nas <- is.na(C_uoas)
  if (any(nas)) {
    # give unique ID's to uoa's in C but not Q
    n_Q_uoas <- length(unique(Q_uoas))
    C_uoas[nas] <- paste0(n_Q_uoas + seq_len(sum(nas)), "*")
  }
  
  # add rows if necessary
  add_C_uoas <- setdiff(unique(C_uoas), unique(Q_uoas))
  add_Q_uoas <- setdiff(unique(Q_uoas), unique(C_uoas))
  if (length(add_C_uoas) > 0) {
    psi <- rbind(psi, matrix(0, nrow = sum(C_uoas %in% add_C_uoas), ncol = ncol(psi)))
  }
  if (length(add_Q_uoas) > 0) {
    phi <- rbind(matrix(0, nrow = sum(Q_uoas %in% add_Q_uoas), ncol = ncol(phi)), phi)
  }
  
  ## form matrix of estimating equations
  n <- nrow(psi)
  a11_inv <- .get_a11_inverse(object)
  a21 <- .get_a21(object)
  mat <- psi - phi %*% a11_inv %*% t(a21)
  
  return(mat)
}

setClass("oldPsiTildeDAClass",
         contains = "DirectAdjusted")

estfun.oldPsiTildeDAClass <- function(object) {
  valid_classes <- c("glm", "lmrob", "svyglm", "lm")
  base_class <- match(object@.S3Class, valid_classes)
  if (all(is.na(base_class))) {
    stop(paste("ITT effect model must have been fitted using a function from the",
               "`flexida`, `stats`, `robustbase`, or `survey` package"))
  }
  psi <- getS3method("estfun", valid_classes[min(base_class, na.rm = TRUE)])(object)
  
  return(psi)
}

#' The following submatrix computation functions are as they appeared in the
#' `flexida` package prior to Feb. 2023
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

.get_a21 <- function(x) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }
  
  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }
  
  # Get contribution to the estimating equation from the ITT effect model
  w <- if (is.null(x$weights)) 1 else x$weights
  
  nq <- nrow(x$model)
  extra_nc <- sum(apply(is.na(sl@keys), 1, any))
  n <- nq + extra_nc
  
  damod_mm <- stats::model.matrix(formula(x),
                                  stats::model.frame(x, na.action = na.pass))
  msk <- (apply(!is.na(sl@prediction_gradient), 1, all) &
            apply(!is.na(damod_mm), 1, all))
  
  out <- crossprod(damod_mm[msk, , drop = FALSE] * w,
                   sl@prediction_gradient[msk, , drop = FALSE])
  
  return(out)
}

.get_a22_inverse <- function(x) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }
  
  # Get expected information per sandwich_infrastructure vignette
  w <- if (is.null(x$weights)) 1 else x$weights
  out <- solve(crossprod(stats::model.matrix(x) * sqrt(w)))
  
  return(out)
}

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
  uoas <- sl@keys
  cluster_cols <- colnames(uoas)
  
  # Replace NA's for rows not in the experimental design with a unique cluster ID
  if (ncol(uoas) == 1) {
    uoas <- uoas[, 1]
  } else {
    uoas <- Reduce(function(...) paste(..., sep = "_"), uoas[, cluster_cols])
    uoas[vapply(strsplit(uoas, "_"), function(x) all(x == "NA"), logical(1))] <- NA_character_
  }
  
  nuoas <- length(unique(uoas))
  nas <- is.na(uoas)
  if (any(nas)) {
    uoas[nas] <- paste0(nuoas - 1 + seq_len(sum(nas)), "*")
  }
  uoas <- factor(uoas)
  nuoas <- length(unique(uoas))
  
  # if the covariance model only uses one cluster, produce warning
  if (nuoas == 1) {
    warning(paste("Covariance adjustment model has meat matrix numerically",
                  "indistinguishable from 0"))
    return(
      matrix(0, nrow = ncol(stats::model.matrix(cmod)), ncol = ncol(stats::model.matrix(cmod)))
    )
  }
  dots$cluster <- uoas
  dots$x <- cmod
  
  out <- do.call(sandwich::meatCL, dots) * nc
  return(out)
}

.get_b12 <- function(x) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }
  
  sl <- x$model$`(offset)`
  if (!is(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }
  
  cluster_cols <- var_names(x@Design, "u")
  keys <- sl@keys
  
  if (ncol(keys) == 1) {
    uoas <- keys[, 1]
  } else {
    uoas <- Reduce(function(...) paste(..., sep = "_"), keys)
    uoas[grepl("NA", uoas)] <- NA_character_
  }
  
  message(paste(sum(!is.na(uoas)),
                "rows in the covariance adjustment model",
                "data joined to the ITT effect model data\n"))
  
  
  # Check number of overlapping clusters n_QC; if n_QC <= 1, return 0 (but
  # similarly to .get_b11(), throw a warning when only one cluster overlaps)
  uoas_overlap <- length(unique(uoas))
  if (uoas_overlap == 1) {
    if (!is.na(unique(uoas))) {
      warning(paste("Covariance matrix between covariance adjustment and ITT effect",
                    "model estimating equations is numerically indistinguishable",
                    "from 0"))
    }
    return(
      matrix(0,
             nrow = dim(stats::model.matrix(sl@fitted_covariance_model))[2],
             ncol = dim(stats::model.matrix(x))[2])
    )
  }
  
  # Sum est eqns to cluster level; since non-overlapping rows are NA in `keys`,
  # `by` call excludes them from being summed
  cmod_estfun <- sandwich::estfun(sl@fitted_covariance_model)
  cmod_aggfun <- ifelse(dim(cmod_estfun)[2] > 1, colSums, sum)
  cmod_eqns <- Reduce(rbind, by(cmod_estfun, uoas, cmod_aggfun))
  
  # get rows from overlapping clusters in experimental data
  Q_uoas <- stats::expand.model.frame(x, cluster_cols, na.expand = TRUE)[cluster_cols]
  if (ncol(Q_uoas) == 1) {
    Q_uoas <- Q_uoas[, 1]
  } else {
    Q_uoas <- Reduce(function(...) paste(..., sep = "_"), Q_uoas)
    Q_uoas[grepl("NA", Q_uoas)] <- NA_integer_
  }
  
  msk <- Q_uoas %in% unique(uoas[!is.na(uoas)])
  damod_estfun <- sandwich::estfun(x)[msk, , drop = FALSE]
  damod_aggfun <- ifelse(dim(damod_estfun)[2] > 1, colSums, sum)
  damod_eqns <- Reduce(rbind, by(damod_estfun, Q_uoas[msk], damod_aggfun))
  
  matmul_func <- if (length(unique(uoas[!is.na(uoas)])) == 1) tcrossprod else crossprod
  return(matmul_func(cmod_eqns, damod_eqns))
}

.get_b22 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }
  
  nq <- nrow(sandwich::estfun(x))
  
  # Create cluster ID matrix depending on cluster argument (or its absence)
  dots <- list(...)
  uoas <- stats::expand.model.frame(x,
                                    var_names(x@Design, "u"))[, var_names(x@Design, "u"),
                                                              drop = FALSE]
  
  if (ncol(uoas) == 1) {
    uoas <- factor(uoas[,1])
  } else {
    uoas <- factor(Reduce(function(...) paste(..., sep = "_"), uoas))
  }
  dots$cluster <- uoas
  dots$x <- x
  
  out <- do.call(sandwich::meatCL, dots) * nq
  
  return(out)
}

#' This function computes the new meat matrix derived in the notes. We compare
#' the crossproduct of the unit-of-assignment level estimating equations (cluster
#' the output of `estfun` if need be) to the output of this function.
compute_meat_matrix_no_estfun <- function(object) {
  cmod <- object$model$`(offset)`@fitted_covariance_model
  phi <- estfun(cmod)
  psi <- estfun(object)
  keys <- object$model$`(offset)`@keys
  nc <- nrow(keys)

  Q_uoas <- stats::expand.model.frame(object, colnames(keys))[, colnames(keys), drop = FALSE]
  if (ncol(Q_uoas) == 1) {
    Q_uoas <- Q_uoas[, 1]
  } else {
    Q_uoas <- apply(Q_uoas, 1, function(...) paste(..., collapse = "_"))
  }
  n_Q_uoas <- length(unique(Q_uoas))

  C_uoas <- apply(keys, 1, function(...) paste(..., collapse = "_"))
  C_uoas[vapply(strsplit(C_uoas, "_"), function(x) all(x == "NA"), logical(1))] <- NA_character_
  nas <- is.na(C_uoas)
  if (any(nas)) {
    C_uoas[nas] <- paste0(n_Q_uoas + seq_len(sum(nas)), "*")
  }

  add_C_uoas <- setdiff(unique(C_uoas), unique(Q_uoas))
  add_Q_uoas <- setdiff(unique(Q_uoas), unique(C_uoas))
  if (length(add_C_uoas) > 0) {
    psi <- rbind(psi, matrix(0, nrow = sum(C_uoas %in% add_C_uoas), ncol = ncol(psi)))
  }
  if (length(add_Q_uoas) > 0) {
    phi <- rbind(matrix(0, nrow = sum(Q_uoas %in% add_Q_uoas), ncol = ncol(phi)), phi)
  }
  psi_aggfun <- ifelse(dim(psi)[2] > 1, colSums, sum)
  phi_aggfun <- ifelse(dim(phi)[2] > 1, colSums, sum)
  psi_uoa_order <- factor(c(Q_uoas, C_uoas[C_uoas %in% add_C_uoas]),
                          levels = unique(c(Q_uoas, C_uoas[C_uoas %in% add_C_uoas])))
  phi_uoa_order <- factor(c(C_uoas, Q_uoas[Q_uoas %in% add_Q_uoas]),
                          levels = unique(c(C_uoas, Q_uoas[Q_uoas %in% add_Q_uoas])))
  psi <- Reduce(rbind, by(psi, psi_uoa_order, psi_aggfun))
  phi <- Reduce(rbind, by(phi, phi_uoa_order, phi_aggfun))

  n <- nrow(psi)
  a11_inv <- .get_a11_inverse(object) * nc # sandwich::bread(object$model$`(offset)`@fitted_covariance_model)
  a21 <- .get_a21(object) / n
  m11 <- crossprod(phi) / nc
  m12 <- crossprod(phi, psi) / sqrt(n * nc)
  m22 <- crossprod(psi) / n
  
  return(
    n * (m22 - sqrt(n / nc) * (t(m12) %*% a11_inv %*% t(a21) + a21 %*% a11_inv %*% m12) +
      n / nc * (a21 %*% a11_inv %*% m11 %*% a11_inv %*% t(a21)))
  )
}

#' This function computes the full model-based sandwich variance estimate using
#' the new `estfun` method
new_vcovMB_CR0 <- function(object, ...) {
  cluster <- var_names(object@Design, "u")
  # Get the unit of assignment ID's in Q given the manual cluster argument or the Design info
  Q_uoas <- stats::expand.model.frame(object, cluster, na.expand = TRUE)[, cluster, drop = FALSE]
  Q_uoas <- apply(Q_uoas, 1, function(...) paste(..., collapse = "_"))
  names(Q_uoas) <- NULL

  # Get the unit of assignment ID's in C
  ca <- object$model$`(offset)`
  C_uoas <- ca@keys
  C_uoas <- apply(C_uoas, 1, function(...) paste(..., collapse = "_"))
  names(C_uoas) <- NULL
  C_uoas[vapply(strsplit(C_uoas, "_"), function(x) all(x == "NA"), logical(1))] <- NA_character_

  all_uoas <- Q_uoas
  nas <- is.na(C_uoas)
  if (any(nas)) {
    # give unique ID's to units of assignment in C but not Q, and concatenate with ID's in Q
    n_Q_uoas <- length(unique(Q_uoas))
    new_C_uoas <- paste0(n_Q_uoas + seq_len(sum(nas)), "*")
    all_uoas <- c(Q_uoas, new_C_uoas)
  }
  
  uoa_ids <- factor(all_uoas, levels = unique(all_uoas))
  n <- length(uoa_ids)

  a22inv <- n * .get_a22_inverse(object) # `sandwich::bread()` but with correct scaling
  meat <- crossprod(Reduce(rbind, by(estfun(object), uoa_ids, colSums))) / n # same as `sandwich::meatCL()`
  vmat <- (1 / n) * a22inv %*% meat %*%  a22inv

  return(vmat)
}

#' This function computes the model-based sandwich variance estimate using
#' the submatrix computation functions that were used in the `flexida` package
#' prior to Feb. 2023.
old_vcovMB_CR0 <- function(x, ...) {
  # compute blocks
  a21 <- .get_a21(x)
  a11inv <- .get_a11_inverse(x)
  b12 <- .get_b12(x)
  
  a22inv <- .get_a22_inverse(x)
  b22 <- .get_b22(x, type = "HC0", ...)
  b11 <- .get_b11(x,  type = "HC0", ...)
  
  meat <- (
    b22 -
      a21 %*% a11inv %*% b12 -
      t(b12) %*% t(a11inv) %*% t(a21) +
      a21 %*% a11inv %*% b11 %*% t(a11inv) %*% t(a21)
  )
  vmat <- a22inv %*% meat %*% a22inv
  
  return(vmat)
}

# FORMULATE 4 DIRECT ADJUSTMENT MODELS TO TEST:
# 1) linear ITT effect model when there's overlap between C and Q + no clustering
# 2) linear ITT effect model when there's no overlap + no clustering
# 3) linear ITT effect model with overlap + clustering
# 4) nonlinear covariance adjustment model
# For each model, verify:
#  1) crossprod of new est. eqns. returns what we expect
#  2) new vcovDA estimates are the same estimates output by the `flexida`
#     package prior to Feb. 2023
#  3) new vcovDA estimates are the same as `sandwich::sandwich(newDA, meat. = sandwich::vcovCL)`
#     with appropriate scaling for the `bread` function and an appropriate
#     `cluster` argument
set.seed(2300)

# TEST MODEL 1
n <- 1000
newdata <- cbind(matrix(rnorm(2 * n), ncol = 2), rbinom(n, 1, 0.5))
beta <- c(1, -0.5, 0.25)
newdata <- as.data.frame(cbind(seq_len(n), newdata, newdata %*% beta + 0.1 * rnorm(n)))
colnames(newdata) <- c("uoa_id", "x1", "x2", "z", "y")

cmod1 <- lm(y ~ x1 + x2, newdata)
des1 <- rct_design(z ~ unitid(uoa_id), data = newdata)
damod1 <- lmitt(y ~ assigned(), data = newdata, design = des1, weights = ate(des1),
                offset = cov_adj(cmod1))
new_psi_tilde_da1 <- new("newPsiTildeDAClass", damod1)
old_psi_tilde_da1 <- new("oldPsiTildeDAClass", damod1)

testthat::expect_equal(crossprod(estfun(new_psi_tilde_da1)),
                       compute_meat_matrix_no_estfun(old_psi_tilde_da1))
testthat::expect_equal(new_vcov1 <- new_vcovMB_CR0(new_psi_tilde_da1),
                       old_vcovMB_CR0(old_psi_tilde_da1, cadjust = FALSE))
testthat::expect_equal(new_vcov1,
                       sandwich::sandwich(new_psi_tilde_da1, adjust = FALSE))
new_vcov1

# TEST MODEL 2
cmod_data <- cbind(matrix(rnorm(4 * n), ncol = 2), rep(0, n))
cmod_data <- as.data.frame(cbind(seq(n + 1, 3 * n), cmod_data, cmod_data %*% beta + 0.1 * rnorm(2 * n)))
colnames(cmod_data) <- c("uoa_id", "x1", "x2", "z", "y")

cmod2 <- lm(y ~ x1 + x2, cmod_data)
des2 <- rct_design(z ~ unitid(uoa_id), data = newdata)
damod2 <- lmitt(y ~ assigned(), data = newdata, design = des2, weights = ate(des2),
                offset = cov_adj(cmod2))
new_psi_tilde_da2 <- new("newPsiTildeDAClass", damod2)
old_psi_tilde_da2 <- new("oldPsiTildeDAClass", damod2)

testthat::expect_equal(crossprod(estfun(new_psi_tilde_da2)),
                       compute_meat_matrix_no_estfun(old_psi_tilde_da2))
testthat::expect_equal(new_vcov2 <- new_vcovMB_CR0(new_psi_tilde_da2),
                       old_vcovMB_CR0(old_psi_tilde_da2, cadjust = FALSE))
testthat::expect_equal(new_vcov2,
                       sandwich::sandwich(new_psi_tilde_da2,
                                          bread. = function(x) {
                                            sandwich::bread(x) / n * 3 * n
                                          },
                                          meat. = sandwich::meatCL,
                                          cluster = seq_len(3 * n),
                                          cadjust = FALSE))
new_vcov2

# TEST MODEL 3
n_clusters <- 10
newdata <- cbind(matrix(rnorm(2 * n), ncol = 2),
                 rep(rbinom(n_clusters, 1, 0.5), each = as.integer(n / n_clusters)))
newdata <- as.data.frame(
  cbind(rep(seq_len(10), each = n / 10), newdata, newdata %*% beta + 0.1 * rnorm(n))
)
colnames(newdata) <- c("uoa_id", "x1", "x2", "z", "y")

cmod3 <- lm(y ~ x1 + x2, newdata)
des3 <- rct_design(z ~ cluster(uoa_id), data = newdata)
damod3 <- lmitt(y ~ assigned(), data = newdata, design = des3, weights = ate(des3),
                offset = cov_adj(cmod3))
new_psi_tilde_da3 <- new("newPsiTildeDAClass", damod3)
old_psi_tilde_da3 <- new("oldPsiTildeDAClass", damod3)

testthat::expect_equal(
  crossprod(Reduce(rbind, by(estfun(new_psi_tilde_da3), newdata$uoa_id, colSums))),
  compute_meat_matrix_no_estfun(old_psi_tilde_da3)
)
testthat::expect_equal(new_vcov3 <- new_vcovMB_CR0(new_psi_tilde_da3),
                       old_vcovMB_CR0(old_psi_tilde_da3, cadjust = FALSE))
testthat::expect_equal(
  new_vcov3,
  sandwich::sandwich(new_psi_tilde_da3,
                     meat. = sandwich::meatCL,
                     cluster = newdata$uoa_id,
                     cadjust = FALSE)
)
new_vcov3

# TEST MODEL 4
newdata <- cbind(matrix(rnorm(2 * n), ncol = 2), rbinom(n, 1, 0.5))
beta <- c(1, -0.5, 0.25)
newdata <- as.data.frame(cbind(seq_len(n), newdata,
                               rbinom(n = n, size = 1, p = 1 / (1 + exp(-(newdata %*% beta + 0.1 * rnorm(n)))))))
colnames(newdata) <- c("uoa_id", "x1", "x2", "z", "y")

cmod4 <- glm(y ~ x1 + x2, newdata, family = binomial())
des4 <- rct_design(z ~ unitid(uoa_id), data = newdata)
damod4 <- lmitt(y ~ assigned(), data = newdata, design = des4, weights = ate(des4),
                offset = cov_adj(cmod4))
new_psi_tilde_da4 <- new("newPsiTildeDAClass", damod4)
old_psi_tilde_da4 <- new("oldPsiTildeDAClass", damod4)

testthat::expect_equal(crossprod(estfun(new_psi_tilde_da4)),
                       compute_meat_matrix_no_estfun(old_psi_tilde_da4))
testthat::expect_equal(new_vcov4 <- new_vcovMB_CR0(new_psi_tilde_da4),
                       old_vcovMB_CR0(old_psi_tilde_da4, cadjust = FALSE))
testthat::expect_equal(new_vcov4,
                       sandwich::sandwich(new_psi_tilde_da4, adjust = FALSE))
new_vcov4

cat("Successfully tested estfun.DirectAdjusted", sep = "\n")
