#' @include Design.R SandwichLayer.R
NULL

#' @title Compute covariance-adjusted cluster-robust sandwich variance estimates
#' @param x A fitted \code{DirectAdjusted} model object
#' @param type A string indicating the desired variance estimator. Currently
#' accepts "CR0"
#' @param ... Arguments to be passed to the internal variance estimation function.
#' One argument a user may want to manually override is the `cluster` argument.
#' Users may be interested in clustering standard errors at levels different than
#' the unit of assignment level specified in the \code{DirectAdjusted} model's
#' \code{Design}. In this case, they may specify column names in the ITT effect/covariance
#' adjustment model datasets corresponding to a different clustering level.
#' \code{.get_b11()} and \code{.get_b12()} tolerate NA values for manually
#' provided cluster ID's, since the covariance adjustment model may be fit to a
#' sample larger than the quasiexperimental sample. In this case, each unit with
#' an NA cluster ID will be treated as an independent observation. \code{.get_b22()},
#' however, will fail if the new cluster levels contain NA's because all units
#' in the quasiexperimental sample should have well-defined design information.
#' @return A \eqn{2\times 2} matrix where the dimensions are
#' given by the intercept and treatment variable terms in the ITT effect model
#' @export
#' @rdname var_estimators
vcovDA <- function(x, type = c("CR0"), ...) {
  type <- match.arg(type)

  var_func <- switch(
    type,
    "CR0" = .vcovMB_CR0
  )
  args <- list(...)
  args$x <- x
  args$cluster <- .make_uoa_ids(x, ...)
  
  est <- do.call(var_func, args)
  return(est)
}

#' Make unit of assignment ID's that align with the output of `estfun.DirectAdjusted`
#' @keywords internal
.make_uoa_ids <- function(x, cluster = NULL, ...) {
  if (!inherits(cluster, "character") & !inherits(x, "DirectAdjusted")) {
    stop(paste("Cannot deduce units of assignment for clustering without a",
               "Design object (stored in a DirectAdjusted object) or a `cluster`",
               "argument specifying the columns with the units of assignment"))
  }
  
  # Must be a DirectAdjusted object for this logic to occur
  if (!inherits(cluster, "character")) {
    cluster <- var_names(x@Design, "u")
  }
  
  # Get the unit of assignment ID's in Q given the manual cluster argument or the Design info
  Q_uoas <- tryCatch(
    stats::expand.model.frame(x, cluster, na.expand = TRUE)[, cluster, drop = FALSE],
    error = function(e) {
      mf <- eval(x$call$data, envir = environment(x))
      missing_cols <- setdiff(cluster, colnames(mf))
      stop(paste("Could not find unit of assignment columns",
                 paste(missing_cols, collapse = ", "), "in ITT effect model data"),
           call. = FALSE)
    })
  Q_uoas <- apply(Q_uoas, 1, function(...) paste(..., collapse = "_"))
  names(Q_uoas) <- NULL

  # If there's no covariance adjustment info, return the ID's found in Q
  ca <- x$model$`(offset)`
  if (!inherits(x, "DirectAdjusted") | !inherits(ca, "SandwichLayer")) {
    return(factor(Q_uoas, levels = unique(Q_uoas)))
  }

  # Get the unit of assignment ID's in C
  all_uoas <- Q_uoas
  C_uoas <- ca@keys
  any_C_uoa_is_na <- apply(is.na(C_uoas), 1, any)
  all_C_uoa_is_na <- apply(is.na(C_uoas), 1, all)
  C_uoas <- apply(C_uoas, 1, function(...) paste(..., collapse = "_"))
  
  if (sum(any_C_uoa_is_na) - sum(all_C_uoa_is_na) > 0) {
    warning(paste("Some rows in the covariance adjustment model dataset have",
                  "NA's for some but not all clustering columns. Rows sharing",
                  "the same non-NA cluster ID's will be clustered together.",
                  "If this is not intended, provide unique non-NA cluster ID's",
                  "for these rows."))
  }
  if (sum(all_C_uoa_is_na) > 0) {
    warning(paste("Some or all rows in the covariance adjustment model dataset",
                  "are found to have NA's for the given clustering columns.",
                  "This is taken to mean these observations should be treated",
                  "as IID. To avoid this warning, provide unique non-NA cluster",
                  "ID's for each row."))
    # give unique ID's to units of assignment in C but not Q, and concatenate with ID's in Q
    C_uoas[all_C_uoa_is_na] <- NA_character_
    n_Q_uoas <- length(unique(Q_uoas))
    C_uoas[all_C_uoa_is_na] <- paste0(n_Q_uoas + seq_len(sum(all_C_uoa_is_na)), "*")
    all_uoas <- c(Q_uoas, C_uoas)
  }
  
  return(factor(all_uoas, levels = unique(all_uoas)))
}

#' Model-based standard errors with HC0 adjustment
#' @keywords internal
#' @rdname var_estimators
.vcovMB_CR0 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  args <- list(...)
  if ("type" %in% names(args)) {
    stop(paste("Cannot override the `type` argument for meat",
               "matrix computations"))
  }
  args$x <- x
  n <- length(args$cluster)
  
  a22inv <- sandwich::bread(x)
  meat <- do.call(sandwich::meatCL, args)
  vmat <- (1 / n) * a22inv %*% meat %*% a22inv

  return(vmat)
}

#' @title (Internal) Compute variance blocks
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
#' ITT effect model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a22_inverse <- function(x) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  # Get expected information per sandwich_infrastructure vignette
  w <- if (is.null(x$weights)) 1 else x$weights
  mm <- stats::model.matrix(x)
  out <- solve(crossprod(mm * sqrt(w)))

  return(out)
}

#' @details The \bold{A11 block} is the \eqn{p\times p} matrix corresponding to
#'   the unscaled inverse of the observed Fisher information of the covariance
#'   adjustment model. The observed information is given by the estimate of the negative
#'   Jacobian of the model's estimating equations. The unscaled version provided
#'   here divides by the number of observations used to fit the covariance
#'   adjustment model.
#' @return \code{.get_a11_inverse()}: A \eqn{p\times p} matrix where the
#'   dimensions are given by the number of terms in the covariance adjustment model
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

#' @details The \bold{A21 block} is the block of the sandwich variance estimator
#'   corresponding to the gradient of the ITT effect model with respect
#'   to the covariates. Some of the information needed for this calculation is
#'   stored in the \code{DirectAdjusted} object's \code{SandwichLayer} offset. This
#'   block is the crossproduct of the prediction gradient and the gradient of
#'   the conditional mean vector for the ITT effect model summed to the
#'   cluster level. In other words, we take this matrix to be \deqn{\sum(d\psi_i
#'   / d\alpha) = -\sum(w_i/\phi) * (d\mu(\eta_i) / d\eta_i) *
#'   (d\upsilon(\zeta_i) / d\zeta_i) * (x_i c_i)x_i'} where \eqn{\mu} and
#'   \eqn{\eta_i} are the conditional mean function and linear predictor for the
#'   ith cluster in the ITT effect model, and \eqn{\upsilon} and
#'   \eqn{\zeta_i} are the conditional mean function and linear predictor for
#'   the ith cluster in the covariance adjustment model.
#' @return \code{.get_a12()}: A \eqn{2\times p} matrix where the number of
#'   rows are given by intercept and treatment variable terms in the ITT effect
#'   model, and the number of columns are given by the number of terms
#'   in the covariance adjustment model
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

  # Get contribution to the estimating equation from the ITT effect model
  w <- if (is.null(x$weights)) 1 else x$weights

  damod_mm <- stats::model.matrix(formula(x),
                                  stats::model.frame(x, na.action = na.pass))
  msk <- (apply(!is.na(sl@prediction_gradient), 1, all) &
            apply(!is.na(damod_mm), 1, all))

  out <- crossprod(damod_mm[msk, , drop = FALSE] * w,
                   sl@prediction_gradient[msk, , drop = FALSE])

  return(out)
}
