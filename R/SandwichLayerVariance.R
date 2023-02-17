#' @include Design.R SandwichLayer.R
NULL

#' @title Compute covariance-adjusted cluster-robust sandwich variance estimates
#' @param object A \code{DirectAdjusted} model
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
vcovDA <- function(object, type = c("CR0"), ...) {
  type <- match.arg(type)

  var_func <- switch(
    type,
    "CR0" = .vcovMB_CR0
  )

  est <- var_func(object, ...)
  return(est)
}

#' @keywords internal
#' @rdname var_estimators
.vcovMB_CR0 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  m <- match.call()
  if ("type" %in% names(m)) {
    stop(paste("Cannot override the `type` argument for meat",
               "matrix computations"))
  }

  # compute blocks
  a22inv <- .get_a22_inverse(x)
  b22 <- .get_b22(x, type = "HC0", ...)
  
  if (!inherits(x$model$`(offset)`, "SandwichLayer")) {
    meat <- b22
  } else {
    a21 <- .get_a21(x)
    a11inv <- .get_a11_inverse(x)
    b12 <- .get_b12(x, ...)
    b11 <- .get_b11(x,  type = "HC0", ...)
    
    meat <- (
      b22 -
        a21 %*% a11inv %*% b12 -
        t(b12) %*% t(a11inv) %*% t(a21) +
        a21 %*% a11inv %*% b11 %*% t(a11inv) %*% t(a21)
    )
  }

  vmat <- a22inv %*% meat %*% a22inv

  return(vmat)
}

#' @title (Internal) Compute variance blocks
#' @details The \bold{B12 block} is the covariance matrix of the cluster-level
#'   estimating equations for the covariance adjustment and ITT effect models. It
#'   has a row for each term in the covariance adjustment model and a column for each term
#'   in the ITT effect model. For any row that does not appear in both
#'   the experimental design and the covariance adjustment model data, its contribution to
#'   this matrix will be 0. Thus, if there is no overlap between the two
#'   datasets, this will return a matrix of 0's.
#' @param x A \code{DirectAdjusted} model
#' @return \code{.get_b12()}: A \eqn{p\times 2} matrix where the number of rows
#'   are given by the number of terms in the covariance adjustment model and the number of
#'   columns correspond to intercept and treatment variable terms in the ITT
#'   effect model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_b12 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  sl <- x$model$`(offset)`
  if (!is(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  # If cluster argument is NULL, use the `keys` dataframe created at initialization,
  # otherwise recreate given the desired clustering columns
  dots <- list(...)
  if (is.null(dots$cluster)) {
    cluster_cols <- var_names(x@Design, "u")
    keys <- sl@keys
  } else if (inherits(dots$cluster, "character")) {
    cluster_cols <- dots$cluster
    wide_frame <- tryCatch(
      stats::expand.model.frame(sl@fitted_covariance_model, cluster_cols,
                                na.expand = TRUE)[cluster_cols],
      error = function(e) {
        data <- eval(sl@fitted_covariance_model$call$data,
                     envir = environment(formula(sl@fitted_covariance_model)))
        stop(paste("The columns",
                   paste(setdiff(cluster_cols, colnames(data)), collapse = ", "),
                   "are missing from the covariance adjustment model dataset"),
             call. = FALSE)
      })
    
    # Check to see if provided column names overlap with design
    # missing_des_cols <- setdiff(cov_adj_cluster_cols, colnames(x@Design@structure))
    missing_des_cols <- setdiff(cluster_cols, colnames(x@Design@structure))
    if (length(missing_des_cols) > 0) {
      stop(paste("The following columns in the `cluster` argument cannot be found",
                 "in the DirectAdjusted object's Design:",
                 paste(missing_des_cols, collapse = ", ")))
    }
    
    # Re-create keys dataframe with the new clustering columns
    keys <- as.data.frame(
      sapply(cluster_cols, function(col) {
        match(wide_frame[[col]], unique(x@Design@structure[[col]]), incomparables = NA)
      })
    )
  } else {
    stop(paste("If overriding `cluster` argument for meat matrix calculations,",
               "must provide a character vector specifying column names that",
               "exist in both the ITT effect and covariance model datasets"))
  }

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
  out <- solve(crossprod(stats::model.matrix(x) * sqrt(w)))

  return(out)
}

#' @param ... Arguments to be passed to sandwich::meatCL
#' @details The \bold{B22 block} refers to a clustered estimate of the covariance
#'   matrix for the ITT effect model. The \code{stats} package offers family objects
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
#'   ITT effect model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_b22 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  nq <- nrow(sandwich::estfun(x))

  # Create cluster ID matrix depending on cluster argument (or its absence)
  dots <- list(...)
  if (is.null(dots$cluster)) {
    uoas <- stats::expand.model.frame(x,
                                      var_names(x@Design, "u"))[, var_names(x@Design, "u"),
                                                                drop = FALSE]
  } else if (inherits(dots$cluster, "character")) {
    uoas <- tryCatch(
      stats::expand.model.frame(x, dots$cluster)[, dots$cluster, drop = FALSE],
      error = function(e) {
        data <- eval(x$call$data,
                     envir = environment(formula(x)))
        stop(paste("The columns",
                   paste(setdiff(dots$cluster, colnames(data)), collapse = ", "),
                   "are missing from the ITT effect model dataset"),
             call. = FALSE)
      })
  } else {
    stop(paste("If overriding `cluster` argument for meat matrix calculations,",
               "must provide a character vector specifying column names in the",
               "ITT effect model dataset"))
  }

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

#' @details The \bold{B11 block} is the block of the sandwich variance estimator
#'   corresponding to the variance-covariance matrix of the covariance model
#'   coefficient estimates. The estimates returned here are potentially
#'   clustered (either by the clustering in the experimental design or by
#'   a manually provided `cluster` argument) if clustering information can be
#'   retrieved from the covariance adjustment model data.
#' @return \code{.get_b11()}: A \eqn{p\times p} matrix where the dimensions are
#'   given by the number of terms in the covariance adjustment model including an
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
    uoas <- sl@keys
    cluster_cols <- colnames(uoas)
  } else if (inherits(dots$cluster, "character")) {
    cluster_cols <- dots$cluster
    uoas <- tryCatch(
      stats::expand.model.frame(cmod, cluster_cols)[, cluster_cols, drop = FALSE],
      error = function(e) {
        data <- eval(cmod$call$data,
                     envir = environment(formula(cmod)))
        stop(paste("The columns",
                   paste(setdiff(cluster_cols, colnames(data)), collapse = ", "),
                   "are missing from the covariance adjustment model dataset"),
             call. = FALSE)
      })

    # check for NA's in the clustering columns
    nas <- rowSums(is.na(uoas[, cluster_cols, drop = FALSE]))
    if (any(nas == length(cluster_cols))) {
      warning(paste("Some or all rows in the covariance adjustment model dataset",
                    "are found to have NA's for the given clustering columns.",
                    "This is taken to mean these observations should be treated",
                    "as IID. To avoid this warning, provide unique non-NA cluster",
                    "ID's for each row."))
    } else if (any(nas > 0 & nas < length(cluster_cols))) {
      warning(paste("Some rows in the covariance adjustment model dataset have",
                    "NA's for some but not all clustering columns. Rows sharing",
                    "the same non-NA cluster ID's will be clustered together.",
                    "If this is not intended, provide unique non-NA cluster ID's",
                    "for these rows."))
    }
  } else {
    stop(paste("If overriding `cluster` argument for meat matrix calculations,",
               "must provide a character vector specifying column names in the",
               "covariance adjustment model dataset"))
  }

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
