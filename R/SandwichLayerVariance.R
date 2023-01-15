#' @include Design.R SandwichLayer.R
NULL

#' @title Compute directly adjusted cluster-robust sandwich variance estimates
#' @param object A \code{DirectAdjusted} model
#' @param type A string indicating the desired variance estimator. Currently
#' accepts "CR0"
#' @param ... Arguments to be passed to the internal variance estimation function.
#' One argument a user may want to manually override is the `cluster` argument.
#' Users may be interested in clustering standard errors at levels different than
#' the unit of assignment level specified in the \code{DirectAdjusted} model's
#' \code{Design}. In this case, they may specify column names in the ITT/covariance
#' adjustment model datasets corresponding to a different clustering level, or
#' dataframes, lists, matrices, or vectors indicating the clusters units pertain
#' to.\n\n\code{.get_b11()} and \code{.get_b12()} tolerate NA values for manually
#' provided cluster ID's, since the covariance adjustment model may be fit to a
#' sample larger than the quasiexperimental sample. In this case, each unit with
#' an NA cluster ID will be treated as an independent observation. \code{.get_b22()},
#' however, will fail if the new cluster levels contain NA's because all units
#' in the quasiexperimental sample should have well-defined design information.
#' 
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

#' @title Compute directly-adjusted CR0 cluster-robust sandwich variance estimates
#' under model-based assumptions
#' @param x A \code{DirectAdjusted} model
#' @param ... Arguments to be passed to sandwich::meatCL
#' @return \code{.vcovMB_CR0()}: A \eqn{2\times 2} matrix where the dimensions are
#' given by the intercept and treatment variable terms in the direct adjustment
#' model
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
  
  m <- match.call()
  if ("type" %in% names(m)) {
    stop(paste("Cannot override the `type` argument for meat",
               "matrix computations"))
  }

  # compute blocks
  a21 <- .get_a21(x)
  a11inv <- .get_a11_inverse(x)
  b12 <- .get_b12(x, ...)

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
.get_b12 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  sl <- x$model$`(offset)`
  if (!is(sl, "SandwichLayer")) {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  cluster_cols <- var_names(x@Design, "u")
  trt_col <- var_names(x@Design, "t")

  # If cluster argument is NULL, use the `keys` dataframe created at initialization,
  # otherwise recreate given the desired clustering columns
  dots <- list(...)
  if (is.null(dots$cluster)) {
    keys <- sl@keys
    cov_adj_cluster_cols <- itt_cluster_cols <- cluster_cols
  } else {
    if (inherits(dots$cluster, "character")) {
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
      cov_adj_cluster_cols <- .check_cluster_col_nas(cluster_cols, wide_frame)
    } else if (inherits(dots$cluster, "data.frame")) {
      cluster_cols <- colnames(dots$cluster)
      wide_frame <- dots$cluster
    } else if (inherits(dots$cluster, c("matrix", "list"))) {
      wide_frame <- as.data.frame(dots$cluster)
      cluster_cols <- colnames(wide_frame)
    } else if (inherits(dots$cluster, c("integer", "numeric", "factor"))) {
      m <- match.call(definition = sys.function(sys.parent(2L)),
                      call = sys.call(sys.parent(2L)))
      cluster_cols <- deparse(m$cluster)
      wide_frame <- as.data.frame(dots$cluster)
      colnames(wide_frame) <- cluster_cols
    } else {
      stop(paste("If overriding `cluster` argument for meat matrix calculations,",
                 "must provide a data frame, matrix, list, numeric/factor vector, or",
                 "a character vector specifying column names that exist in both the ITT and",
                 "covariance model datasets"))
    }
    
    # Check to see if provided column names overlap with design
    missing_des_cols <- setdiff(cluster_cols, colnames(x@Design@structure))
    if (length(missing_des_cols) > 0) {
      stop(paste("The following columns in the `cluster` argument cannot be found",
                 "in the DirectAdjusted object's Design:",
                 paste(missing_des_cols, collapse = ", ")))
    }
    itt_cluster_cols <- .check_cluster_col_nas(cluster_cols,
                                               x@Design@structure,
                                               "ITT effect")
    
    # Re-create keys dataframe with the new clustering columns
    keys <- .merge_preserve_order(wide_frame[cov_adj_cluster_cols],
                                  unique(x@Design@structure[c(cov_adj_cluster_cols, trt_col)]),
                                  all.x = TRUE,
                                  sort = FALSE)
    keys[is.na(keys[, trt_col]), cov_adj_cluster_cols] <- NA
    keys <- keys[, cov_adj_cluster_cols, drop = FALSE]
  }

  if (ncol(keys) == 1) {
    uoas <- keys[, 1]
  } else {
    uoas <- Reduce(function(...) paste(..., sep = "_"), keys)
    uoas[grepl("NA", uoas)] <- NA_character_
  }

  # Check number of overlapping clusters n_QC; if n_QC == 0, return 0
  no_uoas_overlap <- all(is.na(unique(uoas)))
  if (no_uoas_overlap) {
    return(matrix(0,
                  nrow = dim(stats::model.matrix(sl@fitted_covariance_model))[2],
                  ncol = dim(stats::model.matrix(x))[2]))
  }

  # Sum est eqns to cluster level; since non-overlapping rows are NA in `keys`,
  # `by` call excludes them from being summed
  cmod_estfun <- sandwich::estfun(sl@fitted_covariance_model)
  cmod_aggfun <- ifelse(dim(cmod_estfun)[2] > 1, colSums, sum)
  cmod_eqns <- Reduce(rbind, by(cmod_estfun, uoas, cmod_aggfun))

  # get rows from overlapping clusters in experimental data
  Q_uoas <- stats::expand.model.frame(x, itt_cluster_cols, na.expand = TRUE)[itt_cluster_cols]
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
                   "are missing from the ITT model dataset"),
             call. = FALSE)
      })
  } else if (inherits(dots$cluster, "data.frame")) {
    uoas <- dots$cluster
  } else if (inherits(dots$cluster, c("matrix", "list", "integer", "numeric", "factor"))) {
    uoas <- as.data.frame(dots$cluster)
  } else {
    stop(paste("If overriding `cluster` argument for meat matrix calculations,",
               "must provide a data frame, matrix, list, numeric/factor vector, or",
               "a character vector specifying column names in the ITT model dataset"))
  }

  # NOTE: if user passes in matrix with multiple columns, they are concatenated
  # and forced to be one factor column (as we do in other situations, since we
  # assume nested clusters only)
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
    uoas <- sl@keys
  } else if (inherits(dots$cluster, "character")) {
    uoas <- tryCatch(
      stats::expand.model.frame(cmod, dots$cluster)[, dots$cluster, drop = FALSE],
      error = function(e) {
        data <- eval(cmod$call$data,
                     envir = environment(formula(cmod)))
        stop(paste("The columns",
                   paste(setdiff(dots$cluster, colnames(data)), collapse = ", "),
                   "are missing from the covariance adjustment model dataset"),
             call. = FALSE)
      })
  } else if (inherits(dots$cluster, "data.frame")) {
    uoas <- dots$cluster
  } else if (inherits(dots$cluster, c("matrix", "list", "integer", "numeric", "factor"))) {
    uoas <- as.data.frame(dots$cluster)
  } else {
    stop(paste("If overriding `cluster` argument for meat matrix calculations,",
               "must provide a data frame, matrix, list, numeric/factor vector, or",
               "a character vector specifying column names in the covariance model dataset"))
  }

  # Replace NA's for rows not in the experimental design with a unique cluster ID
  if (ncol(uoas) == 1) {
    uoas <- uoas[, 1]
  } else {
    uoas <- Reduce(function(...) paste(..., sep = "_"), uoas)
    uoas[grepl("NA", uoas)] <- NA_character_
  }
  nuoas <- length(unique(uoas))
  nas <- is.na(uoas)
  if (any(nas)) {
    uoas[nas] <- paste0(nuoas - 1 + seq_len(sum(nas)), "*")
  }
  uoas <- factor(uoas)
  nuoas <- length(unique(uoas))

  # if the covariance model only uses one cluster, treat units as independent
  if (nuoas == 1) {
    uoas <- seq_len(length(uoas))
  }
  dots$cluster <- uoas
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

#' @title (Internal) Check the NA status of clustering columns in a dataframe
#' @details \code{.check_cluster_col_nas} checks whether entire columns in a
#' dataframe contain NA values. If a column only contains NA's, the function
#' produces a warning that the given indicating the column cannot be used for
#' clustering.
#' @param cluster_cols vector of column names
#' @param df dataframe whose columns are checked for NA's
#' @param model_type string taking either the value "covariance adjustment",
#' referring to the covariance adjustment model, or "ITT effect", referring to
#' the direct adjusted model. This value is piped into any potential warning
#' message.
#' @return A vector of column names for the valid clustering columns
#' @keywords internal
.check_cluster_col_nas <- function(cluster_cols,
                                   df,
                                   model = c("covariance adjustment", "ITT effect")) {
  model <- match.arg(model)
  all_nas <- sapply(df[, cluster_cols, drop = FALSE], function(col) all(is.na(col)))
  if (any(all_nas)) {
    all_na_cols <- names(which(all_nas))
    msg <- if (length(all_na_cols) == length(all_nas)) {
      paste("This is taken to mean the observations should be treated",
            "as IID. To avoid this warning, provide unique non-NA cluster",
            "ID's for each row.")
    } else {
      paste("Only",
            paste(setdiff(cluster_cols, all_na_cols), collapse = ", "),
            "will be used to cluster the", model, "model.")
    }
    
    cluster_cols <- setdiff(cluster_cols, all_na_cols)
    warning(paste("The columns", paste(all_na_cols, collapse = ", "),
                  "are all NA's in the", model, "model dataset.",
                  msg))
  }
  
  return(cluster_cols)
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
