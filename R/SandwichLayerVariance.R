#' @include Design.R SandwichLayer.R
NULL

#' @title Compute covariance-adjusted cluster-robust sandwich variance estimates
#' @details Supported \code{type} include:
#'
#' - \code{"CR0"}, \code{"MB0"}, \code{"HC0"} are synonyms for ...
#' - Others...
#'
#' To create your own \code{type}, simply define a function \code{.vcov_XXX}.
#' \code{type = "XXX"} will now use your method. Your method should return a
#' matrix of appropriate dimension, with \code{attribute} \code{type = "XXX"}.
#' @param x a fitted \code{DirectAdjusted} model object
#' @param type A string indicating the desired variance estimator. Currently
#'   accepts "CR0", "MB0", or "HC0"
#' @param cluster Defaults to NULL, which means unit of assignment columns
#'   indicated in the Design will be used to generate clustered covariance
#'   estimates. A non-NULL argument to `cluster` specifies a string or character
#'   vector of column names appearing in both the covariance adjustment and
#'   quasiexperimental samples that should be used for clustering covariance
#'   estimates.
#' @param ... Arguments to be passed to the internal variance estimation
#'   function.
#' @return A \eqn{2\times 2} matrix where the dimensions are given by the
#'   intercept and treatment variable terms in the ITT effect model
#' @export
#' @rdname var_estimators
vcovDA <- function(x, type = "CR0", cluster = NULL, ...) {
  if (!exists(paste0(".vcov_", type))) {
    stop(paste0("covariance function .vcov_", type,
                " not defined.\n"))
  }
  var_func <- get(paste0(".vcov_", type))
  args <- list(...)
  args$x <- x
  args$by <- cluster # if cov_adj() was not fit with a "by" argument, this is passed to .order_samples() to order rows of estfun() output
  args$cluster <- .make_uoa_ids(x, cluster, ...) # passed on to meatCL to aggregate SE's at the level given by `cluster`

  est <- do.call(var_func, args)

  return(est)
}

#' Model-based standard errors with HC0 adjustment
#' @keywords internal
#' @rdname var_estimators
.vcov_CR0 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
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
  vmat <- (1 / n) * a22inv %*% meat %*% t(a22inv)

  # NA any invalid estimates due to degrees of freedom checks
  vmat <- .check_df_moderator_estimates(vmat, x, args$cluster)

  attr(vmat, "type") <- "CR0"
  return(vmat)
}

#' @rdname var_estimators
.vcov_HC0 <- function(x, ...) {
  out <- .vcov_CR0(x, ...)
  attr(out, "type") <- "HC0"
  return(out)
}

#' @rdname var_estimators
.vcov_MB0 <- function(x, ...) {
  out <- .vcov_CR0(x, ...)
  attr(out, "type") <- "MB0"
  return(out)
}

#' Model-based standard errors with HC1 adjustment
#' @keywords internal
#' @rdname var_estimators
.vcov_CR1 <- function(x, ...) {
  args <- list(...)
  args$x <- x
  args$cadjust <- FALSE

  vmat <- do.call(.vcov_CR0, args)
  n <- length(args$cluster)
  g <- length(unique(args$cluster))
  k <- ncol(estfun(x))

  vmat <- g / (g-1) * (n-1) / (n - k) * vmat # Hansen (2022) provides this generalization

  attr(vmat, "type") <- "CR1"
  return(vmat)
}

#' @rdname var_estimators
.vcov_HC1 <- function(x, ...) {
  out <- .vcov_CR1(x, ...)
  attr(out, "type") <- "HC1"
  return(out)
}

#' @rdname var_estimators
.vcov_MB1 <- function(x, ...) {
  out <- .vcov_CR1(x, ...)
  attr(out, "type") <- "MB1"
  return(out)
}

#' @title NA vcovDA subgroup estimates that have insufficient degrees of freedom
#' @param vmat variance-covariance matrix corresponding to `model`
#' @param model \code{DirectAdjusted} object
#' @param cluster character or factor vector providing cluster ID's for the
#'  observations used to fit `model`
#' @param model_data dataframe or name corresponding to the data used to fit `model`
#' @param envir environment to get `model_data` from if it is a quote object name
#' @return `vmat` with NA's in the entries lacking sufficient degrees of freedom
.check_df_moderator_estimates <- function(vmat, model, cluster, model_data = quote(data),
                                          envir = environment(formula(model))) {
  if (!inherits(model, "DirectAdjusted")) {
    stop("`model` must be a DirectAdjusted object")
  }

  if (length(model@moderator) == 0) {
    return(vmat)
  }

  if (inherits(model_data, "name")) {
    model_data <- get(as.character(model_data), envir)
  } else if (!inherits(model_data, "data.frame")) {
    stop("`data` must be a dataframe or a quoted object name")
  }

  # For categorical moderators, count the clusters contributing to estimation
  # for each level of the moderator variable; for continuous moderators, just
  # count the number of clusters. The moderator variable (or any level of the
  # moderator variable) must have at least three clusters contributing to
  # estimation for valid SE estimation.
  mod_vars <- model.matrix(as.formula(paste0("~-1+", model@moderator)), model_data)
  if (ncol(mod_vars) > 1) {
    mod_counts <- sweep(rowsum(mod_vars, cluster), 1,
                        rowsum(rep(1, nrow(mod_vars)), cluster), FUN = "/")
    valid_mods <- colSums(mod_counts != 0) > 2
  } else {
    valid_mods <- setNames(length(unique(cluster)) > 2, model@moderator)
  }
    

  # Replace SE's for moderator variable/levels with <= 2 clusters with NA's
  if (any(!valid_mods)) {
    invalid_mods <- gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", names(valid_mods)[!valid_mods]))
    warning(paste("The following subgroups do not have sufficient degrees of",
                  "freedom for standard error estimates and will be returned",
                  "as NA:",
                  paste(names(valid_mods)[!valid_mods], collapse = ", ")),
            call. = FALSE)
    # find cells of the covariance matrix that correspond to txt/mod group
    # interactions
    if (!(valid_mods[1])) invalid_mods <- c("\\(Intercept\\)", invalid_mods)
    dims_to_na <- which(grepl(paste0("(", paste(invalid_mods, collapse = "|"), ")"),
                              colnames(vmat),
                              perl = TRUE))

    vmat[dims_to_na, ] <- NA_real_
    vmat[, dims_to_na] <- NA_real_
  }

  return(vmat)
}

#' @details The \bold{B12 block} is the covariance matrix of the cluster-level
#'   estimating equations for the covariance adjustment and ITT effect models. It
#'   has a row for each term in the covariance adjustment model and a column for each term
#'   in the ITT effect model. For any row that does not appear in both
#'   the experimental design and the covariance adjustment model data, its contribution to
#'   this matrix will be 0. Thus, if there is no overlap between the two
#'   datasets, this will return a matrix of 0's.
#' @param x a fitted \code{DirectAdjusted} model
#' @param ... arguments to pass to internal functions
#' @return \code{.get_b12()}: A \eqn{p\times 2} matrix where the number of rows
#'   are given by the number of terms in the covariance adjustment model and the number of
#'   columns correspond to intercept and treatment variable terms in the ITT
#'   effect model
#' @keywords internal
#' @noRd
.get_b12 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  if (x@.S3Class[1] != "lm") {
    stop("x must be an `lm` object")
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
        unique(x@Design@structure[[col]])[
          match(wide_frame[[col]], unique(x@Design@structure[[col]]), incomparables = NA)
        ]
      })
    )
  } else {
    stop(paste("If overriding `cluster` argument for meat matrix calculations,",
               "must provide a character vector specifying column names that",
               "exist in both the ITT effect and covariance model datasets"))
  }

  Q_uoas <- stats::expand.model.frame(x, cluster_cols, na.expand = TRUE)[cluster_cols]
  Q_uoas <- apply(Q_uoas, 1, function(...) paste(..., collapse = "_"))
  C_uoas <- apply(keys[, cluster_cols, drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  C_uoas_in_Q <- C_uoas %in% unique(Q_uoas)

  message(paste(sum(C_uoas_in_Q),
                "rows in the covariance adjustment model",
                "data joined to the ITT effect model data\n"))

  # Check number of overlapping clusters n_QC; if n_QC <= 1, return 0 (but
  # similarly to .get_b11(), throw a warning when only one cluster overlaps)
  uoas_overlap <- length(unique(C_uoas[C_uoas_in_Q]))
  if (uoas_overlap <= 1) {
    if (uoas_overlap == 1) {
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
  cmod_estfun <- sandwich::estfun(sl@fitted_covariance_model)[C_uoas_in_Q,]
  cmod_aggfun <- ifelse(dim(cmod_estfun)[2] > 1, colSums, sum)
  cmod_eqns <- Reduce(rbind, by(cmod_estfun, C_uoas[C_uoas_in_Q], cmod_aggfun))

  # get rows from overlapping clusters in experimental data
  Q_uoas_in_C <- Q_uoas %in% unique(C_uoas[C_uoas %in% Q_uoas])
  damod_estfun <- sandwich::estfun(as(x, "lm"))[Q_uoas_in_C, , drop = FALSE]
  damod_aggfun <- ifelse(dim(damod_estfun)[2] > 1, colSums, sum)
  damod_eqns <- Reduce(rbind, by(damod_estfun, Q_uoas[Q_uoas_in_C], damod_aggfun))

  return(crossprod(cmod_eqns, damod_eqns))
}

#' @title (Internal) Compute variance blocks
#' @param x a fitted \code{DirectAdjusted} model
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

#' @param x a fitted \code{DirectAdjusted} model
#' @param ... arguments to pass to internal functions
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
#' @noRd
.get_b22 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  nq <- nrow(sandwich::estfun(x))

  # Create cluster ID matrix depending on cluster argument (or its absence)
  dots <- list(...)
  if (is.null(dots$cluster)) {
    uoas <- propertee::.expand.model.frame.DA(x,
                         var_names(x@Design, "u"))[, var_names(x@Design, "u"),
                                                   drop = FALSE]
  } else if (inherits(dots$cluster, "character")) {
    uoas <- tryCatch(
      .expand.model.frame.DA(x, dots$cluster)[,
                                                        dots$cluster,
                                                        drop = FALSE],
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
  dots$x <- as(x, "lm")

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
  if(nc==0) nc <- length(cmod$resid)

  out <- sandwich::bread(cmod) / nc
  return(out)
}

#' @param x a fitted \code{DirectAdjusted} model
#' @param ... arguments to pass to internal functions
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
#' @noRd
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
    keys <- sl@keys
    cluster_cols <- setdiff(colnames(keys), "in_Q")
  } else if (inherits(dots$cluster, "character")) {
    cluster_cols <- dots$cluster
    keys <- tryCatch(
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
    nas <- rowSums(is.na(keys[, cluster_cols, drop = FALSE]))
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
  uoas <- Reduce(function(...) paste(..., sep = "_"), keys[, cluster_cols, drop = FALSE])

  any_nas <- apply(keys[, cluster_cols, drop = FALSE], 1, function(x) any(is.na(x)))
  uoas[any_nas] <- NA_character_

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
               "to propagate covariance adjustment model error"))
  }

  # Get contribution to the estimating equation from the ITT effect model
  w <- if (is.null(x$weights)) 1 else x$weights

  damod_mm <- stats::model.matrix(
    formula(x), stats::model.frame(x, na.action = na.pass))
  msk <- (apply(!is.na(sl@prediction_gradient), 1, all) &
            apply(!is.na(damod_mm), 1, all))

  out <- crossprod(damod_mm[msk, x$qr$pivot[1L:x$rank], drop = FALSE] * w,
                   sl@prediction_gradient[msk, , drop = FALSE])

  return(out)
}
