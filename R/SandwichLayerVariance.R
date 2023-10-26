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

  # For each moderator variable (whether it's been dichotomized or not), count
  # the number of clusters with at least one member of each value
  mod_vars <- model.matrix(as.formula(paste0("~-1+", model@moderator)),
                           model_data)
  mod_counts <- apply(
    mod_vars,
    2,
    function(col) {
      n_vals <- length(unique(col))
      if (n_vals == 2) {
        tapply(cluster, col, function(cluster_ids) length(unique(cluster_ids)))
      } else {
        3 # necessarily enough degrees of freedom if there are at least 3 values
      }
    },
    simplify = FALSE)

  valid_mods <- vapply(mod_counts,
                       function(counts) all(counts > 2),
                       logical(1L))
  if (any(!valid_mods)) {
    invalid_mods <- names(valid_mods)[!valid_mods]
    warning(paste("The following subgroups do not have sufficient degrees of",
                  "freedom for standard error estimates and will be returned",
                  "as NA:",
                  paste(invalid_mods, collapse = ", ")),
            call. = FALSE)
    # find cells of the covariance matrix that correspond to txt/mod group
    # interactions
    invalid_cols <- paste0(paste0(var_names(model@Design, "t"), "."), "_", invalid_mods)
    dims_to_na <- apply(
      vapply(invalid_cols, grepl, logical(length(colnames(vmat))),
             colnames(vmat), fixed = TRUE),
      1,
      any
    )

    if (all(!dims_to_na)) {
      warning(paste("Could not find dimensions of `vmat` corresponding to",
                    "degenerate standard error estimates. Degenerate standard",
                    "error estimates will not be returned as NA"))
    }

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


#' Design-based standard errors with HC0 adjustment
#' @param x a fitted \code{DirectAdjusted} model
#' @keywords internal
#' @rdname var_estimators
.vcov_DB0 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }
  
  #if (!x@lmitt_fitted){
  #  stop("x must have been fitted using lmitt.formula")
  #}
  # x@lmitt_fitted is false if someone created x using as.lmitt
  
  args <- list(...)
  if ("type" %in% names(args)) {
    stop(paste("Cannot override the `type` argument for meat",
               "matrix computations"))
  }
  args$x <- x
  args$db <- TRUE
  n <- length(args$cluster)
  
  if (x@absorbed_intercepts) {
    a22inv <- sandwich::bread(x)
    meat <- do.call(sandwich::meatCL, args)
    
    vmat <- (1 / n) * a22inv %*% meat %*% a22inv
    name <- grep("z", colnames(x$model))
    vmat <- as.matrix(vmat[name, name])
  }
  else {
    if (!is.null(x$call$offset)){
      vmat <- .get_DB_covadj(x)
      if (is.na(vmat))
        stop(paste("Design-based standard errors cannot be computed for ITT effect",
                   "models with covariance adjustment if small strata are present"))
    }
    if (length(x@moderator) > 0){
      stop(paste("Design-based standard errors cannot be computed for ITT effect",
                 "models with moderators"))
    }
    ainv <- .get_DB_a_inverse(x)
    meat <- .get_DB_meat(x)
    vmat <- as.matrix((ainv %*% meat %*% t(ainv))[3,3])
    if (is.na(vmat)) # if small strata are present
      vmat <- .get_DB_small_strata(x)
  }
  name <- colnames(x$model)[grep("z", colnames(x$model))]
  colnames(vmat) <- name
  rownames(vmat) <- name
  
  attr(vmat, "type") <- "DB0"
  return(vmat)
}

#' @title (Internal) Compute design-based variance estimate
#' @param x a fitted \code{DirectAdjusted} model
#' @details Calculate bread matrix for design-based variance estimate for 
#'  ITT effect models with covariance adjustment and without absorbed effects
#'  when only moderate or large strata are present
#' @return design-based estimation of variance 
#' @keywords internal
.get_DB_covadj <- function(x, ...){
  a11inv <- .get_a11_inverse(x)
  a21 <- .get_a21(x)
  a22inv <- .get_a22_inverse(x)
  C <- matrix(c(1,1,0,1), nrow = 2, byrow = TRUE)
  
  bread1 <- a22inv %*% a21 %*% a11inv
  bread2 <- - a22inv %*% C
  
  return(0)
}

#' @title (Internal) Compute design-based variance blocks
#' @param x a fitted \code{DirectAdjusted} model
#' @details Calculate bread matrix for design-based variance estimate for 
#'  ITT effect models without covariance adjustment and without absorbed effects
#'  when only moderate or large strata are present
#' @return inverse of bread matrix 
#' @keywords internal
.get_DB_a_inverse <- function(x, ...){
  res <- .aggregate_individuals(x)
  data <- res[[1]]
  bid <- data[, res[[2]]] # block ids
  zobs <- data[, res[[3]]] # observed zs
  
  nbk <- sapply(c(0,1), function(z) as.vector(table(data[zobs == z, res[[2]]])))
  pbk <- nbk / rowSums(nbk) # assignment probabilities, B by K
  
  A <- matrix(c(rep(0,6), 1,-1,1), nrow = 3, byrow = TRUE)
  ws_agg <- aggregate(data$.w, by = list(bid), FUN = sum)[,2]
  A[1,1] <- sum(pbk[,1] * ws_agg)
  A[2,2] <- sum(pbk[,2] * ws_agg)
  return(solve(A))
}

#' @title (Internal) Compute design-based variance blocks
#' @param x a fitted \code{DirectAdjusted} model
#' @details Calculate meat matrix for design-based variance estimate for 
#'  ITT effect models without covariance adjustment and without absorbed effects
#'  when only moderate or large strata are present
#' @return meat matrix 
#' @keywords internal
.get_DB_meat <- function(x, ...){
  res <- .aggregate_individuals(x)
  data <- res[[1]]
  
  ws <- data$.w
  yobs <- data$.wy # observed ys
  bid <- data[, res[[2]]] # block ids
  zobs <- data[, res[[3]]] # observed zs
  
  rho <- c(sum((1-zobs) * ws * yobs) / sum((1-zobs) * ws),
           sum(zobs * ws * yobs) / sum(zobs * ws))
  
  nbk <- sapply(c(0,1), function(z) as.vector(table(data[zobs==z, res[[2]]])))
  # number of units in each block and treatment, B by K
  pbk <- nbk / rowSums(nbk) # assignment probabilities, B by K
  delbk <- (nbk - 1) / (rowSums(nbk) - 1) * pbk # second assignment probabilities
  
  B <- matrix(0, nrow = 3, ncol = 3)
  wy0_agg <- aggregate((1-zobs)*ws^2*(yobs-rho[1])^2, by = list(bid), FUN = sum)[,2]
  wy1_agg <- aggregate(zobs*ws^2*(yobs-rho[2])^2, by = list(bid), FUN = sum)[,2]
  
  wyy0_agg <- aggregate((ws*(yobs-rho[1]))[zobs == 0], 
                        by = list(bid[zobs == 0]), FUN = .prod_sum)[,2]
  wyy1_agg <- aggregate((ws*(yobs-rho[2]))[zobs == 1], 
                        by = list(bid[zobs == 1]), FUN = .prod_sum)[,2]
  
  B[1,1] <- sum((1-pbk[,1]) * wy0_agg) + sum((1-pbk[,1]^2/delbk[,1]) * wyy0_agg)
  B[2,2] <- sum((1-pbk[,2]) * wy1_agg) + sum((1-pbk[,2]^2/delbk[,2]) * wyy1_agg)
  B[1,2] <- - (B[1,1] + B[2,2]) / 2
  B[2,1] <- B[1,2]
  return(B)
}

#' Aggregate individual-level weights and outcomes to cluster-level
#' @param x a fitted \code{DirectAdjusted} model
#' @return a list of a data frame of weights, outcomes, treatments, and block ids;
#'    the name of the treatment id column; the name of the block id column
#' @keywords internal
.aggregate_individuals <- function(x, ...){
  ws <- if (is.null(x$weights)) 1 else x$weights
  data_temp <- x$call$data
  name_y <- as.character(x$terms[[2]]) # the column of y
  data_temp <- cbind(data_temp, .w = ws, .w0 = ws / ate(x@Design, data=x$call$data),
                     .wy = ws * data_temp[, name_y])
  
  design_obj <- x@Design
  name_trt <- colnames(design_obj@structure)[design_obj@column_index == "t"]
  name_blk <- colnames(design_obj@structure)[design_obj@column_index == "b"]
  name_clu <- colnames(design_obj@structure)[design_obj@column_index == "u"]
  
  data_temp <- .combine_block_ID(data_temp, name_blk)
  name_blk <- name_blk[1]
  
  eq_clu <- paste(name_clu, collapse = " + ")
  form <- paste("cbind(", name_trt, ", ", name_blk, ") ~ ", eq_clu, sep = "")
  form2 <- paste("cbind(.wy, .w, .w0) ~ ", eq_clu, sep = "")
  
  # aggregate z and bid by mean, aggregate w, w0, and w*y by sum
  data_aggr <- cbind(
    do.call("aggregate", list(as.formula(form), FUN = mean, data = data_temp)),
    do.call("aggregate", list(as.formula(form2), FUN = sum, data = data_temp))
  )
  data_aggr$.wy <- data_aggr$.wy / data_aggr$.w
  data_aggr <- data_aggr[order(data_aggr[, name_blk]), ]
  return(list(data = data_aggr, block = name_blk, z = name_trt))
}

#' Combine multiple block IDs to one column 
#' @details
#' the returned data frame has a column named as ids[1] that
#' contains unique numbers based on the combinations of the values
#' in the multiple block ID columns
#' @param df a data frame 
#' @param ids a vector of block IDs, column names of df
#' @return data frame df with a column that contains unique block number IDs
#' @keywords internal
.combine_block_ID <- function(df, ids){
  df[,ids[1]] <- apply(df[, ids, drop = FALSE], 1, paste, collapse = "_")
  unique_ids <- data.frame(unique(df[,ids[1]]))
  colnames(unique_ids) <- ids[1]
  unique_ids$.ID <- seq_len(nrow(unique_ids))
  
  df <- merge(df, unique_ids, by = ids[1], all.x = TRUE)
  df[,ids[1]] <- df$.ID
  df$.ID <- NULL
  return(df)
}

#' Calculate sum of x[i] * x[j] with i not equal to j
#' @param x a numeric vector
#' @keywords internal
.prod_sum <- function(x){
  return(sum(x * sum(x)) - sum(x^2))
}

#' Design-based variance estimate for ITT effect models 
#' without covariance adjustment and without absorbed effects
#' when small strata are present
#' @param x a fitted \code{DirectAdjusted} model
#' @return the design-based variance estimate
#' @keywords internal
.get_DB_small_strata <- function(x, ...){
  res <- .aggregate_individuals(x)
  data <- res[[1]]
  block <- res[[2]]
  
  ws <- data$.w
  yobs <- data$.wy  # observed ys
  zobs <- data[, res[[3]]]  # observed zs
  bid <- data[, block] # block ids
  
  nbk <- sapply(c(0,1), function(z) as.vector(table(data[zobs == z, block])))
  # number of units in each block and treatment, B by K
  nbk_all <- nbk[bid, ]  # n by K
  B <- length(unique(bid))  # number of blocks
  K <- ncol(nbk)  # number of treatments
  
  gammas <- nbk_all * 
    matrix(ws, nrow = nrow(nbk_all), ncol = ncol(nbk_all), byrow = FALSE)
  # gamma for variance estimation, n by K
  gamsbk <- list()  # s^2_b,j, sample variance
  
  nu1 <- matrix(nrow=B, ncol=K-1)  # nu1_b,0k
  varest <- matrix(nrow=1, ncol=K-1)  # variance estimators
  
  for (k in 1:K){
    indk <- zobs == (k-1)
    thetak <- sum(ws[indk] * yobs[indk]) / sum(ws[indk])  # ratio estimate rho_z
    gammas[indk,k] <- gammas[indk,k] * (yobs[indk] - thetak)
    gamsbk[[k]] <- aggregate(gammas[indk,k], by = list(data[indk,block]), FUN = var)
  }
  gamsbk <- merge(gamsbk[[1]], gamsbk[[2]], by = "Group.1")[,2:(K+1)]  # B by K
  gamsbk[is.na(gamsbk)] <- 0
  
  for (k in 2:K){
    for (b in 1:B){
      indb0 <- (zobs == 0) & (bid == b)
      indbk <- (zobs == k-1) & (bid == b)  # units in block b, treatment k
      
      if (gamsbk[b,1] == 0 | gamsbk[b,k] == 0){ # small block
        nu1[b,k-1] <- 
          .get_ms_contrast(gammas[indbk,k], gammas[indb0,1]) -
          sum((nbk[b,c(1,k)]-1) / nbk[b,c(1,k)] * gamsbk[b,c(1,k)])
      }
      else{ # large block
        nu1[b,k-1] <- sum(gamsbk[b,c(1,k)] / nbk[b,c(1,k)])
      }
    }
    varest[1,k-1] <- sum(nu1[,k-1]) / sum(data$.w0)^2  # variance estimators
  }
  return(varest)
}

#' Calculate the mean of squared contrasts of entries of two vectors
#' @param a a numeric vector
#' @param b a numeric vector
#' @return the mean of squared contrasts between elements of the input vectors
#' @keywords internal
.get_ms_contrast <- function(a,b){
  t = 0
  for (i in 1:length(b))
    t = t + sum((a - b[i])^2)
  return (t/length(a)/length(b))
}

#' This function calculates grave{phi}
#' @keywords internal
#' @param x a fitted \code{DirectAdjusted} model
#' @param ... arguments to pass to internal functions
.get_phi_tilde <- function(x, ...){
  ws <- x$weights
  n <- length(ws) # number of units in Q
  design_obj <- x@Design
  df <- x$call$data
  
  # treatment assignments
  name_trt <- colnames(design_obj@structure)[design_obj@column_index == "t"]
  assignment <- df[[name_trt]]
  k <- length(unique(assignment))
  z_ind <- sapply(c(0,1), function(val) as.integer(assignment == val))
  
  # stratum ids 
  name_blk <- colnames(design_obj@structure)[design_obj@column_index == "b"]
  df <- .combine_block_ID(df, name_blk)
  name_blk <- name_blk[1]
  stratum <- df[[name_blk]]
  s <- length(unique(stratum)) # number of blocks
  b_ind <- sapply(unique(stratum), function(val) as.integer(stratum == val)) 
  
  # nuisance parameters p
  wb <- matrix(replicate(s, ws), ncol = s) * b_ind
  p <- t(z_ind) %*% wb / (matrix(1, nrow = k, ncol = n) %*% wb)
  p1 <- p[2, ]
  
  # calculate phi tilde
  phitilde <- matrix(nrow = n, ncol = s)
  for (j in 1:s){
    phitilde[,j] <- ws * (z_ind[,2] - p1[j]) * b_ind[,j]
  }
  return(phitilde)
}

#' This function calculates the product of two matrices A_{pp}^{-1} A_{tau p}^T
#' @keywords internal
#' @param x a fitted \code{DirectAdjusted} model
#' @param ... arguments to pass to internal functions
#' @return An \eqn{s\times k} matrix A_{pp}^{-1} A_{tau p}^T
.get_appinv_atp <- function(x, ...){
  ws <- x$weights
  # estimated treatment effect tau1 <- x$coefficients
  n <- length(ws) # number of units in Q(?)
  design_obj <- x@Design
  df <- x$call$data
  
  name_trt <- colnames(design_obj@structure)[design_obj@column_index == "t"]
  assignment <- df[[name_trt]]
  k <- length(unique(assignment))
  
  # stratum ids 
  name_blk <- colnames(design_obj@structure)[design_obj@column_index == "b"]
  df <- .combine_block_ID(df, name_blk)
  name_blk <- name_blk[1]
  stratum <- df[[name_blk]]
  s <- length(unique(stratum)) # number of blocks
  b_ind <- sapply(unique(stratum), function(val) as.integer(stratum == val)) 
  
  app <- list()
  atp <- list()
  for (i in 1:s){
    atp[[i]] <- ws * (x$residuals) * b_ind[,i]
    # atp[[i]] <- .get_upsilon(x) / term2 * b_ind[,i]
    app[[i]] <- ws * b_ind[,i]
  }
  # w_ipv <- ate(x@Design, data = x$call$data) inverse probability weights
  
  mat <- matrix(0, nrow = (k-1)*s, ncol = k-1)
  for (i in 1:s){
    # mat[i,1] <- sum(app[[i]] * w_ipv) / sum(w_ipv)
    mat[i,1] <- sum(atp[[i]]) / sum(app[[i]])
  }
  return(mat)
}
