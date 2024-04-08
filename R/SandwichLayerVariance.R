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
  args$cluster <- .make_uoa_ids(x, vcov_type = substr(type, 1, 2), cluster, ...) # passed on to meatCL to aggregate SE's at the level given by `cluster`

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

  bread. <- sandwich::bread(x)
  meat. <- do.call(sandwich::meatCL, args)
  vmat <- (1 / n) * bread. %*% meat. %*% t(bread.)

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
#' @keywords internal
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
  mod_vars <- mod_vars[rownames(mod_vars) %in% rownames(stats::model.frame(model)),,drop=FALSE]
  if (ncol(mod_vars) > 1) {
    mod_counts <- sweep(rowsum(mod_vars, cluster), 1,
                        rowsum(rep(1, nrow(mod_vars)), cluster), FUN = "/")
    valid_mods <- colSums(mod_counts != 0) > 2
  } else {
    valid_mods <- stats::setNames(length(unique(cluster)) > 2, model@moderator)
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
#' @details \eqn{A_{22}^{-1}} is the inverse observed Fisher information of the
#' ITT effect estimating equations scaled by \eqn{n_{\mathcal{Q}}}.
#' @param x a fitted \code{DirectAdjusted} object
#' @param ... arguments passed to methods
#' @return \code{.get_a22_inverse()}/\code{.get_tilde_a22_inverse()}: A
#' \eqn{k\times k} matrix where k denotes the number of parameters in the ITT
#' effect model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a22_inverse <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }
  
  out <- utils::getS3method("bread", "lm")(x)

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
    uoas <- .expand.model.frame.DA(x,
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

#' @details \eqn{A_{11}^{-1}} is the inverse of the gradient of the covariance
#'   adjustment model estimating equations scaled by \eqn{n_{\mathcal{C}}^{-1}}.
#' @return \code{.get_a11_inverse()}: A \eqn{p\times p} matrix where the
#'   \eqn{p} is the dimension of the covariance adjustment model including an
#'   intercept
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

  out <- sandwich::bread(sl@fitted_covariance_model)
  
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

#' @details \eqn{A_{21}} is the gradient of the ITT effect estimating equations
#'   scaled by \eqn{n_{\mathcal{Q}}^{-1}} taken with respect to the covariance
#'   adjustment model parameters. This matrix is the crossproduct of the
#'   prediction gradient for the units of observation in \eqn{\mathcal{Q}} and
#'   the model matrix of the ITT effect estimating eqations.
#' @param x a fitted \code{DirectAdjusted} model
#' @return \code{.get_a21()}/\code{.get_tilde_a21()}: A \eqn{k\times p} matrix
#'   where the number of rows are given by the dimension of the ITT effect
#'   estimating equations and the number of columns are given by the number of
#'   terms in the covariance adjustment model
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
  # scale by nq
  nq <- sum(msk)

  return(out / nq)
}

##' @details \eqn{\tilde{A}_{22}^{-1}} is the inverse observed Fisher
##' information of the ITT effect estimating equations scaled by \eqn{n}. This
##' function wraps around the function \code{.get_a22_inverse()} that produces
##' \eqn{A_{22}^{-1}}, where \eqn{A_{22}=\frac{n}{n_{\mathcal{Q}}}\tilde{A}_{22}}.
##' @inheritDotParams .get_a22_inverse
##' @inherit .get_a22_inverse return
##' @keywords internal
##' @rdname sandwich_elements_calc
.get_tilde_a22_inverse <- function(x, ...) {
  out <- .get_a22_inverse(x, ...)
  
  if (!inherits(ca <- x$model$`(offset)`, "SandwichLayer")) {
    return(out)
  }

  nq <- nrow(stats::model.frame(x))
  nc_not_q <- sum(!ca@keys$in_Q)
  n <- nq + nc_not_q
  
  out <- out * n / nq
  
  return(out)
}

##' @details \eqn{\tilde{A}_{21}} is the gradient of the ITT effect estimating
##'   equations scaled by \eqn{n^{-1}} taken with respect to the covariance
##'   adjustment model parameters. This function wraps around \code{.get_a21()},
##'   which produces \eqn{A_{21}}, where \eqn{A_{21} = \frac{n_{\mathcal{Q}}}{n}
##'   \tilde{A}_{21}}.
##' @inheritDotParams .get_a21
##' @inherit .get_a21 return
##' @keywords internal
##' @rdname sandwich_elements_calc
.get_tilde_a21 <- function(x) {
  out <- .get_a21(x)
  
  nq <- nrow(stats::model.frame(x))
  sl <- x$model$`(offset)`
  nc_not_q <- sum(!sl@keys$in_Q)
  n <- nq + nc_not_q
  
  out <- nq / n * out
}

#' @title Design-based standard errors with HC0 adjustment
#' @param x a fitted \code{DirectAdjusted} model
#' @keywords internal
#' @rdname var_estimators
.vcov_DB0 <- function(x, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }
  
  if (!x@lmitt_fitted){
    # x@lmitt_fitted is false if someone created x using as.lmitt
    stop("x must have been fitted using lmitt.formula")
  }
  
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
      vmat <- .get_DB_covadj_se(x)
      if (is.na(vmat))
        stop(paste("Design-based standard errors cannot be computed for ITT effect",
                   "models with covariance adjustment if small strata are present"))
    }
    else if (length(x@moderator) > 0){
      stop(paste("Design-based standard errors cannot be computed for ITT effect",
                 "models with moderators"))
    }
    else {
      vmat <- .get_DB_se(x)
    }
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
.get_DB_covadj_se <- function(x, ...){
  design_obj <- x@Design
  data <- x$call$data
  name_clu <- colnames(design_obj@structure)[design_obj@column_index == "u"]
  cid <- .combine_block_ID(data, name_clu)[, name_clu[1]]
  
  res <- .aggregate_individuals(x)
  data <- res[[1]]
  bid <- data[, res[[2]]] # block ids
  zobs <- data[, res[[3]]] # observed zs
  nbk <- design_table(design=design_obj, x="treatment",y="block")
  
  name_y <- as.character(x$terms[[2]]) # name of the outcome column
  X1 <- model.matrix(x$call$offset@fitted_covariance_model) # design matrix
  p <- ncol(X1)
  n <- nrow(X1)
  
  wc <- x$call$offset@fitted_covariance_model$weight
  if (is.null(wc)) wc <- 1
  X1 <- matrix(rep(wc * (x$call$data[,name_y] - x$offset), p),
               ncol = p) * X1  # n by p, wic * residual * xi
  wi <- x$weights
  X2 <- wi * x$residuals  # n by 1, wi[z] * (residual - rhoz) * z
  
  XX <- cbind(X1, X2)
  XX <- aggregate(XX, by = list(cid), FUN = sum)[, 2:(p+2)]
  const <- sqrt(nbk[,1] * nbk[,2] / rowSums(nbk))
  XX <- sweep(XX, 1, const[bid], '*')
  
  V00 <- .cov_mat_est(XX[zobs==0,], bid[zobs==0])
  V11 <- .cov_mat_est(XX[zobs==1,], bid[zobs==1])
  V01 <- .cov01_est(XX, zobs, bid)
  
  a11inv <- .get_a11_inverse(x)
  a21 <- .get_a21(x)
  a22inv <- .get_a22_inverse(x)
  C <- matrix(c(1,1,0,1), nrow = 2, byrow = TRUE)
  
  bread1 <- a22inv %*% a21 %*% a11inv / n
  bread2 <- - a22inv %*% C / n
  
  signs1 <- ifelse(t(t(bread1[2,])) %*% bread1[2,] > 0, 1, 0)
  signs2 <- ifelse(t(t(bread1[2,])) %*% bread2[2,] > 0, 1, 0)
  signs3 <- ifelse(t(t(bread2[2,])) %*% bread2[2,] > 0, 1, 0)
  
  idl <- (p+2):(2*p+1)
  meat1u <- V00[1:p, 1:p] + V01[1:p, 1:p] + t(V01[1:p, 1:p]) + V11[1:p, 1:p]
  meat1l <- V00[idl, 1:p] + V01[idl, 1:p] + t(V01[idl, 1:p]) + V11[idl, 1:p]
  term1 <- bread1 %*% (meat1u * signs1 + meat1l * (1 - signs1)) %*% t(bread1)
  
  meat2u <- cbind(V00[1:p, p+1], V01[1:p, p+1]) + cbind(V01[p+1, 1:p], V11[1:p, p+1])
  meat2l <- cbind(V00[idl, p+1], V01[idl, p+1]) + cbind(V01[2*p+2, 1:p], V11[idl, p+1])
  term2 <- bread2 %*% t(meat2u * signs2 + meat2l * (1 - signs2)) %*% t(bread1)
  
  meat3u <- matrix(c(V00[p+1, p+1], V01[p+1, p+1], V01[p+1, p+1], V11[p+1, p+1]), ncol = 2)
  meat3l <- matrix(c(V00[2*p+2, p+1], V01[2*p+2, p+1], V01[2*p+2, p+1], V11[2*p+2, p+1]), ncol = 2)
  term3 <- bread2 %*% (meat3u * signs3 + meat3l * (1 - signs3)) %*% t(bread2)
  
  vmat <- term1 + 2*term2 + term3
  return(as.matrix(vmat[2,2]))
}

#' @title (Internal) Design-based variance estimate helper function
#' @keywords internal
.cov_mat_est <- function(XXz, bidz){
  cov0 <- tapply(1:nrow(XXz), bidz, function(s) cov(XXz[s,]))
  covuu <- tapply(1:nrow(XXz), bidz, function(s) .add_vec(XXz[s,]))
  covll <- tapply(1:nrow(XXz), bidz, function(s) .add_vec(XXz[s,], upper=FALSE))
    
  cov0u <- lapply(1:length(cov0), 
                  function(s) if (is.na(cov0[[s]][1,1])) covuu[[s]] else cov0[[s]])
  cov0l <- lapply(1:length(cov0), 
                  function(s) if (is.na(cov0[[s]][1,1])) covll[[s]] else cov0[[s]])
  
  V00u <- Reduce('+', cov0u)
  V00l <- Reduce('+', cov0l)
  return(rbind(V00u, V00l))
}

#' @title (Internal) Design-based variance estimate helper function
#' @keywords internal
.add_mat_diag <- function(A, B){
  d <- nrow(A)
  A <- matrix(rep(diag(A), d), nrow = d, byrow = FALSE)
  B <- matrix(rep(diag(B), d), nrow = d, byrow = TRUE)
  #return(sqrt(abs(A * B)))
  return((A + B) / 2)
}

#' @title (Internal) Design-based variance estimate helper function
#' @keywords internal
.add_vec <- function(a, upper = TRUE){
  if (nrow(a) > 1) return(0)
  a <- as.numeric(a)
  d <- length(a)
  A <- matrix(rep(a, d), nrow = d, byrow = FALSE)
  B <- matrix(rep(a, d), nrow = d, byrow = TRUE)
  #if (upper) return(A*B + sqrt(abs(A*B)))
  #else return(A*B - sqrt(abs(A*B)))
  if (upper) return((A + B)^2 / 2)
  else return(- (A - B)^2 / 2)
}

#' @title (Internal) Design-based variance estimate helper function
#' @keywords internal
.cov01_est <- function(XX, zobs, bid){
  cov0 <- tapply(1:nrow(XX[zobs==0,]), bid[zobs==0], function(s) cov(XX[zobs==0,][s,]))
  cov1 <- tapply(1:nrow(XX[zobs==1,]), bid[zobs==1], function(s) cov(XX[zobs==1,][s,]))
  cov01 <- lapply(1:length(cov0), function(s) .add_mat_diag(cov0[[s]], cov1[[s]]))
  
  cov01u <- lapply(1:length(cov0), 
                   function(s) if (!is.na(cov01[[s]][1,1])) cov01[[s]]
                   else .add_mat_sqdif(XX, zobs, bid, s))
  cov01l <- lapply(1:length(cov0), 
                   function(s) if (!is.na(cov01[[s]][1,1])) -cov01[[s]]
                   else -.add_mat_sqdif(XX, zobs, bid, s, upper=FALSE))
  V01u <- Reduce('+', cov01u)
  V01l <- Reduce('+', cov01l)
  return(rbind(V01u, V01l))
}

#' @title (Internal) Design-based variance estimate helper function
#' @keywords internal
.add_mat_sqdif <- function(X, zobs, bid, b, upper = TRUE){
  A <- X[zobs==0 & bid==b, ]
  if (upper) B <- X[zobs==1 & bid==b, ] else B <- -X[zobs==1 & bid==b, ]
  cov01u <- matrix(0, nrow=ncol(A), ncol=ncol(A))
  for (j in 1:ncol(B)){
    for (h in 1:nrow(B)){
      cov01u[,j] <- cov01u[,j] + colSums((A + B[h,j])^2)
    }
  }
  cov01u <- cov01u / nrow(A) / nrow(B)
  return(cov01u)
}

#' Design-based variance estimate for ITT effect models 
#' without covariance adjustment and without absorbed effects
#' @param x a fitted \code{DirectAdjusted} model
#' @return the design-based variance estimate
#' @keywords internal
.get_DB_se <- function(x, ...){
  res <- .aggregate_individuals(x)
  data <- res[[1]]
  block <- res[[2]]
  
  ws <- data$.w
  yobs <- data$.wy  # observed ys
  zobs <- data[, res[[3]]]  # observed zs
  bid <- data[, block] # block ids
  
  nbk <- design_table(design=x@Design, x="treatment",y="block")
  # number of units by block and treatment, B by K
  small_blocks <- apply(nbk, 1, min) == 1 # small block indicator
  
  nb <- rowSums(nbk)
  rho <- c(sum(ws * yobs * (1-zobs)) / sum(ws * (1-zobs)),
             sum(ws * yobs * zobs) / sum(ws * zobs))
  
  gammas <- (nbk[bid,1]*(1-zobs) + nbk[bid,2]*zobs) * ws # pseudo outcomes
  gamsbk <- list()  # s^2_b,j, sample variance
  for (k in 1:2){
    indk <- zobs == (k-1)
    gammas[indk] <- gammas[indk] * (yobs[indk] - rho[k])
    gamsbk[[k]] <- aggregate(gammas[indk], by = list(data[indk,block]), FUN = var)
  }
  gamsbk <- merge(gamsbk[[1]], gamsbk[[2]], by = "Group.1")[,2:3]
  gamsbk[is.na(gamsbk)] <- 0
  gamsb <- aggregate(gammas, by = list(bid), FUN = var)[,2] # B vector
  
  nu1 <- rowSums(gamsbk / nbk)
  nu2 <- 2 / nbk[,1] / nbk[,2] * choose(nb,2) * gamsb - 
    (1/nbk[,1] + 1/nbk[,2]) * ((nbk[,1]-1) *gamsbk[,1] + (nbk[,2]-1) *gamsbk[,2])
  varest <- (sum(nu1[!small_blocks]) + sum(nu2[small_blocks])) / sum(data$.w0)^2
  return(as.matrix(varest))
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
  
  nbk <- design_table(design=x@Design, x="treatment",y="block")
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
  # Calculate sum of x[i] * x[j] with i not equal to j
  # x is a numeric vector
  .prod_sum <- function(x){
    return(sum(x * sum(x)) - sum(x^2))
  }
  res <- .aggregate_individuals(x)
  data <- res[[1]]
  
  ws <- data$.w
  yobs <- data$.wy # observed ys
  bid <- data[, res[[2]]] # block ids
  zobs <- data[, res[[3]]] # observed zs
  
  rho <- c(sum((1-zobs) * ws * yobs) / sum((1-zobs) * ws),
           sum(zobs * ws * yobs) / sum(zobs * ws))
  
  nbk <- design_table(design=x@Design, x="treatment",y="block")
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
