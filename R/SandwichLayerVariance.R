#' @include Design.R SandwichLayer.R
NULL

#' @title Variance/Covariance for \code{teeMod} objects
#'
#' @description Compute robust sandwich variance estimates with optional covariance
#'   adjustment
#'
#' @details Supported \code{type} include:
#'
#' - \code{"MB0"}, \code{"HC0"}, and \code{"CR0"} for model-based HC0 standard errors
#' - \code{"MB1"}, \code{"HC1"}, and \code{"CR1"} for model-based standard errors
#' with HC1 corrections based on the direct adjustment estimate i.e.,
#' \eqn{n/(n - 2)} for \code{"MB1"} and \code{"HC1"}, and for \code{"CR1"},
#' \eqn{g\cdot(n-1)/((g-1)\cdot(n-2))}, where \eqn{g} is the number of clusters
#' in the direct adjustment sample.
#'
#' To create your own \code{type}, simply define a function \code{.vcov_XXX}.
#' \code{type = "XXX"} will now use your method. Your method should return a
#' matrix of appropriate dimension, with \code{attribute} \code{type = "XXX"}.
#'
#' @param x a fitted \code{teeMod} model
#' @param type a string indicating the desired variance estimator. See Details
#'   for supported variance estimators
#' @param cluster a string or character vector of column names indicating
#'   columns to cluster standard errors by. With prior covariance adjustment,
#'   columns must appear in both the covariance adjustment and direct adjustment
#'   samples. Default is NULL, which uses unit of assignment columns in the
#'   \code{Design} slot of the \code{teeMod} model.
#' @param ... arguments to be passed to the internal variance estimation
#'   function.
#' @return A \eqn{2\times 2} matrix corresponding to an intercept and the
#'   treatment variable in the direct adjustment model
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

#' @keywords internal
.vcov_CR0 <- function(x, ...) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
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

#' @keywords internal
.vcov_HC0 <- function(x, ...) {
  out <- .vcov_CR0(x, ...)
  attr(out, "type") <- "HC0"
  return(out)
}

#' @keywords internal
.vcov_MB0 <- function(x, ...) {
  out <- .vcov_CR0(x, ...)
  attr(out, "type") <- "MB0"
  return(out)
}

#' @keywords internal
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

#' @keywords internal
.vcov_HC1 <- function(x, ...) {
  out <- .vcov_CR1(x, ...)
  attr(out, "type") <- "HC1"
  return(out)
}

#' @keywords internal
.vcov_MB1 <- function(x, ...) {
  out <- .vcov_CR1(x, ...)
  attr(out, "type") <- "MB1"
  return(out)
}

#' @title (Internal) Replace standard errors for moderator effect estimates with insufficient
#' degrees of freedom with \code{NA}
#' @param vmat output of \code{.vcov_XXX()} called with input to \code{model}
#'  argument below as the first argument
#' @param model a fitted \code{teeMod} model
#' @param cluster a character or factor vector denoting clusters the units of
#'  observation used to fit \code{model} belong to. Most conveniently, the
#'  output of \code{.make_uoa_ids()}
#' @param model_data dataframe or name of dataframe used to fit \code{model}
#' @param envir environment to get \code{model_data} from if \code{model_data} has class \code{name}
#' @return A variance-covariance matrix with NA's for estimates lacking sufficient
#' degrees of freedom
#' @keywords internal
.check_df_moderator_estimates <- function(vmat, model, cluster, model_data = quote(data),
                                          envir = environment(formula(model))) {
  if (!inherits(model, "teeMod")) {
    stop("`model` must be a teeMod object")
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

#' @details \code{.get_b12()}: \eqn{B_{12}} is the covariance of the unit of assignment- or cluster-level
#'   estimating equations for the covariance adjustment and direct adjustment models.
#'   It has a row for each term in the former and a column for each term in the latter.
#'   Any unit of assignment or cluster that does not appear in both model fits
#'   makes no contribution to this matrix. If there is no overlap between the two
#'   datasets, this function will return a matrix of zeros of appropriate dimension.
#' @param x a fitted \code{teeMod} model
#' @param ... arguments to pass to internal functions, such as \code{cluster}
#' @return \code{.get_b12()}: A \eqn{p\times 2} matrix
#' @keywords internal
#' @noRd
.get_b12 <- function(x, ...) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }

  if (x@.S3Class[1] != "lm") {
    stop("x must be an `lm` object")
  }

  sl <- x$model$`(offset)`
  if (!is(sl, "SandwichLayer")) {
    stop(paste("teeMod model must have an offset of class `SandwichLayer`",
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
                 "in the teeMod object's Design:",
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
               "exist in both the direct adjustment and covariance adjustment samples"))
  }

  Q_uoas <- stats::expand.model.frame(x, cluster_cols, na.expand = TRUE)[cluster_cols]
  Q_uoas <- apply(Q_uoas, 1, function(...) paste(..., collapse = "_"))
  C_uoas <- apply(keys[, cluster_cols, drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  C_uoas_in_Q <- C_uoas %in% unique(Q_uoas)

  message(paste(sum(C_uoas_in_Q),
                "rows in the covariance adjustment data",
                "joined to the direct adjustment data\n"))

  # Check number of overlapping clusters n_QC; if n_QC <= 1, return 0 (but
  # similarly to .get_b11(), throw a warning when only one cluster overlaps)
  uoas_overlap <- length(unique(C_uoas[C_uoas_in_Q]))
  if (uoas_overlap <= 1) {
    if (uoas_overlap == 1) {
      warning(paste("Covariance matrix between covariance adjustment and direct adjustment",
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

#' @title (Internal) Estimate components of the sandwich covariance matrix
#' returned by \code{vcovDA()}
#' @details \code{.get_a22_inverse()}/\code{.get_tilde_a22_inverse()}: \eqn{A_{22}^{-1}} is the "bread" of the
#' sandwich covariance matrix returned by \code{vcovDA()} whether one has fit
#' a prior covariance adjustment model or not.
#' @param x a fitted \code{teeMod} object
#' @param ... arguments passed to \code{bread} method
#' @return \code{.get_a22_inverse()}/\code{.get_tilde_a22_inverse()}: A
#' \eqn{2\times 2} matrix corresponding to an intercept and the treatment
#' variable in the direct adjustment model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a22_inverse <- function(x, ...) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }

  out <- utils::getS3method("bread", "lm")(x)

  return(out)
}

#' @param x a fitted \code{teeMod} model
#' @param ... arguments to pass to internal functions
#' @details \code{.get_b22()}: \eqn{B_{22}} is the covariance of the unit of
#'   assignment- or cluster-level estimating equations for the direct adjustment
#'   model. In the absence of covariance adjustment, this is the meat of the
#'   sandwich covariance matrix for the direct adjustment model.
#' @return \code{.get_b22()}: A \eqn{2\times 2} matrix corresponding to an
#' intercept and the treatment variable in the direct adjustment model
#' @keywords internal
#' @noRd
.get_b22 <- function(x, ...) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }

  nq <- nrow(sandwich::estfun(x))

  # Create cluster ID matrix depending on cluster argument (or its absence)
  dots <- list(...)
  if (is.null(dots$cluster)) {
    uoas <- .expand.model.frame_teeMod(x,
                         var_names(x@Design, "u"))[, var_names(x@Design, "u"),
                                                   drop = FALSE]
  } else if (inherits(dots$cluster, "character")) {
    uoas <- tryCatch(
      .expand.model.frame_teeMod(x, dots$cluster)[,
                                                        dots$cluster,
                                                        drop = FALSE],
      error = function(e) {
        data <- eval(x$call$data,
                     envir = environment(formula(x)))
        stop(paste("The columns",
                   paste(setdiff(dots$cluster, colnames(data)), collapse = ", "),
                   "are missing from the direct adjustment data"),
             call. = FALSE)
      })
  } else {
    stop(paste("If overriding `cluster` argument for meat matrix calculations,",
               "must provide a character vector specifying column names in the",
               "direct adjustment data"))
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

#' @details \code{.get_a11_inverse()}: \eqn{A_{11}^{-1}} is the "bread" of
#'   the sandwich covariance matrix for the covariance adjustment model. This
#'   matrix contributes to the meat matrix of the direct adjustment sandwich
#'   covariance matrix.
#' @param x a fitted \code{teeMod} model
#' @return \code{.get_a11_inverse()}: A \eqn{p\times p} matrix where \eqn{p} is
#'   the dimension of the covariance adjustment model, including an intercept
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a11_inverse <- function(x) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("teeMod model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  out <- sandwich::bread(sl@fitted_covariance_model)

  return(out)
}

#' @param x a fitted \code{teeMod} model
#' @param ... arguments to pass to internal functions
#' @details \code{.get_b11()}: \eqn{B_{11}} is the covariance of the unit of
#'   assignment- or cluster-level estimating equations for the covariance
#'   adjustment model. This matrix contributes to the meat matrix of the direct
#'   adjustment sandwich covariance matrix.
#' @return \code{.get_b11()}: A \eqn{p\times p} matrix where \eqn{p} is
#'   the dimension of the covariance adjustment model, including an intercept
#' @keywords internal
#' @noRd
.get_b11 <- function(x, ...) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("teeMod model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }

  cmod <- sl@fitted_covariance_model
  nc <- sum(summary(cmod)$df[1L:2L])
  if (nc==0) nc <- length(cmod$fitted.values)
  if (nc==0) stop("can't determine extent of covariance model fitting data")

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

#' @details \code{.get_a21()}/\code{.get_tilde_a21()}: \eqn{A_{21}} is the
#'   gradient of the estimating equations for the direct adjustment model taken
#'   with respect to the covariance adjustment model parameters. This matrix is
#'   the crossproduct of the prediction gradient for the units of observation in
#'   \eqn{\mathcal{Q}} and the model matrix of the direct adjustment model.
#' @param x a fitted \code{teeMod} model
#' @return \code{.get_a21()}/\code{.get_tilde_a21()}: A \eqn{2\times p} matrix
#'   where the number of rows are given by the intercept and the treatment
#'   variable in the direct adjustment model, and the number of columns are given
#'   by the dimension of the covariance adjustment model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a21 <- function(x) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("teeMod model must have an offset of class `SandwichLayer`",
               "to propagate covariance adjustment model error"))
  }

  # Get contribution to the estimating equation from the direct adjustment model
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
