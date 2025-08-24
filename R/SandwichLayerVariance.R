#' @include StudySpecification.R SandwichLayer.R block_center_residuals.R
NULL

#' @title Variance/Covariance for \code{teeMod} objects
#'
#' @description Compute robust sandwich variance estimates with optional
#'   covariance adjustment
#'
#' @details Variance estimates will be clustered on the basis of the columns
#'  provided to \code{cluster} (or obtained by the default behavior). As a result,
#'  providing \code{"HCx"} or \code{"CRx"} to \code{type} will produce the same
#'  variance estimate given that \code{cluster} remains the same.
#'
#' With prior covariance adjustment, unless the \code{data} argument of the covariance
#'  model fit is the same as the \code{data} argument for fitting \code{x} and the
#'  \code{StudySpecification} of \code{x} has been created with a formula of the
#'  form `trt_col ~ 1`, the column(s) provided to \code{cluster} must appear in the
#'  dataframes in both \code{data} arguments, even if the clustering structure does
#'  not exist, per se, in the covariance adjustment sample. For instance, in a
#'  finely stratified randomized trial, one might desire standard errors clustered
#'  at the block level, but the covariance adjustment model may include auxiliary
#'  units that did not participate in the trial. In this case, in the \code{data}
#'  argument of the fitted covariance model, the column(s) passed to \code{cluster}
#'  should have the block ID's for rows overlapping with the \code{data} argument
#'  used for fitting \code{x}, and NA's for any auxiliary units. \code{vcov_tee()}
#'  will treat each row with an NA as its own cluster.
#'  
#' For ITT effect estimates without covariance adjustment, \code{type}
#'  corresponds to the variance estimate desired. Supported options include:
#'
#' - \code{"MB0"}, \code{"HC0"}, and \code{"CR0"} for model-based HC/CR0 standard errors
#' - \code{"MB1"}, \code{"HC1"}, and \code{"CR1"} for model-based HC/CR1 standard errors
#' (for \code{"MB1"} and \code{"HC1"}, this is \eqn{n/(n - 2)}, and for \code{"CR1"},
#' this is \eqn{g\cdot(n-1)/((g-1)\cdot(n-2))}, where \eqn{g} is the number of
#' clusters in the sample used for fitting \code{x})
#' - \code{"MB2"}, \code{"HC2"}, and \code{"CR2"} for model-based HC/CR2 standard errors
#' - \code{"DB0"} for design-based HC0 variance estimates
#' 
#' The \code{type} argument does not correspond to existing variance estimators in
#'   the literature in the case of prior covariance adjustment. It specifies the
#'   bias correction to the residuals of \code{x}, but the residuals of the covariance
#'  model are corrected separately based on the \code{cov_adj_rcorrect} argument.
#'  The \code{cov_adj_rcorrect} argument takes the same options as \code{type}
#'  except \code{"DB0"}. When the covariance model includes rows in the treatment
#'  condition for fitting, the residuals of \code{x} are further corrected by
#'  having the values of \code{offset} replaced by predictions that use coefficient
#'  estimates that leave out rows in the same cluster (as defined by the
#'  \code{cluster} argument).
#'
#' @param x a fitted \code{teeMod} model
#' @param type a string indicating the desired bias correction for the residuals
#'   of \code{x}. Default makes no bias correction. See Details for supported types
#' @param cluster a vector indicating the columns that define clusters. The default
#'   is the unit of assignment columns in the \code{StudySpecification} stored
#'   in \code{x}. These columns should appear in the dataframe used for fitting
#'   \code{x} as well as the dataframe passed to the covariance model fit in the case of
#'   prior covariance adjustment. See Details
#' @param ... arguments to be passed to the internal variance estimation
#'   function, such as \code{cov_adj_rcorrect} and \code{loco_residuals}. If
#'   \code{x} has a \code{SandwichLayer} object in its offset, The former
#'   specifies the bias correction to the residuals of the covariance model, and
#'   the latter indicates whether the offset should be replaced with predictions
#'   from leave-one-cluster-out fits of the covariance adjustment model.
#'   See Details
#' @return A variance-covariance matrix with row and column entries for the estimated
#'   coefficients in \code{x}, the marginal mean outcome in the control condition,
#'   the marginal mean \code{offset} in the control condition (if an \code{offset}
#'   is provided), and if a moderator variable is specified in the formula for \code{x},
#'   the mean interaction in the control condition of the outcome and \code{offset}
#'   with the moderator variable
#' @export
#' @rdname var_estimators
vcov_tee <- function(x, type = NULL, cluster = NULL, ...) {
  if (is.null(type)) type <- "HC0"

  vcov_type <- substr(type, 1, 2)
  if (type != "DB0" & !(vcov_type %in% c("HC", "CR", "MB"))) {
    stop("covariance function not defined.\n")
  }
  var_func <- get(paste0(".vcov_", ifelse(vcov_type == "DB", type, vcov_type)))
  
  args <- list(...)
  args$x <- x
  
  if (vcov_type == "DB") {
    args$vcov_type <- vcov_type
  } else {
    args$type <- type
  }
  
  if (is.null(args$cov_adj_rcorrect) & inherits(x, "teeMod") &
      inherits(x$model$`(offset)`, "SandwichLayer")) args$cov_adj_rcorrect <- "HC0"

  if (is.null(args$by)) args$by <- cluster # if cov_adj() was not fit with a "by" argument, this is passed to .order_samples() to order rows of estfun() output
  args$cluster_cols <- cluster
  args$cluster <- .make_uoa_ids(x, vcov_type = vcov_type, cluster, ...) # passed on to meatCL to aggregate SE's at the level given by `cluster`

  est <- do.call(var_func, args)

  return(est)
}

#' @keywords internal
.vcov_CR <- function(x, ...) {
  if (!inherits(x, c("teeMod", "mmm"))) {
    stop("x must be a teeMod or mmm object")
  }

  args <- list(...)
  args$x <- x
  args$itt_rcorrect <- args$type
  args$type <- "HC0" # this ensures sandwich::meatCL returns meat w/o correcting after we've corrected
  n <- length(args$cluster)
  args$cls <- args$cluster # pass this down to estfun.teeMod and its internal calls to try and speed things up

  bread. <- sandwich::bread(x, ...)
  meat. <- do.call(sandwich::meatCL, args)
  vmat <- (1 / n) * bread. %*% meat. %*% t(bread.)

  # NA any invalid estimates due to degrees of freedom checks
  if (inherits(x, "teeMod")) {
    vmat <- .check_df_moderator_estimates(vmat, x, args$cluster_cols)
  } else {
    start_ix <- 0
    for (mod_ix in seq_along(x)) {
      mod <- x[[mod_ix]]
      vmat_ix <- start_ix + seq_along(mod$coefficients)
      vmat[vmat_ix, vmat_ix] <- .check_df_moderator_estimates(vmat, mod, args$cluster_cols)[vmat_ix, vmat_ix]
      start_ix <- start_ix + length(mod$coefficients)
    }
    vmat[apply(is.na(vmat), 1, any),] <- NA_real_
    vmat[,apply(is.na(vmat), 2, any)] <- NA_real_
  }

  attr(vmat, "type") <- match.call()$type
  attr(vmat, "cov_adj_rcorrect") <- args$cov_adj_rcorrect
  return(vmat)
}

#' @keywords internal
.vcov_HC <- function(x, ...) {
  out <- .vcov_CR(x, ...)
  return(out)
}

#' @keywords internal
.vcov_MB <- function(x, ...) {
  out <- .vcov_CR(x, ...)
  return(out)
}

#' @title Get the degrees of freedom of a sandwich variance estimate associated with a teeMod fit
#' @keywords internal
.get_dof <- function(x, vcov_type, ell, cluster = NULL, cls = NULL, ...) {
  if (!inherits(x, "teeMod")) stop("x must be a teeMod object")
  if (is.null(cluster)) cluster <- var_names(x@StudySpecification, "u")
  if (is.null(cls)) cls <- .make_uoa_ids(x, substr(vcov_type, 1, 2), cluster)
  
  if (!(vcov_type %in% paste0(c("MB", "HC", "CR"), 2))) {
    # if no HC2 correction, use dof = G - 1 (see pg. 331 of Cameron and Miller,
    # 2015). this is correct for both absorb=FALSE and absorb=TRUE
    dof <- length(unique(cls)) - 1
  } else {
    # use dof from IK for HC2 corrections
    k <- ncol(x$qr$qr)
    if (length(ell) == 1) {
      ell <- replace(numeric(k), ell, 1)
    } else if (length(ell) != k) {
      stop("ell must have the same length as the rank of the fitted model")
    } else if (length(setdiff(unique(ell), c(0, 1))) != 0) {
      stop("ell must have only zeros or ones")
    }
    dof <- .compute_IK_dof(x, ell, cluster, length(unique(x$model[,1])) == 2)
  }
  
  return(dof)
}

#' @title Compute the degrees of freedom of a sandwich standard error with HC2 correction
#' @importFrom utils combn
#' @importFrom stats na.action
#' @keywords internal
#' @references Guido W. Imbens and Michael Kolesár. "Robust Standard Errors in
#' Small Samples: Some Practical Advice". In: *The Review of Economics and
#' Statistics* 98.4 (Oct. 2016), pp. 701-712.
.compute_IK_dof <- function(tm,
                            ell,
                            cluster = NULL,
                            bin_y = FALSE,
                            exclude = na.action(tm),
                            tol = 1e-09) {
  if (is.null(cluster)) cluster <- var_names(tm@StudySpecification, "u")
  cls <- .sanitize_Q_ids(tm, cluster)$cluster
  cls <- cls[setdiff(seq_along(cls), exclude)]
  lmitt_data <- get("data", environment(formula(tm)))
  if (!is.null(sbst <- tm@lmitt_call$subset)) lmitt_data <- subset(lmitt_data, eval(sbst, envir = lmitt_data))
  wc <- tm@lmitt_call$weights
  if (inherits(wc, "call")) {
    if (exists(deparse1(wc[[1L]]))) {
      if (length(wc) == 1) {
        wc$specification <- get("specification", environment(formula(tm)))
      } else if (!inherits(eval(wc[[2L]], environment(formula(tm))), "StudySpecification")) {
        wc$specification <- get("specification", environment(formula(tm)))
      }
      if (is.null(wc$data)) wc$data <- lmitt_data
    }
  }
  trts <- if (is.null(tm@lmitt_call$dichotomy) &
              (is.null(wc <- eval(wc, envir = environment(formula(tm)))) |
               !inherits(wc, "WeightedStudySpecification"))) {
    treatment(tm@StudySpecification, newdata = lmitt_data)[,1]
  } else {
    if (length(dich <- wc@dichotomy) == 0) {
      treatment(tm@StudySpecification, newdata = lmitt_data)[,1]
    } else {
      eval(dich, envir = lmitt_data)
    }
  }
  
  ## estimate rho
  r <- stats::residuals(tm, type = "response")
  if (length(cls) == length(unique(cls)) | bin_y) {
    # for nonclustered data rho = 0. for binary data, use a working model of independence
    # so that rho = 0
    rho <- 0
  } else {
    # don't need to subtract off squared residuals because combn only returns combos
    # of distinct indices
    rho <- sum(
      tapply(r, cls, function(rs) {
        combs <- combn(length(rs), 2)
        sum(rs[combs[1,]] * rs[combs[2,]])
      })
    ) / (sum(tapply(rep(1, length(cls)), cls, sum)^2) - length(cls))
    rho <- max(rho, 0) # IK code doesn't allow negative rho
  }
  
  ## estimate sigma
  # working correlation matrix is heteroskedastic for binary data
  sig2 <- if (bin_y) tm$fitted.values * (1-tm$fitted.values) else mean(r^2)
  
  # function for getting the inverse symmetric square root 
  iss <- function(s, cls, trts = NULL) {
    if (length(cls) == length(unique(cls)) | bin_y) {
      # no clustering or independent working correlation matrix
      return(matrix(1 / sqrt(1 - stats::hatvalues(tm)[cls == s]), 1, 1))
    } else {
      return(cluster_iss(tm, cluster_unit = s, cluster_ids = cls, assigned_trt = trts))
    }
  }
  
  # this is IK code with: 1) our fast CR2 correction implemented,
  # 2) calls to the collapse package replaced with calls to rowsum(), and 3)
  # removing the intercept column when there's absorption because the corrections
  # for those cases in cluster_iss() are based on having no intercept column
  piv <- if (tm@absorbed_intercepts & !bin_y &
             length(cls) != length(unique(cls))) seq(2,tm$rank) else seq_len(tm$rank)
  Q <- qr.Q(tm$qr)[, piv,drop=FALSE]
  R <- qr.R(tm$qr)[piv, piv,drop=FALSE]
  if (length(cls) == length(unique(cls)) | bin_y) {
    AQ <- Q / sqrt(1-stats::hatvalues(tm))
  } else {
    AQ <- Reduce(
      rbind,
      lapply(
        unique(cls),
        function(s, cls, trts = NULL) {
          Qs <- Q[cls == s,,drop=FALSE]
          iss(s, cls, trts) %*% Qs
        },
        cls, trts
      )
    )
  }
  
  # if there are insufficient degrees of freedom, AQ will be numerically 0.
  # this is also the case in the dfadjust package, but that package does not
  # proactively address numerical zeros, instead allowing the division to wash
  # out and return 1 DOF. we address that here by rounding off AQ, and if all
  # entries are 0, we actively return 1 dof
  if (isTRUE(
    all.equal(AQ, matrix(0, nrow = nrow(AQ), ncol = ncol(AQ)), tolerance = tol)
  )) {
    return(1)
  }
  a <- drop(AQ %*% backsolve(R, ell, transpose = TRUE))
  as <- rowsum((a*sqrt(sig2))^2, cls)[,1] # accommodate heteroskedasticity
  B <- rowsum(a * sqrt(sig2) * Q, cls)
  D <- rowsum(a, cls)[,1]
  Fm <- rowsum(Q, cls)
  GG <- (diag(as, nrow = length(as)) - tcrossprod(B)) + rho * 
    tcrossprod(diag(D, length(D)) - tcrossprod(B, Fm))
  sum(diag(GG))^2/sum(GG^2)
}

#' @title Use properties of idempotent matrices to cheaply compute inverse symmetric
#' square roots of cluster-specific subsets of projection matrices
#' @param exclude index of units to exclude from computing the correction; for
#' example, if they're NA's
#' @importFrom stats model.frame na.action
#' @keywords internal
cluster_iss <- function(tm,
                        cluster_unit,
                        cluster_ids = NULL,
                        cluster_var = NULL,
                        exclude = na.action(tm),
                        tol = 1e-09,
                        ...) {
  if (!inherits(tm, "teeMod")) stop("Must provide a teeMod object")
  dots <- list(...)
  if (is.null(cluster_ids)) {
    if (is.null(dots$cluster)) cluster <- var_names(tm@StudySpecification, "u")
    cluster_ids <- .sanitize_Q_ids(tm, cluster)$cluster
  }
  ix <- setdiff(which(cluster_ids == cluster_unit), exclude)
  
  lmitt_data <- get("data", environment(formula(tm)))
  if (!is.null(sbst <- tm@lmitt_call$subset)) lmitt_data <- subset(lmitt_data, eval(sbst, envir = lmitt_data))
  if (is.null(trts <- dots$assigned_trt)) {
    wc <- tm@lmitt_call$weights
    if (inherits(wc, "call")) {
      if (exists(deparse1(wc[[1L]]))) {
        if (length(wc) == 1) {
          wc$specification <- get("specification", environment(formula(tm)))
        } else if (!inherits(eval(wc[[2L]], environment(formula(tm))), "StudySpecification")) {
          wc$specification <- get("specification", environment(formula(tm)))
        }
        if (is.null(wc$data)) wc$data <- lmitt_data
      }
    }
    trts <- if (is.null(tm@lmitt_call$dichotomy) &
                (is.null(wc <- eval(wc, envir = environment(formula(tm)))) |
                 !inherits(wc, "WeightedStudySpecification"))) {
      treatment(tm@StudySpecification, newdata = lmitt_data)[,1]
    } else {
      if (length(dich <- wc@dichotomy) == 0) {
        treatment(tm@StudySpecification, newdata = lmitt_data)[,1]
      } else {
        eval(dich, envir = lmitt_data)
      }
    }
  }

  trtg <- trts[ix]
  tt <- stats::delete.response(stats::terms(tm))
  mf <- call("model.frame",
             tt,
             lmitt_data,
             subset = sbst,
             na.action = na.pass,
             xlev = tm$xlevels)
  mf <- eval(mf)
  K <- tm$qr$rank
  piv <- tm$qr$pivot[seq_len(K)]
  A <- stats::model.matrix(tt, mf, contrasts.arg = tm$contrasts)[,piv,drop=FALSE]
  Ag <- A[ix,,drop=FALSE]
  inv <- chol2inv(tm$qr$qr[piv,piv])
  if (is.null(wts <- stats::weights(tm))) wts <- rep(1, nrow(A))
  wg <- wts[ix]

  cg <- NULL
  Mgg <- NULL
  if (length(unique(trtg)) == 1) {
    # first, check for a moderator variable
    if (length(tm@moderator) > 0) {
      xvar <- lmitt_data[[tm@moderator]]
      # make sure the moderator variable is invariant within the cluster
      if (length(x <- unique(xvar[ix])) == 1) {
        if (tm@absorbed_intercepts) {
          # with an invariant moderator variable, the 1st row of the cluster model
          # matrix is the same as the rest
          cg <- sum(wg) * drop(crossprod(Ag[1,], inv) %*% Ag[1,])
          Mgg <- tcrossprod(sqrt(wg)) / sum(wg)
        } else {
          if (inherits(xvar, c("factor", "character"))) {
            # control clusters and treated clusters have different values of cg
            if (unique(trtg) == 0) {
              xcols <- which(colnames(Ag) %in% paste0(tm@moderator, x))
            } else {
              xcols <- which(colnames(Ag) %in% paste0(
                c(paste0(var_names(tm@StudySpecification, "t"), "._"), ""),
                tm@moderator, x)
              )
            }
            cg <- sum(inv[c(1, xcols), c(1, xcols)]) * sum(wg)
            Mgg <- tcrossprod(sqrt(wg)) / sum(wg)
          } else {
            if (unique(trtg) == 0) {
              cg <- (inv[1,1] + Ag[1,3] * inv[1,3] * 2 + Ag[1,3]^2 * inv[3,3]) * sum(wg)
              Mgg <- tcrossprod(sqrt(wg)) / sum(wg) 
            } else {
              ul <- sum(inv[1:2, 1:2])
              br <- sum(inv[3:4, 3:4])
              cg <- (ul + br * Ag[1,3]^2 + (sum(inv) - ul - br) * Ag[1,3]) * sum(wg)
              Mgg <- tcrossprod(sqrt(wg)) / sum(wg) 
            }
          }
        }
      }
    } else {
      if (tm@absorbed_intercepts) {
        cg <- sum(wg * Ag[,2]^2)
        # the 1st row of the cluster model matrix is the same as the rest
        cg <- sum(wg) * (1 / sum(wts) + Ag[1,2]^2 / sum(wts * A[,2]^2))
        Mgg <- tcrossprod(sqrt(wg)) / sum(wg)
      } else {
        # control clusters and treated clusters have different values of cg
        if (unique(trtg) == 0) {
          cg <- sum(wg) / sum(wts[setdiff(which(A[,2] == 0), exclude)])
          Mgg <- tcrossprod(sqrt(wg)) / sum(wg)
        } else {
          cg <- sum(wg) / sum(wts[setdiff(which(A[,2] == 1), exclude)])
          Mgg <- tcrossprod(sqrt(wg)) / sum(wg)
        }
      }
    }
  } else {
    if (length(tm@moderator) == 0) {
      # when clusters contain both treated and control units, a sufficient
      # condition for simple construction of cg and Mgg is a constant ratio of weighted
      # size of the treatment group to the weighted size of the control group
      # across blocks
      rts <- tapply(data.frame(w = wts, t = trts),
                    cluster_ids,
                    function(df) {
                      wga <- rowsum(df[,1], df[,2])
                      wga[2,]/wga[1,]
                    })
      if (all(rts == mean(rts))) {
        Mgg <- matrix(0, nrow = length(trtg), ncol = length(trtg))
        Mgg[trtg == 0,trtg == 0] <- tcrossprod(sqrt(wg[trtg == 0])) / sum(wg[trtg == 0])
        Mgg[trtg == 1,trtg == 1] <- tcrossprod(sqrt(wg[trtg == 1])) / sum(wg[trtg == 1])
        # under this condition, cg's below are the same whether we use treated or control units
        if (tm@absorbed_intercepts) {
          cg <- sum(wg[trtg == 0]) * (1 / sum(wts) + Ag[trtg == 0,2][1]^2 / sum(wts * A[,2]^2))
        } else {
          cg <- sum(wg[trtg == 0]) / sum(wts[trts == 0])
        }
      }
    }
  }
  
  # cg == 1 indicates our speedup can't be accommodated, so fall back to
  # original CR2 computation (use all.equal because, numerically, cg may just be
  # close to 1--and use within isTRUE per the all.equal documentation)
  if ((!is.null(cg) && isTRUE(all.equal(cg, 1))) | is.null(cg) | is.null(Mgg)) {
    Pgg <- Ag %*% inv %*% t(Ag * wg)
    eg <- eigen(diag(nrow = length(ix)) - Pgg)
    return(eg$vec %*% diag((eg$val >= tol) * 1/sqrt(pmax(eg$val, tol))) %*%
             solve(eg$vec))
  } else {
    return(diag(nrow = length(ix)) + (-1 + sqrt(1 + cg / (1-cg))) * Mgg)
  }
}

#' @title (Internal) Replace standard errors for moderator effect estimates with
#'   insufficient degrees of freedom with \code{NA}
#' @param vmat output of \code{.vcov_XXX()} called with input to \code{model}
#'   argument below as the first argument
#' @param model a fitted \code{teeMod} model
#' @param cluster_cols a character vector indicating the column(s) defining
#'   cluster ID's
#' @param model_data dataframe or name of dataframe used to fit \code{model}
#' @param envir environment to get \code{model_data} from if \code{model_data}
#'   has class \code{name}
#' @return A variance-covariance matrix with NA's for estimates lacking
#'   sufficient degrees of freedom
#' @keywords internal
.check_df_moderator_estimates <- function(vmat, model, cluster_cols, model_data = quote(data),
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
    stop("`model_data` must be a dataframe or a quoted object name")
  }
  
  if (!inherits(model_data, "data.frame")) {
    stop(paste("Could not find argument passed to `model_data` in the given `envir`"))
  }
  
  # set these attributes so model.matrix() includes rows with NA's, which will match
  # the length of `cluster` generated below
  attr(model_data, "terms") <- NULL
  attr(model_data, "na.action") <- "na.pass"

  # For categorical moderators, count the clusters contributing to estimation
  # for each level of the moderator variable; for continuous moderators, just
  # count the number of clusters. The moderator variable (or any level of the
  # moderator variable) must have at least three clusters contributing to
  # estimation for valid SE estimation.
  mod_vars <- model.matrix(as.formula(paste0("~-1+", model@moderator)), model_data)
  cluster <- .sanitize_Q_ids(model, cluster_cols)$cluster
  if (ncol(mod_vars) > 1) {
    # Since model_vars and cluster may include rows with NA's that weren't used
    # to fit the model, we need to create in_model_fit to determine rows that
    # were. The approach below uses the row.names in the na.action to indicate
    # which were not
    # in_model_fit <- as.numeric(
    #   apply(mapply(function(c) is.na(c),
    #                eval(attr(terms(as.formula(model$call$formula)), "variables"), env = model_data)),
    #         1,
    #         function(r) sum(r) == 0)
    # )
    in_model_fit <- rep(1, length(cluster))
    in_model_fit[model$na.action] <- 0
    mod_counts <- rowsum(mod_vars * in_model_fit, cluster, na.rm = TRUE)
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

#' @details \code{.get_b12()}: \eqn{B_{12}} is the covariance of the unit of
#'   assignment- or cluster-level estimating equations for the covariance
#'   adjustment and direct adjustment models. It has a row for each term in the
#'   former and a column for each term in the latter. Any unit of assignment or
#'   cluster that does not appear in both model fits makes no contribution to
#'   this matrix. If there is no overlap between the two datasets, this function
#'   will return a matrix of zeros of appropriate dimension.
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
    cluster_cols <- var_names(x@StudySpecification, "u")
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

    # Check to see if provided column names overlap with specification
    # missing_spec_cols <- setdiff(cov_adj_cluster_cols, colnames(x@StudySpecification@structure))
    missing_spec_cols <- setdiff(cluster_cols, colnames(x@StudySpecification@structure))
    if (length(missing_spec_cols) > 0) {
      stop(paste("The following columns in the `cluster` argument cannot be found",
                 "in the teeMod object's StudySpecification:",
                 paste(missing_spec_cols, collapse = ", ")))
    }

    # Re-create keys dataframe with the new clustering columns
    keys <- as.data.frame(
      sapply(cluster_cols, function(col) {
        unique(x@StudySpecification@structure[[col]])[
          match(wide_frame[[col]], unique(x@StudySpecification@structure[[col]]), incomparables = NA)
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
#'   returned by \code{vcov_tee()}
#' @details \code{.get_a22_inverse()}/\code{.get_tilde_a22_inverse()}:
#'   \eqn{A_{22}^{-1}} is the "bread" of the sandwich covariance matrix returned
#'   by \code{vcov_tee()} whether one has fit a prior covariance adjustment
#'   model or not.
#' @param x a fitted \code{teeMod} object
#' @param ... arguments passed to \code{bread} method
#' @return \code{.get_a22_inverse()}/\code{.get_tilde_a22_inverse()}: A
#'   \eqn{2\times 2} matrix corresponding to an intercept and the treatment
#'   variable in the direct adjustment model
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a22_inverse <- function(x, ...) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }
  mc <- match.call()
  vcov_type <- match.arg(mc$vcov_type, c("MB", "HC", "CR", "DB")) # this sets the default to model-based
  
  teeMod_bread <- utils::getS3method("bread", "lm")(x)
  if (vcov_type == "DB") {
    return(teeMod_bread)
  }

  cm_mod <- x@ctrl_means_model
  ctrl_means_bread <- bread(cm_mod)
  if (!inherits(cm_mod, "mlm")) {
    dimnames(ctrl_means_bread) <- lapply(
      dimnames(ctrl_means_bread),
      function(nms, x) paste(formula(x)[[2L]], nms, sep = ":"),
      x = x)
  }

  out <- matrix(0, nrow = nrow(teeMod_bread) + nrow(ctrl_means_bread),
                ncol = ncol(teeMod_bread) + ncol(ctrl_means_bread),
                dimnames = list(c(rownames(teeMod_bread), rownames(ctrl_means_bread)),
                                  c(colnames(teeMod_bread), colnames(ctrl_means_bread))))
  out[1:nrow(teeMod_bread), 1:ncol(teeMod_bread)] <- teeMod_bread
  out[(nrow(teeMod_bread)+1):nrow(out), (ncol(teeMod_bread)+1):ncol(out)] <- ctrl_means_bread

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
                         var_names(x@StudySpecification, "u"))[, var_names(x@StudySpecification, "u"),
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

#' @details \code{.get_a11_inverse()}: \eqn{A_{11}^{-1}} is the "bread" of the
#'   sandwich covariance matrix for the covariance adjustment model. This matrix
#'   contributes to the meat matrix of the direct adjustment sandwich covariance
#'   matrix.
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

  # Replace NA's for rows not in the experimental specification with a unique cluster ID
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
#'   variable in the direct adjustment model, and the number of columns are
#'   given by the dimension of the covariance adjustment model
#' @importFrom stats weights formula
#' @keywords internal
#' @rdname sandwich_elements_calc
.get_a21 <- function(x, ...) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }

  sl <- x$model$`(offset)`
  if (!inherits(sl, "SandwichLayer")) {
    stop(paste("teeMod model must have an offset of class `SandwichLayer`",
               "to propagate covariance adjustment model error"))
  }
  
  # Get contribution to the estimating equation from the direct adjustment model
  damod_mm <- stats::model.matrix(
    formula(x), stats::model.frame(x, na.action = na.pass))
  msk <- (apply(!is.na(sl@prediction_gradient), 1, all) &
            apply(!is.na(damod_mm), 1, all))
  if (!is.null(x$na.action)) class(x$na.action) <- "exclude"
  w <- if (is.null(w <- stats::weights(x))) numeric(length(msk)) + 1 else replace(w, is.na(w), 0)
  wZ <- damod_mm[, x$qr$pivot[1L:x$rank], drop = FALSE] * w
  
  mc <- match.call()
  vcov_type <- match.arg(mc$vcov_type, c("MB", "HC", "CR", "DB")) # this sets the default to model-based
  if (vcov_type != "DB") {
    # get gradient for ctrl means model
    cm <- x@ctrl_means_model
    cm_mm <- stats::model.matrix(formula(cm), stats::model.frame(cm, na.action = na.pass))
    cm_grad <- matrix(0, nrow = nrow(cm_mm), ncol = ncol(cm_mm),
                      dimnames = list(NULL, paste(formula(x)[[2]], colnames(cm_mm), sep = ":")))
    colnames(cm_mm) <- paste("cov_adj", colnames(cm_mm), sep = ":")
    if (inherits(cm, "mlm")) cm_grad <- cbind(cm_grad, cm_mm)
    cm_wts <- replace(cm_wts <- stats::weights(cm), is.na(cm_wts), 0)
    wZ <- cbind(wZ, -cm_grad * cm_wts)
  }

  out <- crossprod(wZ[msk,], sl@prediction_gradient[msk, , drop = FALSE])
  # scale by nq and keep it consistent with other nq calculations with na.action = na.pass
  nq <- nrow(stats::model.frame(x, na.action = na.pass))

  return(out / nq)
}

##' @keywords internal
##' @rdname sandwich_elements_calc
.get_tilde_a22_inverse <- function(x, ...) {
  out <- .get_a22_inverse(x, ...)

  if (!inherits(ca <- x$model$`(offset)`, "SandwichLayer")) {
    return(out)
  }

  nq <- nrow(stats::model.frame(x, na.action = na.pass))
  nc_not_q <- sum(!ca@keys$in_Q)
  n <- nq + nc_not_q

  out <- out * n / nq

  return(out)
}

##' @keywords internal
##' @rdname sandwich_elements_calc
.get_tilde_a21 <- function(x) {
  out <- .get_a21(x)

  nq <- nrow(stats::model.frame(x, na.action = na.pass))
  sl <- x$model$`(offset)`
  nc_not_q <- sum(!sl@keys$in_Q)
  n <- nq + nc_not_q

  out <- nq / n * out
}

#' @title (Internal) Design-based variance estimates with HC0 adjustment
#' @param x a fitted \code{teeMod} model
#' @details The design-based variance estimates can be calculated for
#' \code{teeMod} models satisfying the following requirements:
#' - The model uses \code{rct_spec} as \code{StudySpecification}
#' - The model only estimates a main treatment effect
#' - Inverse probability weighting is incorporated
#'
#' @keywords internal
#' @rdname var_estimators
.vcov_DB0 <- function(x, ...) {
  if (!inherits(x, "teeMod")) {
    stop("x must be a teeMod model")
  }

  if (!x@lmitt_fitted){
    # x@lmitt_fitted is false if someone created x using as.lmitt
    stop("x must have been fitted using lmitt.formula")
  }

  if (inherits(os <- x$model$`(offset)`, "SandwichLayer"))
    if (!all(os@keys$in_Q)){
      stop(paste("Design-based standard errors are not supported for teeMod",
                 "models with external sample for covariance adjustment"))
    }

  if (x@StudySpecification@type != "RCT"){
    stop("Design-based standard errors can only be computed for RCT specifications")
  }

  args <- list(...)
  if ("type" %in% names(args)) {
    stop(paste("Cannot override the `type` argument for meat",
               "matrix computations"))
  }
  
  args$x <- if (isTRUE(args$const_effect)) block_center_residuals(x) else x
  args$db <- TRUE
  
  n <- length(args$cluster)

  if (x@absorbed_intercepts) {
    a22inv <- sandwich::bread(x, ...)
    meat <- do.call(sandwich::meatCL, args)
    vmat <- (1 / n) * a22inv %*% meat %*% a22inv
  }
  else {
    if (length(x@moderator) > 0){
      stop(paste("Design-based standard errors are not supported for teeMod",
                 "models with moderators"))
    }
    
    # if model weights does not incorporate IPW, throw a warning
    if (!(inherits(x@lmitt_call$weights, "call") &
          sum(grepl("ate", x@lmitt_call$weights)) > 0)){
      warning(paste("When calculating design-based standard errors,",
                    "please ensure that inverse probability weights are applied.",
                    "This could be done by specifying weights = ate() in",
                    "lmitt() or lm()."))
    }

    if (!is.null(x$call$offset)){
      vmat <- .get_DB_covadj_se(x, ...)
    }
    else {
      vmat <- .get_DB_wo_covadj_se(x)
    }
    
    name <- paste0(var_names(x@StudySpecification, "t"), ".")
    colnames(vmat) <- name
    rownames(vmat) <- name
  }
  attr(vmat, "type") <- "DB0"
  return(vmat)
}

### (If/when DB variance estimators other than DB0 are
###  added, consider replacing this with a stop that
###  invites the user to specify DB0 instead. -BH)
#' @rdname var_estimators
#' @keywords internal
.vcov_DB <- .vcov_DB0

#' @title (Internal) Design-based variance for models with
#'   covariance adjustment
#' @param x a fitted \code{teeMod} model
#' @details Calculate design-based variance estimate for \code{teeMod}
#'   models with covariance adjustment and without absorbed effects
#' @return design-based variance estimate of the main treatment effect
#'   estimate
#' @importFrom stats formula
#' @keywords internal
.get_DB_covadj_se <- function(x, ...){
  if (x@absorbed_intercepts) {
      stop(paste("Design-based standard errors are not supported for\n",
                 "tee models with absorbed intercepts"))
  }
  specification_obj <- x@StudySpecification
  name_trt <- var_names(specification_obj, "t")
  if (name_trt %in% all.vars(stats::formula(x$model$`(offset)`@fitted_covariance_model))) {
    stop(paste("Design-based standard errors are not supported for\n",
               "tee models with treatment in prior covariance adjustment"))
  }

  bread <- .get_DB_covadj_bread(x, ...)

  signs1 <- ifelse(t(t(bread$b1[2,])) %*% bread$b1[2,] > 0, 1, 0)
  signs2 <- ifelse(t(t(bread$b1[2,])) %*% bread$b2[2,] > 0, 1, 0)
  signs3 <- ifelse(t(t(bread$b2[2,])) %*% bread$b2[2,] > 0, 1, 0)

  meat <- .get_DB_covadj_meat(x)

  term1 <- bread$b1 %*% (meat$m1u * signs1 + meat$m1l * (1-signs1)) %*% t(bread$b1)
  term2 <- bread$b2 %*% t(meat$m2u * signs2 + meat$m2l * (1-signs2)) %*% t(bread$b1)
  term3 <- bread$b2 %*% (meat$m3u * signs3 + meat$m3l * (1-signs3)) %*% t(bread$b2)

  vmat <- term1 + 2*term2 + term3
  return(as.matrix(vmat[2,2]))
}

#' @title (Internal) Bread matrix of design-based variance
#' @param x a fitted \code{teeMod} model
#' @details Calculate bread matrix for design-based variance estimate for
#'  \code{teeMod} models with covariance adjustment and without absorbed effects
#' @return a list of bread matrices
#' @keywords internal
.get_DB_covadj_bread <- function(x, ...) {
  a11inv <- .get_a11_inverse(x)
  a21 <- .get_a21(x, vcov_type="DB")
  a22inv <- .get_a22_inverse(x, vcov_type="DB")
  
  # C adjusts for linear transformation of the stacked estimating equation
  # each entry of the estimating equation are transformed to the form: 
  # a_fixed_quantity (under the design-based perspective) * a_treatment_indicator
  C <- matrix(c(1,1,0,1), nrow = 2, byrow = TRUE)

  n <- nrow(x$model)
  bread1 <- a22inv %*% a21 %*% a11inv / n
  bread2 <- - a22inv %*% C / n

  return(list(b1 = bread1, b2 = bread2))
}

#' @title (Internal) Meat matrix of design-based variance
#' @param x a fitted \code{teeMod} model
#' @details Calculate upper and lower bound estimates of meat matrix for
#'   design-based variance estimate for \code{teeMod} models with
#'   covariance adjustment and without absorbed effects
#' @return a list of meat matrix bounds
#' @importFrom stats model.frame
#' @keywords internal
.get_DB_covadj_meat <- function(x, ...) {
  agg <- .aggregate_to_cluster(x)
  data <- agg$data
  bid <- data[, agg$block] # block ids
  zobs <- data[, agg$z] # treatment assignments

  # number of covariates in the covariance model
  p <- ncol(stats::model.matrix(x$model$`(offset)`@fitted_covariance_model))
  XX <- .prepare_spec_matrix(x)
  
  # estimated covariances of entries of \Phi involving indicators z=0 or z=1
  V00 <- .cov_mat_est(XX[zobs==0,], bid[zobs==0])
  V11 <- .cov_mat_est(XX[zobs==1,], bid[zobs==1])
  V01 <- .cov01_est(XX, zobs, bid)
  
  # indices of V
  idl <- (p+2):(2*p+1)
  # upper and lower bounds of the 1st meat matrix term
  meat1u <- V00[1:p, 1:p] + V01[1:p, 1:p] + t(V01[1:p, 1:p]) + V11[1:p, 1:p]
  meat1l <- V00[idl, 1:p] + V01[idl, 1:p] + t(V01[idl, 1:p]) + V11[idl, 1:p]
  
  # upper and lower bounds of the 2nd meat matrix term
  meat2u <- cbind(V00[1:p, p+1], V01[1:p, p+1]) + cbind(V01[p+1, 1:p], V11[1:p, p+1])
  meat2l <- cbind(V00[idl, p+1], V01[idl, p+1]) + cbind(V01[2*p+2, 1:p], V11[idl, p+1])

  # upper and lower bounds of the 3rd meat matrix term
  meat3u <- matrix(c(V00[p+1, p+1], V01[p+1, p+1], V01[p+1, p+1], V11[p+1, p+1]), ncol = 2)
  meat3l <- matrix(c(V00[2*p+2, p+1], V01[2*p+2, p+1], V01[2*p+2, p+1], V11[2*p+2, p+1]), ncol = 2)

  return(list(m1u = meat1u, m1l = meat1l,
              m2u = meat2u, m2l = meat2l,
              m3u = meat3u, m3l = meat3l))
}

#' @title (Internal) Design-based variance for models without
#'   covariance adjustment
#' @param x a fitted \code{teeMod} model
#' @details Calculate bread matrix for design-based variance estimate for
#'   \code{teeMod} models without covariance adjustment and without absorbed
#'   effects
#' @return design-based variance estimate of the main treatment effect
#'   estimate
#' @importFrom stats aggregate var
#' @keywords internal
.get_DB_wo_covadj_se <- function(x, ...){
  if (x@absorbed_intercepts) {
    stop("x should not have absorbed intercepts")
  }
  if (!is.null(x$call$offset)){
    stop("x should not have covariance adjustment")
  }
  # aggregate block ids and outcomes to the cluster level
  agg <- .aggregate_to_cluster(x)
  data <- agg$data
  block <- agg$block

  ws <- data$.w
  yobs <- data$.wy # observed outcomes by cluster
  zobs <- data[, agg$z] # observed treatments by cluster
  bid <- data[, block] # block ids of clusters

  if (length(unique(bid)) == 1) {
    nbk <- t(specification_table(specification=x@StudySpecification, x="treatment"))
  } else {
    nbk <- specification_table(specification=x@StudySpecification, x="treatment",y="block")
  }
  # number of units by block and treatment, B by K
  small_blocks <- apply(nbk, 1, min) == 1 # small block indicator
  nb <- rowSums(nbk) # number of clusters in each block

  rho <- c(sum(ws * yobs * (1-zobs)) / sum(ws * (1-zobs)),
           sum(ws * yobs * zobs) / sum(ws * zobs))

  gammas <- (nbk[bid,1]*(1-zobs) + nbk[bid,2]*zobs) * ws # pseudo outcome gamma
  gamsbk <- list()  # s^2_b,j, sample variance of gamma by block and treatment
  for (k in 1:2){
    indk <- zobs == (k-1)
    gammas[indk] <- gammas[indk] * (yobs[indk] - rho[k])
    gamsbk[[k]] <- stats::aggregate(gammas[indk], by = list(data[indk,block]), FUN = stats::var)
  }
  gamsbk <- merge(gamsbk[[1]], gamsbk[[2]], by = "Group.1")[,2:3]
  gamsbk[is.na(gamsbk)] <- 0
  gamsb <- stats::aggregate(gammas, by = list(bid), FUN = var)[,2] # B vector

  nu1 <- rowSums(gamsbk / nbk) # large block variance estimates
  # small block variance estimates
  nu2 <- 2 / nbk[,1] / nbk[,2] * choose(nb,2) * gamsb -
    (1/nbk[,1] + 1/nbk[,2]) * ((nbk[,1]-1) *gamsbk[,1] + (nbk[,2]-1) *gamsbk[,2])
  varest <- (sum(nu1[!small_blocks]) + sum(nu2[small_blocks])) / sum(data$.w0)^2
  return(as.matrix(varest))
}

#' @title (Internal) Aggregate weights and outcomes to cluster level
#' @param x a fitted \code{teeMod} model
#' @details aggregate individual weights and outcomes to cluster weighted sums
#' @return a list of
#' - a data frame of cluster weights, outcomes, treatments, and block ids;
#' - treatment id column name;
#' - block id column name
#' @keywords internal
.aggregate_to_cluster <- function(x, ...){
  ws <- if (is.null(stats::weights(x))) 1 else stats::weights(x)
  data_temp <- x$call$data
  name_y <- as.character(x$terms[[2]]) # the column of y
  data_temp <- cbind(data_temp, .w = ws, 
                     .w0 = ws / ate(x@StudySpecification, data=data_temp),
                     .wy = ws * data_temp[, name_y])

  specification_obj <- x@StudySpecification
  name_trt <- var_names(specification_obj, "t")
  name_blk <- var_names(specification_obj, "b")
  name_clu <- var_names(specification_obj, "u")

  if (length(name_blk) > 0) {
    data_temp <- .merge_block_id_cols(data_temp, name_blk)
    name_blk <- name_blk[1] # merged block id is stored at column name_blk[1]
  }
  else {
    data_temp$bid <- 1
    name_blk <- 'bid'
  }
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

#' @title (Internal) Merge multiple block IDs
#' @param df a data frame
#' @param ids a vector of block IDs, column names of df
#' @details merge multiple block ID columns by the value combinations and store
#'   the new block ID in the column \code{ids[1]}
#' @return a data frame with a column that contains unique block number IDs
#' @keywords internal
.merge_block_id_cols <- function(df, ids){
  if (!all(ids %in% colnames(df))){
    stop("Some block IDs are not present in the dataframe")
  }
  df[[ids[1]]] <- as.numeric(factor(do.call(paste, c(df[ids], sep = "."))))
  return(df)
}

#' @title (Internal) Helper function for design-based meat matrix
#'   calculation
#' @param x a fitted \code{teeMod} model
#' @return a \eqn{m \times (p+2)} matrix of cluster sums of design-based
#'   estimating equations scaled by \eqn{\sqrt{m_{b0}m_{b1}}/m_{b}}. Here
#'   \eqn{m} is the number of clusters, \eqn{p} is the number of covariates used
#'   in the prior covariance adjustment (excluding intercept)
#' @importFrom stats aggregate
#' @keywords internal
.prepare_spec_matrix <- function(x, ...) {
  specification_obj <- x@StudySpecification
  data <- x$call$data
  name_clu <- var_names(specification_obj, "u")
  # merged cluster id is stored at column name_clu[1]
  cid <- .merge_block_id_cols(data, name_clu)[, name_clu[1]]

  agg <- .aggregate_to_cluster(x)
  bid <- agg$data[, agg$block]  # block ids

  name_y <- as.character(x$terms[[2]]) # name of the outcome column
  # design matrix of the fitted covariance model
  X1 <- model.matrix(x$model$`(offset)`@fitted_covariance_model)
  p <- ncol(X1)

  wc <- x$model$`(offset)`@fitted_covariance_model$weight
  if (is.null(wc)) wc <- 1
  X1 <- matrix(
    rep(wc * (data[, name_y] - x$offset), p), ncol = p
    ) * X1 # n by p, wic * residual * xi

  wi <- if (is.null(stats::weights(x))) 1 else stats::weights(x)
  X2 <- wi * residuals(x)  # n by 1, wi[z] * (residual - rhoz) * z

  XX <- cbind(X1, X2)
  XX <- stats::aggregate(XX, by = list(cid), FUN = sum)[,-1]

  # multiple sqrt(product of #treated divided by block sizes) to XX by group
  if (length(unique(bid)) == 1) {
    nbk <- t(specification_table(specification=x@StudySpecification, x="treatment"))
  } else {
    nbk <- specification_table(specification=x@StudySpecification, x="treatment",y="block")
  }
  const <- sqrt(nbk[, 1] * nbk[, 2] / rowSums(nbk))
  XX <- sweep(XX, 1, const[bid], '*')
  return(XX)
}

#' @title (Internal) Helper function for design-based meat matrix
#'   calculation
#' @details Diagonal elements are estimated by sample variances Off-diagonal
#'   elements are estimated using the Young's elementary inequality
#' @return estimated upper and lower bounds of covariance matrix of estimating
#'   function vectors under either treatment or control
#' @importFrom stats cov
#' @keywords internal
.cov_mat_est <- function(XXz, bidz){
  cov0 <- tapply(1:nrow(XXz), bidz, function(s) stats::cov(XXz[s,]))
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

#' @title (Internal) Helper function for design-based meat matrix calculation
#' @keywords internal
.add_mat_diag <- function(A, B){
  d <- nrow(A)
  A <- matrix(rep(diag(A), d), nrow = d, byrow = FALSE)
  B <- matrix(rep(diag(B), d), nrow = d, byrow = TRUE)
  #return(sqrt(abs(A * B)))
  return((A + B) / 2)
}

#' @title (Internal) Helper function for design-based meat matrix calculation
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

#' @title (Internal) Helper function for design-based meat matrix
#'   calculation
#' @details the Young's elementary inequality is used
#' @return estimated upper and lower bounds of covariance matrix of estimating
#'   function vectors under treatment and under control
#' @importFrom stats cov
#' @keywords internal
.cov01_est <- function(XX, zobs, bid){
  cov0 <- tapply(
    1:nrow(XX[zobs==0,]), bid[zobs==0],
    function(s) stats::cov(XX[zobs==0,][s,])
    )
  cov1 <- tapply(
    1:nrow(XX[zobs==1,]), bid[zobs==1],
    function(s) stats::cov(XX[zobs==1,][s,])
    )
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

#' @title (Internal) Helper function for design-based meat matrix
#'   calculation
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

#' @title (Internal) Calculate grave\{phi\}
#' @keywords internal
#' @param x a fitted \code{teeMod} model
.get_phi_tilde <- function(x, ...){
  specification_obj <- x@StudySpecification
  df <- x$call$data

  # treatment assignments
  name_trt <- var_names(specification_obj, "t")
  assignment <- df[[name_trt]]
  z_ind <- sapply(c(0,1), function(val) as.integer(assignment == val))

  # block ids
  name_blk <- var_names(specification_obj, "b")
  df <- .merge_block_id_cols(df, name_blk)
  name_blk <- name_blk[1]
  stratum <- df[[name_blk]]
  s <- length(unique(stratum)) # number of blocks
  b_ind <- sapply(unique(stratum), function(val) as.integer(stratum == val))

  # nuisance parameters p
  n <- length(stratum) # number of units in Q
  ws <- if (is.null(stats::weights(x))) rep(1,n) else stats::weights(x)
  wb <- matrix(replicate(s, ws), ncol = s) * b_ind
  p <- t(z_ind) %*% wb / (matrix(1, nrow = 2, ncol = n) %*% wb) # 2 is the no. of treatments
  p1 <- p[2, ]

  # calculate phi tilde
  phitilde <- matrix(nrow = n, ncol = s)
  for (j in 1:s){
    phitilde[,j] <- ws * residuals(x) * (z_ind[,2] - p1[j]) * b_ind[,j]
  }
  return(phitilde)
}

#' @title (Internal) Product of \eqn{A_{pp}^{-1} A_{\tau p}^T}
#' @keywords internal
#' @param x a fitted \code{teeMod} model
#' @return An \eqn{s\times k} matrix \eqn{A_{pp}^{-1} A_{\tau p}^T}
.get_appinv_atp <- function(x, ...){
  specification_obj <- x@StudySpecification
  df <- x$call$data

  name_trt <- var_names(specification_obj, "t")
  assignment <- df[[name_trt]]
  k <- length(unique(assignment))

  # stratum ids
  name_blk <- var_names(specification_obj, "b")
  df <- .merge_block_id_cols(df, name_blk)
  name_blk <- name_blk[1]
  stratum <- df[[name_blk]]
  s <- length(unique(stratum)) # number of blocks
  b_ind <- sapply(unique(stratum), function(val) as.integer(stratum == val))
  
  n <- length(stratum) # number of units in the experiment
  ws <- if (is.null(stats::weights(x))) rep(1,n) else stats::weights(x)
  
  app <- list()
  atp <- list()
  for (i in 1:s){
    atp[[i]] <- ws * residuals(x) * b_ind[,i]
    app[[i]] <- ws * b_ind[,i]
  }
  
  mat <- matrix(0, nrow = (k-1)*s, ncol = k-1)
  for (i in 1:s){
    mat[i,1] <- sum(atp[[i]]) / sum(app[[i]])
  }
  return(mat)
}
