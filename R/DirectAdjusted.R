#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R
NULL
# The above ensures that `Design`, `WeightedDesign`, and `vcovDA` are defined
# prior to `DirectAdjusted`

setClass("DirectAdjusted",
         contains = "lm",
         slots = c(Design = "Design",
                   lmitt_fitted = "logical",
                   absorbed_intercepts = "logical",
                   moderator = "character",
                   lmitt_call = "call"))

setValidity("DirectAdjusted", function(object) {
  if (length(object@lmitt_fitted) != 1) {
    return("@lmitt_fitted slot must be a single logical")
  }
  return(TRUE)
})

##' @title Show an DirectAdjusted
##' @param object DirectAdjusted object
##' @return an invisible copy of `object`
##' @export
setMethod("show", "DirectAdjusted", function(object) {
  coeffs <- object$coefficients
  # Display only treatment effects
  if (object@lmitt_fitted) {
    # This should match any coefficients starting with the "txt." or "`txt."
    toprint <- grepl(paste0("^\\`?", var_names(object@Design, "t"), "\\."),
                     names(coeffs))
    print(coeffs[toprint])
  } else {
    print(coeffs)
  }
  invisible(object)
})

##' @title Variance-Covariance matrix of \code{DirectAdjusted} object
##' @details If a \code{DirectAdjusted} object is fit with a \code{SandwichLayer}
##' offset, then its \code{vcov()} method provides a sandwich estimate of the
##' covariance-adjusted variance-covariance matrix. Otherwise, it provides
##' the default OLS estimate of the matrix.
##' @param object DirectAdjusted
##' @param ... Additional arguments to \code{vcovDA()} or \code{stats:::vcov.lm()}.
##' @return Variance-Covariance matrix
##' @exportS3Method
vcov.DirectAdjusted <- function(object, ...) {
  cl <- match.call()

  if (is.null(cl[["type"]])) {
    confint_calls <- grepl("confint.DirectAdjusted", lapply(sys.calls(), "[[", 1))
    if (any(confint_calls)) {
      type <- tryCatch(get("call", sys.frame(which(confint_calls)[1]))$type,
                       error = function(e) NULL)
      cl$type <- type # will not append if type is NULL
    }
  }

  cl$x <- cl$object
  argmatch <- match(c("x", "type", "cluster"), names(cl), nomatch = 0L)
  new_cl <- cl[c(1L, argmatch)]
  new_cl[[1L]] <-  quote(vcovDA)

  vmat <- eval(new_cl, parent.frame())
  return(vmat)
}

##' @title Variance-Covariance matrix
##' @param object DirectAdjusted
##' @param parm a specification of which parameters are to be given confidence
##'   intervals, either a vector of numbers or a vector of names. If missing,
##'   all parameters are considered.
##' @param level the confidence level required.
##' @param ... Add'l arguments
##' @return Variance-Covariance matrix
##' @exportS3Method
confint.DirectAdjusted <- function(object, parm, level = 0.95, ...) {
  call <- match.call()
  call[[1L]] <- quote(stats::confint.lm)

  ci <- eval(call, parent.frame())
  return(ci)
}

##' @title Extract empirical estimating equations from a \code{DirectAdjusted} model fit
##' @param x a fitted \code{DirectAdjusted} object
##' @param ... arguments passed to methods
##' @return An \eqn{n\times k} matrix of empirical estimating equations for the
##' covariance-adjusted ITT effect regression. \eqn{n} represents the number of
##' observations in the union of the samples used to fit the two regressions.
##' \eqn{k} represents the number of parameters in the latter model.\cr\cr
##' Each row represents an observation's contribution to the stacked estimating
##' equations. This contribution, denoted \eqn{\tilde{\psi}_{i}} for the \eqn{i}th
##' observation, is given by \deqn{\tilde{\psi}_{i} = \psi_{i} +
##' \phi_{i}A_{11}^{-1}A_{21}^{T}} where \eqn{\psi_{i}} is the observation's
##' contribution to the ITT effect model fit, \eqn{\phi_{i}} is the observation's
##' contribution to the covariance adjustment model fit, and the \eqn{A} matrices
##' are given by typical sandwich calculations.\cr\cr
##' Note that the formulation of the output matrix \eqn{\tilde{\Psi}} requires
##' information about each observation's contributions to both the covariance
##' adjustment and ITT effect models (where some observations may not contribute
##' to both models). Estimating equations are taken from \code{sandwich::estfun}
##' calls on both fitted models and aligned as closely as possible.
##' The `by` argument in `cov_adj()` can be used to specify a column unrelated to
##' the design that will allow for exact alignment of such matrices. If no `by`
##' argument is provided, clustering information given in the \code{DirectAdjusted}'s
##' \code{Design} object will be used to align rows by unit of assignment, even
##' though no guarantees can be made about aligning the matrices within units
##' of assignment. Regardless of the eventual alignment and initial ordering of the
##' observations in the two matrices, however, when using \code{vcovDA},
##' variance estimates will ultimately be the same due to the clustering passed
##' to any \code{sandwich::meatCL} calls.
##' @rdname da_estfun_methods
##' @exportS3Method
estfun.DirectAdjusted <- function(x, ...) {
  ## if ITT model offset doesn't contain info about covariance model, estimating
  ## equations should be the ITT model estimating equations
  if (is.null(sl <- x$model$`(offset)`) | !inherits(sl, "SandwichLayer")) {
    return(.base_S3class_estfun(x))
  }

  ## otherwise, extract/compute the rest of the relevant matrices/quantities
  estmats <- .align_and_extend_estfuns(x, ...)
  a11_inv <- .get_a11_inverse(x)
  a21 <- .get_a21(x)

  ## get scaling constants
  nq <- nrow(stats::model.frame(x))
  nc <- nrow(stats::model.frame(sl@fitted_covariance_model))
  n <- nrow(estmats[["psi"]])

  ## form matrix of estimating equations
  mat <- estmats[["psi"]] - (nq / sqrt(nc * n)) * estmats[["phi"]] %*% t(a11_inv) %*% t(a21)

  return(mat)
}

##' @title Extract bread matrix from a \code{DirectAdjusted} model fit
##' @details This function is a thin wrapper around \code{.get_tilde_a22_inverse()}.
##' @param x a fitted \code{DirectAdjusted} object
##' @param ... arguments passed to methods
##' @inherit .get_a22_inverse return
##' @rdname da_estfun_methods
##' @exportS3Method
bread.DirectAdjusted <- function(x, ...) .get_tilde_a22_inverse(x, ...)

##' (Internal) Align the dimensions and rows of estimating equations matrices
##' from the ITT effect and covariance adjustment models
##' @param x a fitted \code{DirectAdjusted} object
##' @param by character vector; indicates unit of assignment columns to generate
##' ID's from; default is NULL, which uses the unit of assignment columns specified
##' in the \code{DirectAdjusted} object's \code{Design} slot
##' @param ... arguments passed to methods
##' @return list of two matrices, one being the aligned contributions to the
##' estimating equations for the ITT effect model, and the other being the
##' aligned contributions to the covariance adjustment model
##' @keywords internal
.align_and_extend_estfuns <- function(x, by = NULL, ...) {
  if (!inherits(x, "DirectAdjusted") | !inherits(x$model$`(offset)`, "SandwichLayer")) {
    stop("`x` must be a fitted DirectAdjusted object with a SandwichLayer offset")
  }

  # get the original estimating equations
  psi <- .base_S3class_estfun(x)
  phi <- estfun(x$model$`(offset)`@fitted_covariance_model)

  # the ordering output by `.order_samples()` is explained in that function's
  # documentation
  id_order <- .order_samples(x, by = by, verbose = FALSE)
  n <- length(c(id_order$Q_not_C, id_order$C_in_Q, id_order$C_not_Q))
  nq <- length(c(id_order$Q_not_C, id_order$Q_in_C))
  nc <- length(c(id_order$C_in_Q, id_order$C_not_Q))

  Q_order <- c(as.numeric(names(id_order$Q_not_C)), as.numeric(names(id_order$Q_in_C)))
  aligned_psi <- matrix(0, nrow = n, ncol = ncol(psi),
                        dimnames = list(seq(n), colnames(psi)))
  aligned_psi[1:nq,] <- psi[Q_order,,drop=FALSE]

  C_order <- c(as.numeric(names(id_order$C_in_Q)), as.numeric(names(id_order$C_not_Q)))
  aligned_phi <- matrix(0, nrow = n, ncol = ncol(phi),
                        dimnames = list(seq_len(n), colnames(phi)))
  aligned_phi[(n-nc+1):n,] <- phi[C_order,,drop=FALSE]

  return(list(psi = aligned_psi, phi = aligned_phi))
}

##' (Internal) Call \code{sandwich::estfun} method for a fitted \code{DirectAdjusted}
##' object based on its base S3 class
##' @inheritParams estfun.DirectAdjusted
##' @return S3 method
##' @keywords internal
.base_S3class_estfun <- function(x) {
  ## this vector indicates the hierarchy of `sandwich::estfun` methods to use
  ## to extract ITT model's estimating equations
  valid_classes <- c("glm", "lmrob", "svyglm", "lm")
  base_class <- match(x@.S3Class, valid_classes)
  if (all(is.na(base_class))) {
    stop(paste("ITT effect model must have been fitted using a function from the",
               "`propertee`, `stats`, `robustbase`, or `survey` package"))
  }
  return(getS3method("estfun", valid_classes[min(base_class, na.rm = TRUE)])(x))
}

#' Make unit of assignment ID's to pass to \code{sandwich::meatCL} `cluster`
#' argument
#' @details These ID's align with the output of \code{estfun.DirectAdjusted}. If
#' a \code{by} argument was used for \code{.order_samples}, \code{.make_uoa_ids}
#' will return the values of the columns specified in \code{cluster} associated
#' with that ordering.
#' @param x a fitted \code{DirectAdjusted} object
#' @param vcov_type a string indicating model-based or design-based covariance
#' estimation. Currently, "MB", "CR", and "HC" are the only strings registered as
#' indicating model-based estimation.
#' @param cluster character vector or list; optional. Specifies column names that appear
#' in both the covariance adjustment dataframe C and quasiexperimental dataframe
#' Q. Defaults to NULL, in which case unit of assignment columns indicated in
#' the Design will be used to generate clustered covariance estimates. If there
#' are multiple clustering columns, they are concatenated together for each row
#' and separated by "_".
#' @param ... arguments passed to methods
#' @return A vector of length \eqn{|\mathcal{Q}| + |\mathcal{C} \ \mathcal{Q}|}
#' @keywords internal
.make_uoa_ids <- function(x, vcov_type, cluster = NULL, ...) {
  if (!inherits(x, "DirectAdjusted")) {
    stop("Must be provided a DirectAdjusted object")
  }

  # Must be a DirectAdjusted object for this logic to occur
  if (!inherits(cluster, "character")) {
    cluster <- var_names(x@Design, "u")
  }

  # get observation-level unit of assignment and cluster ID's for observations in Q
  Q_obs <- .sanitize_Q_ids(x, id_col = cluster, ...)
  Q_obs_ids <- Q_obs$cluster

  # for model-based vcov calls on blocked designs when clustering is called for
  # at the assignment level, replace unit of assignment ID's with block ID's
  # for small blocks
  if (vcov_type %in% c("CR", "HC", "MB") & has_blocks(x@Design)) {
    uoa_cols <- var_names(x@Design, "u")
    if (length(setdiff(cluster, uoa_cols)) == 0 & length(setdiff(uoa_cols, cluster)) == 0) {
      des_blocks <- blocks(x@Design)
      # uoa_block_ids <- apply(des_blocks, 1, function(...) paste(..., collapse = ","))
      uoa_block_ids <- apply(des_blocks, 1,
                             function(...) paste(..., collapse = ","))
      small_blocks <- identify_small_blocks(x@Design)
      structure_w_small_blocks <- cbind(
        x@Design@structure,
        small_block = small_blocks[uoa_block_ids],
        block_replace_id = apply(des_blocks, 1,
                                 function(nms, ...) paste(paste(nms, ..., sep = ""), collapse = ","),
                                 nms = colnames(des_blocks))
      )
      Q_obs <- merge(Q_obs, structure_w_small_blocks, by = uoa_cols, all.x = TRUE)
      Q_obs$cluster[Q_obs$small_block] <- Q_obs$block_replace_id[Q_obs$small_block]
      Q_obs_ids <- Q_obs$cluster
    }
  }

  # If there's no covariance adjustment info, return the ID's found in Q
  if (!inherits(ca <- x$model$`(offset)`, "SandwichLayer")) {
    return(factor(Q_obs_ids, levels = unique(Q_obs_ids)))
  }

  # `keys` may have columns in addition to "in_Q" and the uoa columns if a
  # `by` argument was specified in `cov_adj()` or `as.SandwichLayer()`. If it
  # does, use the columns exclusively specified in `by` to produce the order
  by <- setdiff(colnames(ca@keys), c(var_names(x@Design, "u"), "in_Q"))
  if (length(by) == 0) {
    by <- cluster
  }
  id_order <- .order_samples(x, by = by, ...)

  # if no "by" was specified in cov_adj(), cluster variable was used for ordering,
  # so we can take the names of the sorted vector. Otherwise, we need to get
  # ID's associated with the ordering
  if (length(setdiff(cluster, by)) == 0) {
    ids <- Reduce(c, id_order[c("Q_not_C", "Q_in_C", "C_not_Q")])
  } else {
    C_ids <- .sanitize_C_ids(ca, cluster, verbose = FALSE, sorted = FALSE, ...)
    ids <- c(
      Q_obs_ids[as.numeric(names(id_order$Q_not_C))],
      Q_obs_ids[as.numeric(names(id_order$Q_in_C))],
      C_ids[as.numeric(names(id_order$C_not_Q))]
    )
  }

  na_ids <- is.na(ids)
  ids[na_ids] <- apply(
    matrix(sample(c(letters, LETTERS), 8 * sum(na_ids), replace = TRUE), ncol = 8),
    1, function(...) paste(..., collapse = "")
  )

  return(factor(ids, levels = unique(ids)))
}

#' Order observations used to fit a \code{DirectAdjusted} model and its
#' covariance adjustment model
#' @details \code{.order_samples} underpins the ordering for \code{.make_uoa_ids}
#' and \code{estfun.DirectAdjusted}. This function needs to order the rows in
#' in \eqn{\mathcal{Q}\cup\mathcal{C}}, but it also needs to explain how the
#' original matrices of estimating equations should be indexed so the contributions
#' from both sets of equations match. So instead of returning a numeric vector,
#' which one might expect for an ordering function, \code{.order_samples}
#' returns a list of vectors, which is explained in the Return section. Ultimately,
#' the order is given by concatenating the vectors stored in \code{Q_not_C},
#' \code{Q_in_C}, and \code{C_not_q} (\code{Q_in_C} and \code{C_in_Q} are
#' interchangeable in terms of deriving the order). The names of the
#' \code{Q_not_C} and \code{Q_in_C} vectors correspond to row indices of the
#' matrix of estimating equations for the ITT effect model, while the names of
#' the \code{C_in_Q} and \code{C_not_Q} vectors correspond to row indices of
#' the matrix of estimating equations for the covariance adjustment model.
#' When a \code{by} argument is provided to \code{cov_adj}, it is used to
#' deduce the order.
#' @param x a fitted \code{DirectAdjusted} object
#' @param by character vector; indicates unit of assignment columns to generate
#' ID's from; default is NULL, which uses the unit of assignment columns specified
#' in the \code{DirectAdjusted} object's \code{Design} slot
#' @return a list of four named vectors
#' @param ... arguments passed to methods
#' @keywords internal
.order_samples <- function(x, by = NULL, ...) {
  if (!inherits(x, "DirectAdjusted") | !inherits(ca <- x$model$`(offset)`, "SandwichLayer")) {
    stop(paste("x must be a DirectAdjusted object with a SandwichLayer offset or",
               "ca must be a SandwichLayer object to retrieve information about",
               "the covariance adjustment model"))
  }
  ## `keys` may have additional columns beyond "in_Q" and the uoa columns if a
  ## `by` argument was specified in `cov_adj()` or `as.SandwichLayer()`. If it
  ## does, use the columns exclusively specified in `by` to merge.
  if (is.null(by) | length(setdiff(colnames(ca@keys),
                                   c(var_names(x@Design, "u"), "in_Q"))) > 0) {
    by <- setdiff(colnames(ca@keys), c(var_names(x@Design, "u"), "in_Q"))
  }

  if (length(by) == 0) {
    by <- var_names(x@Design, "u")
  }

  # The order, given by the names of the output vector, will be:
  # Q not in C --> Q in C --> C not in Q. The values in the vector correspond to
  # the rows to pull from the original estfun matrices

  # get all ID's in Q
  Q_ids <- .sanitize_Q_ids(x, id_col = by, ...)[, "cluster"]

  # get all ID's in C and replace NA's with unique ID
  C_ids <- .sanitize_C_ids(ca, by, sorted = FALSE, ...)

  # need Q_in_C and C_in_Q to have the same order so contributions are aligned
  Q_in_C <- stats::setNames(Q_ids[which(Q_ids %in% C_ids)], which(Q_ids %in% C_ids))
  Q_in_C <- sort(Q_in_C)
  C_in_Q <- stats::setNames(C_ids[which(ca@keys$in_Q)], which(ca@keys$in_Q))
  C_in_Q <- sort(C_in_Q)

  out <- list(
    Q_not_C = stats::setNames(Q_ids[!(Q_ids %in% C_ids)], which(!(Q_ids %in% C_ids))),
    Q_in_C = Q_in_C,
    C_in_Q = C_in_Q,
    C_not_Q = stats::setNames(C_ids[!ca@keys$in_Q], which(!ca@keys$in_Q))
  )

  return(out)
}

#' @title Return ID's for observations in \eqn{Q}
#' @param x a fitted \code{DirectAdjusted} object
#' @param id_col character vector; optional. Specifies column(s) whose ID's will
#' be returned. The column must exist in the data that created the \code{Design}
#' object. Default is NULL, in which case unit of assignment columns indicated
#' in the design will be used to generate ID's.
#' @param ... arguments passed to methods
#' @return A vector of length \eqn{|Q|}
#' @keywords internal
.sanitize_Q_ids <- function(x, id_col = NULL, ...) {
  # link the units of assignment in the Design with desired cluster ID's
  uoa_cls_df <- .make_uoa_cluster_df(x@Design, id_col)
  uoa_cols <- var_names(x@Design, "u")
  if (nrow(uoa_cls_df) == nrow(x$model)) {
    expand_cols <- unique(c(uoa_cols, id_col))
    by.y <- if (length(expand_cols) == length(uoa_cols)) uoa_cols else c(uoa_cols, "cluster")
  } else {
    expand_cols <- by.y <- uoa_cols
  }

  obs_uoa_ids <- stats::expand.model.frame(x,
                                           expand_cols,
                                           na.expand = TRUE)[, expand_cols, drop = FALSE]

  out <- merge(cbind(obs_uoa_ids, .order_id = seq_len(nrow(obs_uoa_ids))),
               uoa_cls_df, by.x = expand_cols, by.y = by.y, all.x = TRUE, sort = FALSE)

  out <- out[sort(out$.order_id, index.return = TRUE)$ix,]
  out$.order_id <- NULL
  colnames(out) <- unique(c(by.y, "cluster"))
  out$cluster <- as.character(out$cluster)
  rownames(out) <- NULL

  return(out)
}

#' @exportS3Method getCall DirectAdjusted
#' @importFrom stats getCall
getCall.DirectAdjusted <- function(x, ...) {
  return(x$call)
}

#' @export
update.DirectAdjusted <- function(object, ...) {
  stop(paste("DirectAdjusted objects do not support an `update` method.\n",
             "You can use `update` on the formula object passed to `lmitt`",
             "instead."))
}
