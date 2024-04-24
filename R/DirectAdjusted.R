#' @include Design.R WeightedDesign.R DesignAccessors.R SandwichLayerVariance.R confint_lm.R
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

##' @title Compute variance-covariance matrix for fitted \code{DirectAdjusted} model
##' @description
##'   An S3method for \code{stats::vcov} that computes standard errors for
##'   \code{DirectAdjusted} models using \code{vcovDA()}.
##' @details
##'   \code{vcov.DirectAdjusted()} wraps around \code{vcovDA()}, so additional
##'   arguments passed to \code{...} will be passed to the \code{vcovDA()} call.
##'   See documentation for \code{vcovDA()} for information about necessary
##'   arguments.
##' @param object a fitted \code{DirectAdjusted} model
##' @param ... additional arguments to \code{vcovDA()}.
##' @inherit vcovDA return
##' @exportS3Method
vcov.DirectAdjusted <- function(object, ...) {
  cl <- match.call()
  
  cl$x <- cl$object
  argmatch <- match(c("x", "type", "cluster"), names(cl), nomatch = 0L)
  new_cl <- cl[c(1L, argmatch)]
  new_cl[[1L]] <-  quote(vcovDA)

  vmat <- eval(new_cl, parent.frame())
  return(vmat)
}

##' @title Confidence intervals with standard errors provided by \code{vcov.DirectAdjusted()}
##' @description
##'   An S3method for \code{stats::confint} that uses standard errors computed
##'   using \code{vcov.DirectAdjusted()}. Additional arguments passed to this
##'   function, such as \code{cluster} and \code{type}, specify the arguments of
##'   the \code{vcov.DirectAdjusted()} call.
##' @details
##'   Rather than call \code{stats::confint.lm()}, \code{confint.DirectAdjusted()}
##'   calls \code{.confint_lm()}, a function internal to the \code{propertee}
##'   package that ensures additional arguments in the \code{...} of the
##'   \code{confint.DirectAdjusted()} call are passed to the internal \code{vcov()} call.
##' @inheritParams .confint_lm
##' @inherit .confint_lm return
##' @exportS3Method
confint.DirectAdjusted <- function(object, parm, level = 0.95, ...) {
  cl <- match.call()
  cl[[1L]] <- quote(.confint_lm)
  return(eval(cl, parent.frame()))
}

##' @title Extract empirical estimating equations from a \code{DirectAdjusted} model fit
##' @description
##'   An S3method for \code{sandwich::estfun} for producing a matrix of contributions
##'   to the direct adjustment estimating equations.
##' @details
##'   If a prior covariance adjustment model has
##'   been passed to the \code{offset} argument of the \code{DirectAdjusted} model
##'   using \code{cov_adj()}, \code{estfun.DirectAdjusted()} incorporates
##'   contributions to the estimating equations of the covariance adjustment model.\cr\cr
##'   The covariance adjustment sample may not fully overlap with the direct
##'   adjustment sample, in which case \code{estfun.DirectAdjusted()} returns a
##'   matrix with the same number of rows as the number of unique units of observation
##'   used to fit the two models. Uniqueness is determined by matching units of
##'   assignment used to fit the covariance adjustment model to units of assignment
##'   in the \code{DirectAdjusted} model's \code{Design} slot; units of observation
##'   within units of assignment that do not match are additional units that add
##'   to the row count.\cr\cr
##'   The\code{by} argument in \code{cov_adj()} can provide a column or a pair of
##'   columns (a named vector where the name specifies a column in the direct
##'   adjustment sample and the value a column in the covariance adjustment
##'   sample) that uniquely specifies units of observation in each sample. This 
##'   information can be used to align each unit of observation's contributions
##'   to the two sets of estimating equations. If no \code{by} argument is 
##'   provided and units of observation cannot be uniquely specified, contributions
##'   are aligned up to the unit of assignment level. If standard errors are
##'   clustered no finer than that, they will provide the same result as if each
##'   unit of observation's contributions were aligned exactly.
##' 
##' @param x a fitted \code{DirectAdjusted} model
##' @param ... arguments passed to methods
##' @return An \eqn{n\times 2} matrix of empirical
##'  estimating equations for the direct adjustment model fit. See Details for
##'  definition of \eqn{n}.
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
  mat <- estmats[["psi"]] - nq / nc * estmats[["phi"]] %*% t(a11_inv) %*% t(a21)

  return(mat)
}

##' @title Extract bread matrix from a \code{DirectAdjusted} model fit
##' @description
##'   An S3method for \code{sandwich::bread} that extracts the bread of the
##'   direct adjustment model sandwich covariance matrix.
##' 
##' @details This function is a thin wrapper around \code{.get_tilde_a22_inverse()}.
##' @param x a fitted \code{DirectAdjusted} model
##' @param ... arguments passed to methods
##' @inherit vcovDA return
##' @exportS3Method
bread.DirectAdjusted <- function(x, ...) .get_tilde_a22_inverse(x, ...)

##' @title (Internal) Align the dimensions and rows of direct adjustment and covariance
##' adjustment model estimating equations matrices
##' @details
##'   \code{.align_and_extend_estfuns()} first extracts the matrices of contributions
##'   to the empirical estimating equations for the direct adjustment and covariance
##'   adjustment models; then, it pads the matrices with zeros to account for units
##'   of observation that appear in one model-fitting sample but not the other;
##'   finally it orders the matrices so units of observation
##'   (or if unit of observation-level ordering is impossible, units of assignment)
##'   are aligned.
##' @param x a fitted \code{DirectAdjusted} model
##' @param by character vector; indicates unit of assignment columns to generate
##' ID's from; default is NULL, which uses the unit of assignment columns specified
##' in the \code{DirectAdjusted} object's \code{Design} slot
##' @param ... arguments passed to methods
##' @return A list of two matrices, one being the aligned contributions to the
##' estimating equations for the direct adjustment model, and the other being the
##' aligned contributions to the covariance adjustment model.
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

##' @title (Internal) Extract empirical estimating equations from a
##' \code{DirectAdjusted} model using the S3 method associated with its
##' \code{.S3Class} slot
##' @inheritParams estfun.DirectAdjusted
##' @return S3 method
##' @keywords internal
.base_S3class_estfun <- function(x) {
  ## this vector indicates the hierarchy of `sandwich::estfun` methods to use
  ## to extract ITT model's estimating equations
  valid_classes <- c("glm", "lmrob", "svyglm", "lm")
  base_class <- match(x@.S3Class, valid_classes)
  if (all(is.na(base_class))) {
    stop(paste("Direct adjustment model must have been fitted using a function from the",
               "`propertee`, `stats`, `robustbase`, or `survey` package"))
  }
  return(getS3method("estfun", valid_classes[min(base_class, na.rm = TRUE)])(x))
}

#' @title Make ID's to pass to the \code{cluster} argument of \code{vcovDA()}
#' @description
#'   \code{.make_uoa_ids()} returns a factor vector of cluster ID's that align
#'    with the order of the units of observations' contributions in
#'    \code{estfun.DirectAdjusted()}. This is to ensure that when \code{vcovDA()}
#'    calls \code{sandwich::meatCL()}, the \code{cluster} argument aggregates the
#'    correct contributions to estimating equations within clusters.
#' @param x a fitted \code{DirectAdjusted} object
#' @param vcov_type a string indicating model-based or design-based covariance
#' estimation. Currently, "MB", "CR", and "HC" are the only strings registered as
#' indicating model-based estimation.
#' @param cluster character vector or list; optional. Specifies column names that appear
#' in both the covariance adjustment and direct adjustment model dataframes.
#' Defaults to NULL, in which case unit of assignment columns indicated in
#' the Design will be used for clustering. If there are multiple clustering columns,
#' they are concatenated together for each row and separated by "_".
#' @param ... arguments passed to methods
#' @return A vector with length equal to the number of unique units of observation
#' used to fit the two models. See Details of \code{estfun.DirectAdjusted()} for
#' the method for determining uniqueness.
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

#' @title (Internal) Order observations used to fit a \code{DirectAdjusted} model
#' and a prior covariance adjustment model
#' @details \code{.order_samples()} underpins the ordering for \code{.make_uoa_ids()}
#' and \code{estfun.DirectAdjusted()}. This function orders the outputs of those
#' functions, but also informs how the original matrices of contributions to
#' estimating equations need to be indexed to align units of observations'
#' contributions to both sets of estimating equations.\cr\cr When a \code{by}
#' argument is provided to \code{cov_adj()}, it is used to construct the order
#' of \code{.order_samples()}.
#' @param x a fitted \code{DirectAdjusted} model
#' @param by character vector of columns to get ID's for ordering from. Default
#' is NULL, in which case unit of assignment ID's are used for ordering.
#' @param ... arguments passed to methods
#' @return A list of four named vectors. The \code{Q_not_C} element holds the
#' ordering for units of observation in the direct adjustment sample but not the
#' covariance adjustment samples; \code{Q_in_C} and \code{C_in_Q}, the ordering
#' for units in both; and \code{C_not_Q}, the ordering for units in the covariance
#' adjustment sample only. \code{Q_in_C} and \code{C_in_Q} differ in that the
#' names of the \code{Q_in_C} vector correspond to row indices of the original matrix of
#' estimating equations for the direct adjustment model, while the names of
#' \code{C_in_Q} correspond to row indices of the matrix of estimating equations for
#' the covariance adjustment model. Similarly, the names of \code{Q_not_C} and 
#' \code{C_not_Q} correspond to row indices of the direct adjustment and covariance
#' adjustment samples, respectively. Ultimately, the order of \code{.make_uoa_ids()}
#' and \code{estfun.DirectAdjusted()} is given by concatenating the vectors stored
#' in \code{Q_not_C}, \code{Q_in_C}, and \code{C_not_q}.
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

#' @title (Internal) Return ID's used to order observations in the direct adjustment sample
#' @param x a fitted \code{DirectAdjusted} model
#' @param id_col character vector; optional. Specifies column(s) whose ID's will
#' be returned. The column must exist in the data that created the \code{Design}
#' object. Default is NULL, in which case unit of assignment columns indicated
#' in the design will be used to generate ID's.
#' @param ... arguments passed to methods
#' @return A vector with length equal to the number of units of observation in
#' the direct adjustment sample
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
