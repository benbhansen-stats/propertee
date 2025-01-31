#' @include StudySpecification.R WeightedStudySpecification.R StudySpecificationAccessors.R SandwichLayerVariance.R confint_lm.R
NULL
# The above ensures that `StudySpecification`, `WeightedStudySpecification`, and
# `vcov_tee` are defined prior to `teeMod`

setClass("teeMod",
         contains = "lm",
         slots = c(StudySpecification = "StudySpecification",
                   lmitt_fitted = "logical",
                   absorbed_intercepts = "logical",
                   moderator = "character",
                   ctrl_means_model = "lm",
                   lmitt_call = "call"))

setValidity("teeMod", function(object) {
  if (length(object@lmitt_fitted) != 1) {
    return("@lmitt_fitted slot must be a single logical")
  }
  return(TRUE)
})

##' @title Show a \code{teeMod}
##' @description Display information about a \code{teeMod} object
##' @param object \code{teeMod} object, usually a result of a call to
##'   [lmitt()].
##' @return \code{object}, invisibly.
##' @export
setMethod("show", "teeMod", function(object) {
  coeffs <- object$coefficients
  ## Display only treatment effects, intercepts from supplementary
  ## regression of same y but w/o the covadj offset
  if (object@lmitt_fitted) {
    toprint <- (
  ## This should match any coefficients starting with the "txt." or "`txt."
      grepl(paste0("^\\`?", var_names(object@StudySpecification, "t"), "\\."), names(coeffs)) |
  ## This should match the coefficients of the supplementary regression
        grepl(paste0("^", formula(object)[[2L]], ":"), names(coeffs))
    )
    print(coeffs[toprint])
  } else {
    print(coeffs)
  }
  invisible(object)
})

##' @title Compute variance-covariance matrix for fitted \code{teeMod} model
##' @description An S3method for \code{stats::vcov} that computes standard
##'   errors for \code{teeMod} models using \code{vcov_tee()}.
##' @details \code{vcov.teeMod()} wraps around \code{vcov_tee()}, so additional
##'   arguments passed to \code{...} will be passed to the \code{vcov_tee()}
##'   call. See documentation for \code{vcov_tee()} for information about
##'   necessary arguments.
##' @param object a fitted \code{teeMod} model
##' @param ... additional arguments to \code{vcov_tee()}.
##' @inherit vcov_tee return
##' @exportS3Method
vcov.teeMod <- function(object, ...) {
  cl <- match.call()
  cl[[1L]] <-  quote(vcov_tee)
  names(cl)[2L] <- "x"

  vmat <- eval(cl, parent.frame())
  return(vmat)
}

##' @title Confidence intervals with standard errors provided by
##'   \code{vcov.teeMod()}
##' @description An S3method for \code{stats::confint} that uses standard errors
##'   computed using \code{vcov.teeMod()}. Additional arguments passed to this
##'   function, such as \code{cluster} and \code{type}, specify the arguments of
##'   the \code{vcov.teeMod()} call.
##' @details Rather than call \code{stats::confint.lm()},
##'   \code{confint.teeMod()} calls \code{.confint_lm()}, a function internal to
##'   the \code{propertee} package that ensures additional arguments in the
##'   \code{...} of the \code{confint.teeMod()} call are passed to the internal
##'   \code{vcov()} call.
##' @inheritParams .confint_lm
##' @inherit .confint_lm return
##' @exportS3Method
confint.teeMod <- function(object, parm, level = 0.95, ...) {
  cl <- match.call()
  cl[[1L]] <- quote(.confint_lm)
  return(eval(cl, parent.frame()))
}

##' @title Extract empirical estimating equations from a \code{teeMod} model fit
##' @description An S3method for \code{sandwich::estfun} for producing a matrix
##'   of contributions to the direct adjustment estimating equations.
##' @details If a prior covariance adjustment model has been passed to the
##'   \code{offset} argument of the \code{teeMod} model using \code{cov_adj()},
##'   \code{estfun.teeMod()} incorporates contributions to the estimating
##'   equations of the covariance adjustment model.\cr\cr The covariance
##'   adjustment sample may not fully overlap with the direct adjustment sample,
##'   in which case \code{estfun.teeMod()} returns a matrix with the same number
##'   of rows as the number of unique units of observation used to fit the two
##'   models. Uniqueness is determined by matching units of assignment used to
##'   fit the covariance adjustment model to units of assignment in the
##'   \code{teeMod} model's \code{StudySpecification} slot; units of observation
##'   within units of assignment that do not match are additional units that add
##'   to the row count.\cr\cr The\code{by} argument in \code{cov_adj()} can
##'   provide a column or a pair of columns (a named vector where the name
##'   specifies a column in the direct adjustment sample and the value a column
##'   in the covariance adjustment sample) that uniquely specifies units of
##'   observation in each sample. This information can be used to align each
##'   unit of observation's contributions to the two sets of estimating
##'   equations. If no \code{by} argument is provided and units of observation
##'   cannot be uniquely specified, contributions are aligned up to the unit of
##'   assignment level. If standard errors are clustered no finer than that,
##'   they will provide the same result as if each unit of observation's
##'   contributions were aligned exactly.
##'
##' @param x a fitted \code{teeMod} model
##' @param ... arguments passed to methods
##' @return An \eqn{n\times 2} matrix of empirical
##'  estimating equations for the direct adjustment model fit. See Details for
##'  definition of \eqn{n}.
##' @exportS3Method
estfun.teeMod <- function(x, ...) {
  # change model object's na.action to na.exclude so estfun returns NA rows
  if (!is.null(x$na.action)) class(x$na.action) <- "exclude"
  mc <- match.call()
  vcov_type <- match.arg(mc$vcov_type, c("MB", "HC", "CR", "DB")) # this sets the default to model-based

  # MB and DB estfun without covariance adjustment; below, this will be replaced
  # entirely if there's covariance adjustment
  mat <- .base_S3class_estfun(x) - .estfun_DB_blockabsorb(x, ...)
  
  if (vcov_type != "DB") {
    # get estfun for ctrl mean regression 
    cm_ef <- estfun(x@ctrl_means_model)
    q <- ncol(cm_ef)
    cm_ef[is.na(cm_ef)] <- 0
    mat <- cbind(mat, cm_ef)
  } else {
    q <- 0
    cm_ef <- NULL
  }

  ## if ITT model offset doesn't contain info about covariance model, estimating
  ## equations should be the ITT model estimating equations
  if (is.null(sl <- x$model$`(offset)`) | !inherits(sl, "SandwichLayer")) {
    return(mat)
  }

  ## otherwise, extract/compute the rest of the relevant matrices/quantities
  estmats <- .align_and_extend_estfuns(x, cm_ef, ...)
  a11_inv <- .get_a11_inverse(x)
  a21 <- .get_a21(x, ...)

  ## get scaling constants
  nq <- nrow(stats::model.frame(x, na.action = "na.pass")) # this includes NA's as our other routines do
  nc <- nrow(stats::model.frame(sl@fitted_covariance_model))
  n <- nrow(estmats[["psi"]])

  ## form matrix of estimating equations
  mat <- estmats[["psi"]] - nq / nc * estmats[["phi"]] %*% t(a11_inv) %*% t(a21)
  mat[,1:(ncol(estmats[["psi"]]) - q)] <- mat[,1:(ncol(estmats[["psi"]])-q)] - .estfun_DB_blockabsorb(x, ...)
  return(mat)
}

##' @title Extract bread matrix from a \code{teeMod} model fit
##' @description
##'   An S3method for \code{sandwich::bread} that extracts the bread of the
##'   direct adjustment model sandwich covariance matrix.
##'
##' @details This function is a thin wrapper around
##'   \code{.get_tilde_a22_inverse()}.
##' @param x a fitted \code{teeMod} model
##' @param ... arguments passed to methods
##' @inherit vcov_tee return
##' @exportS3Method
bread.teeMod <- function(x, ...) .get_tilde_a22_inverse(x, ...)

##' @title (Internal) Align the dimensions and rows of direct adjustment and
##'   covariance adjustment model estimating equations matrices
##' @details \code{.align_and_extend_estfuns()} first extracts the matrices of
##'   contributions to the empirical estimating equations for the direct
##'   adjustment and covariance adjustment models; then, it pads the matrices
##'   with zeros to account for units of observation that appear in one
##'   model-fitting sample but not the other; finally it orders the matrices so
##'   units of observation (or if unit of observation-level ordering is
##'   impossible, units of assignment) are aligned.
##' @param x a fitted \code{teeMod} model
##' @param by character vector; indicates unit of assignment columns to generate
##'   ID's from; default is NULL, which uses the unit of assignment columns
##'   specified in the \code{teeMod} object's \code{StudySpecification} slot
##' @param ... arguments passed to methods
##' @return A list of two matrices, one being the aligned contributions to the
##'   estimating equations for the direct adjustment model, and the other being
##'   the aligned contributions to the covariance adjustment model.
##' @keywords internal
.align_and_extend_estfuns <- function(x, ctrl_means_ef_mat = NULL, by = NULL, ...) {
  if (!inherits(x, "teeMod") | !inherits(x$model$`(offset)`, "SandwichLayer")) {
    stop("`x` must be a fitted teeMod object with a SandwichLayer offset")
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
  aligned_psi <- matrix(0, nrow = n, ncol = ncol(psi), dimnames = list(seq(n), colnames(psi)))
  aligned_psi[1:nq,] <- psi[Q_order,,drop=FALSE]

  C_order <- c(as.numeric(names(id_order$C_in_Q)), as.numeric(names(id_order$C_not_Q)))
  aligned_phi <- matrix(0, nrow = n, ncol = ncol(phi),
                        dimnames = list(seq_len(n), colnames(phi)))
  aligned_phi[(n-nc+1):n,] <- phi[C_order,,drop=FALSE]
  
  out <- list(psi = aligned_psi, phi = aligned_phi)
  
  if (!is.null(ctrl_means_ef_mat)) {
    aligned_cm_ef <- matrix(0, nrow = n, ncol = ncol(ctrl_means_ef_mat),
                            dimnames = list(seq(n), colnames(ctrl_means_ef_mat)))
    aligned_cm_ef[1:nq,] <- ctrl_means_ef_mat[Q_order,,drop=FALSE]
    out[["psi"]] <- cbind(out[["psi"]], aligned_cm_ef)
  }

  return(out)
}

##' @title (Internal) Extract empirical estimating equations from a
##' \code{teeMod} model using the S3 method associated with its
##' \code{.S3Class} slot
##' @inheritParams estfun.teeMod
##' @return S3 method
##' @keywords internal
.base_S3class_estfun <- function(x) {
  ## this vector indicates the hierarchy of `sandwich::estfun` methods to use
  ## to extract ITT model's estimating equations
  x$coefficients <- replace(x$coefficients, is.na(x$coefficients), 0)
  valid_classes <- c("glm", "lmrob", "svyglm", "lm")
  base_class <- match(x@.S3Class, valid_classes)
  if (all(is.na(base_class))) {
    stop(paste("Direct adjustment model must have been fitted using a function from the",
               "`propertee`, `stats`, `robustbase`, or `survey` package"))
  }
  
  ef <- getS3method("estfun", valid_classes[min(base_class, na.rm = TRUE)])(x)
  ef[is.na(ef)] <- 0
  keep_coef <- if (is.null(aliased <- stats::alias(x)$Complete)) {
    colnames(ef)
  } else {
    if (is.null(nms <- colnames(aliased))) 1 else nms
  }
  return(ef[, keep_coef,drop=FALSE])
}

#' @title Make ID's to pass to the \code{cluster} argument of \code{vcov_tee()}
#' @description \code{.make_uoa_ids()} returns a factor vector of cluster ID's
#'   that align with the order of the units of observations' contributions in
#'   \code{estfun.teeMod()}. This is to ensure that when \code{vcov_tee()} calls
#'   \code{sandwich::meatCL()}, the \code{cluster} argument aggregates the
#'   correct contributions to estimating equations within clusters.
#' @param x a fitted \code{teeMod} object
#' @param vcov_type a string indicating model-based or design-based
#'   covariance estimation. Currently, "MB", "CR", and "HC" are the only strings
#'   registered as indicating model-based estimation.
#' @param cluster character vector or list; optional. Specifies column names
#'   that appear in both the covariance adjustment and direct adjustment model
#'   dataframes. Defaults to NULL, in which case unit of assignment columns
#'   indicated in the StudySpecification will be used for clustering. If there
#'   are multiple clustering columns, they are concatenated together for each
#'   row and separated by "_".
#' @param ... arguments passed to methods
#' @return A vector with length equal to the number of unique units of
#'   observation used to fit the two models. See Details of
#'   \code{estfun.teeMod()} for the method for determining uniqueness.
#' @keywords internal
.make_uoa_ids <- function(x, vcov_type, cluster = NULL, ...) {
  if (!inherits(x, c("teeMod", "mmm"))) {
    stop("Must be provided a teeMod or mmm object")
  } else if (inherits(x, "mmm")) {
    # just use the first model to generate the ID's; this will be the same for
    # mmm objects of models with the same Q (and possibly C) samples
    x <- x[[1L]]
  }
  mc <- match.call()

  # Must be a teeMod object for this logic to occur
  if (!inherits(cluster, "character")) {
    cluster <- var_names(x@StudySpecification, "u")
  }

  # get observation-level unit of assignment and cluster ID's for observations
  # in Q
  Q_obs <- .sanitize_Q_ids(x, id_col = cluster, ...)
  Q_obs_ids <- Q_obs$cluster

  # for model-based vcov calls on blocked specifications when clustering is
  # called for at the assignment level, replace unit of assignment ID's with
  # block ID's for small blocks
  if (vcov_type %in% c("CR", "HC", "MB") & has_blocks(x@StudySpecification)) {
    uoa_cols <- var_names(x@StudySpecification, "u")
    if (length(setdiff(cluster, uoa_cols)) == 0 & length(setdiff(uoa_cols, cluster)) == 0) {
      spec_blocks <- blocks(x@StudySpecification)
      uoa_block_ids <- apply(spec_blocks, 1,
                             function(...) paste(..., collapse = ","))
      small_blocks <- identify_small_blocks(x@StudySpecification)
      structure_w_small_blocks <- cbind(
        x@StudySpecification@structure,
        small_block = small_blocks[uoa_block_ids],
        block_replace_id = apply(spec_blocks, 1,
                                 function(nms, ...) paste(paste(nms, ..., sep = ""), collapse = ","),
                                 nms = colnames(spec_blocks))
      )
      Q_obs <- merge(Q_obs, structure_w_small_blocks, by = uoa_cols, all.x = TRUE)
      na_blocks <- apply(Q_obs[var_names(x@StudySpecification, "b")], 1, function(x) any(is.na(x)))
      Q_obs$cluster[Q_obs$small_block & !na_blocks] <-
        Q_obs$block_replace_id[Q_obs$small_block & !na_blocks]
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
  if (is.null(by <- mc$by)) by <- setdiff(colnames(ca@keys),
                                            c(var_names(x@StudySpecification, "u"), "in_Q"))
  if (length(by) == 0) {
    by <- cluster
  }

  ord_call <- call(".order_samples", quote(x), by = quote(by))
  id_order <- eval(ord_call)

  # if no "by" was specified in cov_adj(), cluster variable was used for ordering,
  # so we can take the names of the sorted vector. Otherwise, we need to get
  # ID's associated with the ordering.
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

#' @title (Internal) Order observations used to fit a \code{teeMod} model and a
#'   prior covariance adjustment model
#' @details \code{.order_samples()} underpins the ordering for
#'   \code{.make_uoa_ids()} and \code{estfun.teeMod()}. This function orders the
#'   outputs of those functions, but also informs how the original matrices of
#'   contributions to estimating equations need to be indexed to align units of
#'   observations' contributions to both sets of estimating equations.\cr\cr
#'   When a \code{by} argument is provided to \code{cov_adj()}, it is used to
#'   construct the order of \code{.order_samples()}.
#' @param x a fitted \code{teeMod} model
#' @param by character vector of columns to get ID's for ordering from. Default
#'   is NULL, in which case unit of assignment ID's are used for ordering.
#' @param ... arguments passed to methods
#' @return A list of four named vectors. The \code{Q_not_C} element holds the
#'   ordering for units of observation in the direct adjustment sample but not
#'   the covariance adjustment samples; \code{Q_in_C} and \code{C_in_Q}, the
#'   ordering for units in both; and \code{C_not_Q}, the ordering for units in
#'   the covariance adjustment sample only. \code{Q_in_C} and \code{C_in_Q}
#'   differ in that the names of the \code{Q_in_C} vector correspond to row
#'   indices of the original matrix of estimating equations for the direct
#'   adjustment model, while the names of \code{C_in_Q} correspond to row
#'   indices of the matrix of estimating equations for the covariance adjustment
#'   model. Similarly, the names of \code{Q_not_C} and \code{C_not_Q} correspond
#'   to row indices of the direct adjustment and covariance adjustment samples,
#'   respectively. Ultimately, the order of \code{.make_uoa_ids()} and
#'   \code{estfun.teeMod()} is given by concatenating the vectors stored in
#'   \code{Q_not_C}, \code{Q_in_C}, and \code{C_not_q}.
#' @keywords internal
.order_samples <- function(x, by = NULL, ...) {
  if (!inherits(x, "teeMod") | !inherits(ca <- x$model$`(offset)`, "SandwichLayer")) {
    stop(paste("x must be a teeMod object with a SandwichLayer offset or",
               "ca must be a SandwichLayer object to retrieve information about",
               "the covariance adjustment model"))
  }
  ## `keys` may have additional columns beyond "in_Q" and the uoa columns if a
  ## `by` argument was specified in `cov_adj()` or `as.SandwichLayer()`. If it
  ## does, use the columns exclusively specified in `by` to merge.
  if (is.null(by) | length(setdiff(colnames(ca@keys),
                                   c(var_names(x@StudySpecification, "u"), "in_Q"))) > 0) {
    by <- setdiff(colnames(ca@keys), c(var_names(x@StudySpecification, "u"), "in_Q"))
  }

  ## Should only hit this if a custom `cluster` argument hasn't been passed to
  ## vcov_tee or no `by` was specified in `cov_adj`
  if (length(by) == 0) {
    by <- var_names(x@StudySpecification, "u")
  }

  # The order, given by the names of the output vector, will be:
  # Q not in C --> Q in C --> C not in Q. The values in the vector correspond to
  # the rows to pull from the original estfun matrices

  # get all ID's in Q
  # Q_ids <- .sanitize_Q_ids(x, id_col = by, ...)[, "cluster"]
  if (x@StudySpecification@unit_of_assignment_type == "none") {
    Q_ids <- rownames(x$model)
  } else {
    Q_ids <- stats::expand.model.frame(x, by)[, by, drop = FALSE]
    Q_ids <- apply(Q_ids, 1, function(...) paste(..., collapse = "_"))
  }

  # get all ID's in C and replace NA's with unique ID
  C_ids <- .sanitize_C_ids(ca, by, sorted = FALSE, ...)

  # need Q_in_C and C_in_Q to have the same order so contributions are aligned
  Q_in_C <- stats::setNames(Q_ids[which(Q_ids %in% C_ids)], which(Q_ids %in% C_ids))
  Q_in_C <- sort(Q_in_C)
  C_in_Q <- stats::setNames(C_ids[which(C_ids %in% Q_ids)], which(C_ids %in% Q_ids))
  C_in_Q <- sort(C_in_Q)

  if (length(Q_in_C) != length(C_in_Q)) {
    stop(paste("Contributions to covariance adjustment and/or effect estimation",
               "are not uniquely specified. Provide a `by` argument to `cov_adj()`",
               "or `vcov.teeMod()`"))
  }

  out <- list(
    Q_not_C = stats::setNames(Q_ids[!(Q_ids %in% C_ids)], which(!(Q_ids %in% C_ids))),
    Q_in_C = Q_in_C,
    C_in_Q = C_in_Q,
    C_not_Q = stats::setNames(C_ids[!(C_ids %in% Q_ids)], which(!(C_ids %in% Q_ids)))
  )

  return(out)
}

#' @title (Internal) Return ID's used to order observations in the direct
#'   adjustment sample
#' @param x a fitted \code{teeMod} model
#' @param id_col character vector; optional. Specifies column(s) whose ID's will
#'   be returned. The column must exist in the data that created the
#'   \code{StudySpecification} object. Default is NULL, in which case unit of
#'   assignment columns indicated in the specification will be used to generate
#'   ID's.
#' @param ... arguments passed to methods
#' @return A vector with length equal to the number of units of observation in
#'   the direct adjustment sample
#' @keywords internal
.sanitize_Q_ids <- function(x, id_col = NULL, ...) {
  # link the units of assignment in the StudySpecification with desired cluster ID's
  uoa_cls_df <- .make_uoa_cluster_df(x@StudySpecification, id_col)
  uoa_cols <- var_names(x@StudySpecification, "u")
  if (nrow(uoa_cls_df) == nrow(x$model)) {
    expand_cols <- unique(c(uoa_cols, id_col))
    by.y <- if (length(expand_cols) == length(uoa_cols)) uoa_cols else c(uoa_cols, "cluster")
  } else {
    expand_cols <- by.y <- uoa_cols
  }

  if (x@StudySpecification@unit_of_assignment_type == "none") {
    moddata <- x$call$data
    moddata$..uoa.. <- rownames(moddata)
    x$call$data <- moddata
  }
  obs_uoa_ids <- stats::expand.model.frame(x,
                                           expand_cols)[, expand_cols, drop = FALSE]

  out <- merge(cbind(obs_uoa_ids, .order_id = seq_len(nrow(obs_uoa_ids))),
               uoa_cls_df, by.x = expand_cols, by.y = by.y, all.x = TRUE, sort = FALSE)

  out <- out[sort(out$.order_id, index.return = TRUE)$ix,]
  out$.order_id <- NULL
  colnames(out) <- unique(c(by.y, "cluster"))
  out$cluster <- as.character(out$cluster)
  rownames(out) <- NULL

  return(out)
}

#' @exportS3Method getCall teeMod
#' @importFrom stats getCall
getCall.teeMod <- function(x, ...) {
  return(x$call)
}

#' @export
update.teeMod <- function(object, ...) {
  stop(paste("teeMod objects do not support an `update` method.\n",
             "You can use `update` on the formula object passed to `lmitt`",
             "instead."))
}

#' @title Design-based estimating equations contributions
#' @param x a fitted \code{teeMod} object
#' @param ... arguments passed to methods
#' @details calculate contributions to empirical estimating equations from a
#'   \code{teeMod} model with absorbed intercepts from the design-based
#'   perspective
#' @return An \eqn{n\times k} matrix
#' @keywords internal
.estfun_DB_blockabsorb <- function (x, by = NULL, ...){
  # return 0 if not asking for a design-based SE or DA does not absorb intercepts
  cl <- match.call()
  vcov_type <- match.arg(cl$vcov_type, c("MB", "HC", "CR", "DB"))
  if (!(db <- vcov_type == "DB") | !x@absorbed_intercepts) {
    return(0)
  }
  
  if (inherits(x$model$`(offset)`, "SandwichLayer")){
    # if the model involves SandwichLayer covariance adjustment
    temp <- .align_and_extend_estfuns(x, by = by, ...)[["psi"]]
  }
  else
    temp <- .base_S3class_estfun(x)

  phitilde <- .get_phi_tilde(x)
  appinv_atp <- .get_appinv_atp(x)

  if (!is.null(x$call$offset)){
    # if the model involves covariance adjustment
    # define the row ordering using .order_samples() and insert rows of 0's into
    # the matrices where necessary while maintaining observation alignment
    id_order <- .order_samples(x, by = by, verbose = FALSE)
    aligned_phitilde <- matrix(
      0, nrow = nrow(temp), ncol = ncol(phitilde),
      dimnames = list(seq_len(nrow(temp)), colnames(phitilde)))

    C_order <- c(as.numeric(names(id_order$C_in_Q)), as.numeric(names(id_order$C_not_Q)))
    aligned_phitilde[which(!is.na(C_order)),] <- phitilde[C_order[!is.na(C_order)], ]
    mat <- aligned_phitilde %*% appinv_atp
  }
  else{
    mat <- phitilde %*% appinv_atp
  }
  mat <- cbind(matrix(0, nrow = nrow(mat), ncol = ncol(temp) - ncol(mat)), mat)
  return(mat)
}
