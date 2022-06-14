#' @include Design.R SandwichLayer.R
NULL

#' @title Generate Cluster-level Estimating Equations Covariance Matrix
#' for Overlapping Rows in Experimental and Covariance Model Data
#' @description `get_overlap_vcov_matrix` returns the covariance matrix of the
#' cluster-level estimating equations for rows that appear in both the experimental
#' and covariance model data
#' @param x DirectAdjusted
#' @return outer product of cluster-level estimating equations for the covariance model
#' and the treatment, matrix should have dimension p x k where p is the column
#' count in the covariance model matrix and k is the column count in the treatment
#' model
#' @export
get_overlap_vcov_matrix <- function(x) {
  sl <- x$model$`(offset)`
  if (class(sl) != "SandwichLayer") {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }
  
  uoanames <- var_names(x@Design, "u")
  zname <- var_names(x@Design, "t")

  # Sum est eqns to cluster level; since non-overlapping rows are NA they are
  # excluded automatically in `by` call
  cmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(sl@fitted_covariance_model),
       lapply(uoanames, function(col) sl@keys[, col]),
       colSums))
  
  # cmod_eqns will be NULL if there's no overlap, so no need to continue
  if (is.null(cmod_eqns)) {
    return(matrix(0,
                  nrow = dim(stats::model.matrix(sl@fitted_covariance_model))[2],
                  ncol = dim(stats::model.matrix(x))[2]))
  }
  
  # get overlapping rows from experimental data joining with `keys`
  .data <- stats::expand.model.frame(x, uoanames)
  dmod_data <- .merge_preserve_order(
    .data,
    merge(unique(sl@keys), x@Design@structure), # merge here to use txt col for finding NA's
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)
  
  msk <- !is.na(dmod_data[, paste0(zname, ".y")])
  dmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(x)[msk, ],
       lapply(uoanames, function(col) dmod_data[msk, col]),
       colSums))

  return(t(cmod_eqns) %*% dmod_eqns)
}

#' @title Get the "bread" matrix for the Direct Adjustment Model
#' @description The "bread" matrix of the direct adjustment model is the inverse
#' of the model's Fisher information matrix when the gradient is taken with
#' respect to the treatment variable
#' @param x DirectAdjusted
#' @return scalar
#' @export
get_da_bread <- function(x) {
  if (!is(x, "DirectAdjusted")) {
    stop("x must be a DirectAdjusted model")
  }

  zname <- var_names(x@Design, "t")

  if ("glm" %in% x@.S3Class) {
    eqns <- x$family$mu.eta(x$linear.predictors) * x$weights * x$model[, zname]
  } else {
    eqns <- x$weights * x$model[, zname] # for lm's
  }

  return(1 / sum(eqns))
}
