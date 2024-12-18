#' @title Adjust residuals for both-sides absorption
#' @param x a fitted \code{teeMod} model
#' @details
#' This function subtracts off the block residual mean function
#' \eqn{\hat \alpha(v_b, \theta)} for each observation from model residuals
#' @return the fitted \code{teeMod} with updated block center residuals.
#' @keywords internal
block_center_residuals <- function(x){
  blks <- stats::expand.model.frame(
    x, var_names(x@StudySpecification, "b")
    )[names(residuals(x)),var_names(x@StudySpecification, "b")]
  n <- length(blks)
  w <- if (is.null(weights(x))) rep(1, n) else weights(x)
  blk_means <- suppressWarnings(rowsum(
    residuals(x) * w, blks, na.rm = TRUE) / rowsum(w, blks, na.rm = TRUE))
  centered_u <- blk_means[blks]
  x$residuals <- residuals(x) - centered_u
  return(x)
}
