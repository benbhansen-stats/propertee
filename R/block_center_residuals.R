#' @title Adjust residuals for both-sides absorption
#' @param x a fitted \code{teeMod} model
#' @return the fitted \code{teeMod} with updated block center residuals
#' @keywords internal
block_center_residuals <- function(x){
  blks <- stats::expand.model.frame(x, var_names(x@Design, "b"))[,var_names(x@Design, "b")]
  blk_means <- rowsum(x$residuals * x$weights, blks) / rowsum(x$weights, blks)
  centered_u <- blk_means[blks]
  x$residuals <- x$residuals - centered_u
  return(x)
}