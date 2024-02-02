##' Group-center akin to Stata's \code{areg}
##'
##' From the Stata documentation: \code{areg} begins by recalculating Y and X
##' and to have mean 0 within the groups specified by \code{absorb()}. The
##' overall mean of each variable is then added back in.
##' @param mm Matrix of variables to center
##' @param grp Group to center on
##' @param wts Optional weights
##' @param grand_mean_center Optional center output at \code{mean(var)}
##' @return Vector of group-centered values
##' @keywords internal
areg.center <- function(mm, grp, wts = NULL, grand_mean_center = FALSE) {
  if (is.null(wts)) {
    wts <- rep(1, nrow(mm))
  }
  mm2 <- mm[!is.na(wts),,drop=FALSE]
  grp_keep <- !is.na(grp)
  group_means <- sweep(rowsum((wts * mm2)[grp_keep,], grp[grp_keep], na.rm = TRUE),
                       1, rowsum(wts[grp_keep], grp[grp_keep], na.rm = TRUE), FUN = "/")
  group_means[is.nan(group_means)] <- 0
  out <- mm2 - group_means[grp,]
  if (grand_mean_center) {
    grand_means <- colSums(wts * mm2, na.rm = TRUE) / sum(wts)
    out <- sweep(out, 2, grand_means, FUN = "+")
  }
  return(out)
}
