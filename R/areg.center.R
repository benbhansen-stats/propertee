##' Group-center akin to Stata's \code{areg}
##'
##' From the Stata documentation: \code{areg} begins by recalculating Y and X
##' and to have mean 0 within the groups specified by \code{absorb()}. The
##' overall mean of each variable is then added back in.
##' @param var Variable to center
##' @param grp Group to center on
##' @param wts Optional weights
##' @return Vector of group-centered values
##' @keywords internal
##' @importFrom stats weighted.mean
areg.center <- function(var, grp, wts = NULL) {
  if (!is.null(wts)) {
    # weighted.mean produces NA if any weights are NA
    if (any(is.na(wts))) {
      var2 <- var[!is.na(wts)]
    } else {
      var2 <- var
    }

    df <- data.frame(var = var, wts = wts)

    group_means <- sapply(split(df, grp), function(x) {
      stats::weighted.mean(x$var, x$wts, na.rm = TRUE)
    })
    # #140 - if all weights are 0, weighted.mean returns NaN
    group_means[is.nan(group_means)] <- 0

    # Expand to original data
    group_means <- group_means[as.character(grp)]


    grand_mean <- stats::weighted.mean(var2, w = wts[!is.na(wts)], na.rm = TRUE)
    out <- var - group_means + grand_mean
  } else {
    out <- var - tapply(var, grp, mean, na.rm = TRUE)[as.character(grp)] +
      mean(var, na.rm = TRUE)
  }
  return(as.vector(out))
}
