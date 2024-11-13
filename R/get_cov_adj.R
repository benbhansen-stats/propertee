#' @include StudySpecification.R
NULL

# (Internal) Extract SandwichLayer from model. Can exist in two places,
# depending on whether model was called as `lm(y ~ z, offset = cov_adj(spec))`
# or `lm(y ~ z + offset(cov_adj(spec))`. Returns `NULL` if it can't find it.
.get_cov_adj <- function(x) {

  # Look for an `offset = ` argument
  offset_from_arg <- x$model$"(offset)"

  if (inherits(offset_from_arg, "PreSandwichLayer") ||
      inherits(offset_from_arg, "SandwichLayer")) {
    return(offset_from_arg)
  }

  # Look for a `+ offset(` term in the formula
  offset_col <- x$model[grepl("offset\\(", names(x$model))]

  which_are_SL <- vapply(offset_col,
                         function(x) {
                           inherits(x, "PreSandwichLayer") ||
                             inherits(x, "SandwichLayer")
                          },
                         logical(1))

  if (sum(which_are_SL) == 1) {
    # Most of the time there will only be one offset, so return it
    return(offset_col[, which_are_SL])
  } else if (sum(which_are_SL) > 1) {
    if (length(unique(as.list(offset_col))) == 1) {
      # All cov_adj are the same, return any of them
      return(offset_col[, 1])
    }
    # At least one cov_Adj is different, error
    stop("Multiple cov_adj() calls found in offset")
  }

  return(NULL)
}
