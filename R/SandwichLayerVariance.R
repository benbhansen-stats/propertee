#' @include Design.R SandwichLayer.R
NULL

# (Internal) Ensure `keys` in a SandwichLayer has NA's for rows where units of
# assignment don't appear in the quasiexperimental sample
.update_keys <- function(x, design) {
  uoanames <- var_names(design, "u")
  zname <- var_names(design, "t")
  check <- .merge_preserve_order(x@keys, design@structure, all.x = TRUE, sort = FALSE)
  convert_na_idx <- which(!is.na(check[, uoanames[1]]) & is.na(check[, zname]))
  x@keys[convert_na_idx, ] <- NA
  return(x)
}

#' @title Generate Cluster-level Estimating Equations Covariance Matrix
#' for Overlapping Rows in Experimental and Covariance Model Data
#' @description `get_overlap_vcov_matrix` returns the covariance matrix of the
#' cluster-level estimating equations for rows that appear in both the experimental
#' and covariance model data
#' @param x DirectAdjusted
#' @return matrix of dim (p+1) x t where p is the number of covariates in the
#' covariance model and t is the number of treatment levels
#' @export
get_overlap_vcov_matrix <- function(x) {
  sl <- x$model$`(offset)`
  if (class(sl) != "SandwichLayer") {
    stop(paste("DirectAdjusted model must have an offset of class `SandwichLayer`",
               "for direct adjustment standard errors"))
  }
  
  uoanames <- var_names(x@Design, "u")
  zname <- var_names(x@Design, "t")

  .update_keys(sl, x@Design)
  
  # Sum est eqns to cluster level; since non-overlapping rows are NA they are
  # excluded automatically in `by` call
  cmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(sl@fitted_covariance_model),
       lapply(uoanames, function(col) sl@keys[, col]),
       colSums))
  
  # look for uoa columns in as.DirectAdjusted data and Design data
  .data <- eval(x@Design@call$data, envir = environment(formula(x@Design@call$formula)))
  # .data <- NULL
  # for (i in seq_len(sys.nframe())) {
  #   try({
  #     .data <- get(x@Design@call$data, sys.frame(i))
  #   }, silent = TRUE)
  #   if (!is.null(.data)) break
  # }
  
  if (is.null(.data)) {
    stop("Could not get quasiexperimental sample from Design")
  }

  # get overlapping rows from experimental data joining with `keys`
  dmod_data <- .merge_preserve_order(
    .data,
    merge(unique(sl@keys), x@Design@structure),
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