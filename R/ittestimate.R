##' @title Estimate Treatment Effect
##' @param design study Design
##' @param data Data for analysis
##' @param outcome string containing name of outcome variable in `data`
##' @param target a function returning a WeightedDesign object, or either "ate"
##'   (default) or "ett"
##' @param covAdjModel optional; covariable adjustment model
##' @param clusterIds optional; list connecting names of cluster variables in
##'   `design` to cluster variables in `data`
##' @param weights optional; manually include weights. If included, weights will
##'   **not** be automatically generated. Instead, the `weights` argument should
##'   include the product of any externally created weights and the result of
##'   `ate()` or `ett()`. For example, `weights = myweights*ate(des)`.
##' @param ... Additional arguments to the model.
##' @return Estimated treatment effect
##' @export
##' @importFrom stats lm predict
ittestimate <- function(design,
                        data,
                        outcome,
                        target = "ate",
                        covAdjModel = NULL,
                        clusterIds = NULL,
                        ...,
                        weights = NULL) {
  if (!outcome %in% colnames(data)) {
    stop("outcome must be a column in data")
  }

  if (!is.null(covAdjModel)) {
    covAdj <- cov_adj(covAdjModel)
  }

  if (!is.null(clusterIds)) {
    design <- update_clusterIds(design, data, clusterIds)
  }

  if (is.null(weights)) {
    wtfn <- switch(target,
                   "ate" = ate,
                   "ett" = ett,
                   stop('invalid target'))
    weights <- wtfn(design, data)
  } else {
    if (length(weights) != nrow(data)) {
      stop("weights must be same length as rows of `data`")
    }
    if (!is.numeric(weights)) {
      stop("weights must be numeric")
    }
  }


  # Expand treatment status
  ctdata <- design@structure[, design@columnIndex %in% c("t", "c")]
  colnames(ctdata)[1] <- "Design_Treatment"
  merged <- merge(data, ctdata, by = colnames(ctdata)[-1])

  if (!is.null(covAdjModel)) {
    model <- lm(merged[, outcome] ~ Design_Treatment + offset(covAdj),
                data = merged, weights = weights, ...)
  } else {
    model <- lm(merged[, outcome] ~ Design_Treatment,
                data = merged, weights = weights, ...)
  }

  return(DirectAdjusted(model, Design = design, target = target))
}
