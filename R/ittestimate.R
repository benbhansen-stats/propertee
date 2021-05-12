##' @title Estimate Treatment Effect
##' @param design study Design
##' @param data Data for analysis
##' @param outcome string containing name of outcome variable in `data`
##' @param target a function returning a WeightedDesign object, or either "ate"
##'   (default) or "ett"
##' @param covAdjModel optional; covariable adjustment model
##' @param clusterIds optional; list connecting names of cluster variables in
##'   `design` to cluster variables in `data`
##' @param ... Additional arguments to the model.
##' @return Estimated treatment effect
##' @export
##' @importFrom stats lm
ittestimate <- function(design,
                        data,
                        outcome,
                        target = ate,
                        covAdjModel = NULL,
                        clusterIds = NULL,
                        ...) {
  if (!outcome %in% colnames(data)) {
    stop("outcome must be a column in data")
  }

  if (is.function(target)) {
    wtfn <- target
  } else if (is.character(target)) {
    wtfn <- switch(target,
                   "ate" = ate,
                   "ett" = ett,
                   stop(paste('invalid target; must be function returning',
                              'a WeightedDesign object, or either "ett" or "ate"')))
  } else {
    stop(paste('invalid target; must be function returning a WeightedDesign',
               'object, or either "ett" or "ate"'))
  }
  weights <- wtfn(design, data)





  # Expand treatment status
  ctdata <- design@structure[, design@columnIndex %in% c("t", "c")]
  colnames(ctdata)[1] <- "Design_Treatment"
  merged <- merge(data, ctdata, by = colnames(ctdata)[-1])

  model <- lm(merged[, outcome] ~ Design_Treatment, data = merged, weights = weights, ...)

  return(model)
}
