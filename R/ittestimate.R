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
##' @importFrom stats lm predict
ittestimate <- function(design,
                        data,
                        outcome,
                        target = "ate",
                        covAdjModel = NULL,
                        clusterIds = NULL,
                        ...) {
  if (!outcome %in% colnames(data)) {
    stop("outcome must be a column in data")
  }

  if (!is.null(covAdjModel)) {
    covAdj <- cov_adj(covAdjModel)
  }

  if ( !is.null(clusterIds)) {
    if (!is.list(clusterIds) ||
          is.null(names(clusterIds)) ||
          any(names(clusterIds) == "")) {
      stop("clusterIds must be named list")
    }
    if (any(duplicated(names(clusterIds))) || any(duplicated(clusterIds))) {
      stop("clusterIds must be unique")
    }

    # Ensure all names and replacements are valid
    missingnames <- !(names(clusterIds) %in% colnames(design@structure))
    if (any(missingnames)) {
      warning(paste("clusterIds labels not found in Design. unknown elements:",
                 paste(names(clusterIds)[missingnames], collapse = ", ")))
    }
    missingdata <-  !(clusterIds %in% colnames(data))
    if (any(missingdata)) {
      warning(paste("clusterIds replacement values not found in data. unknown elements:",
                   paste(clusterIds[missingnames], collapse = ", ")))
    }

    # if we have any names or replacements missing in the design or data,
    # there's a warning, and then don't try to replace that element
    clusterIds <- clusterIds[!missingnames & !missingdata]


    newnames <- vapply(colnames(design@structure), function(x) {
      pos <- names(clusterIds) == x
      if (any(pos)) {
        return(clusterIds[[which(pos)]])
      }
      return(x)
    }, "character")

    colnames(design@structure) <- newnames
  }


  wtfn <- switch(target,
                 "ate" = ate,
                 "ett" = ett,
                 stop('invalid target'))
  weights <- wtfn(design, data)



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
