##' @title Estimate Treatment Effect
##' @param design study Design
##' @param data Data for analysis
##' @param outcome string containing name of outcome variable in \code{data}
##' @param target a function returning a WeightedDesign object, or either "ate"
##'   (default) or "ett"
##' @param cov_adj_model optional; covariable adjustment model
##' @param by optional; list connecting names of units of
##'   assignment/clusters/units variables in \code{design} to units of
##'   assignment/clusters/units variables in \code{data}
##' @param weights optional; manually include weights. If included, the weights
##'   specified in \code{target} will **not** be automatically generated.
##'   Instead, the \code{weights} argument should include the product of any
##'   externally created weights and the result of \code{ate()} or \code{ett()}.
##'   For example, \code{weights = myweights*ate(des)}.
##' @param ... Additional arguments to the model.
##' @return Estimated treatment effect
##' @export
##' @importFrom stats lm predict weights
ittestimate <- function(design,
                        data,
                        outcome,
                        target = "ate",
                        cov_adj_model = NULL,
                        by = NULL,
                        ...,
                        weights = NULL) {

  if (!class(design) %in% c("Design", "WeightedDesign")) {
    stop("design must be Design or WeightedDesign")
  }
  if (!is.data.frame(data)) {
    stop("data must be data.frame")
  }
  if (!is.character(outcome)) {
    stop("outcome must be quoted name of outcome variable in `data`")
  }

  if (is(design, "WeightedDesign")) {
    weights <- design
    design <- design@Design
  }

  if (!outcome %in% colnames(data)) {
    stop("outcome must be a column in data")
  }

  if (!is.null(cov_adj_model)) {
    cov_adj <- cov_adj(cov_adj_model, design = design)
  }

  if (!is.null(by)) {
    design <- .update_by(design, data, by)
  }

  if ("z__" %in% colnames(data)) {
    stop(paste("'z__' is used internally to identify treatment; please",
               "rename column in your data"))
  }

  # Expand treatment status
  ctdata <- design@structure[, design@column_index %in% c("t", "u")]
  colnames(ctdata)[1] <- "z__"
  if (is.logical(ctdata$z__)) {
    # To avoid `z__` being renamed below, convert logical to numeric
    ctdata$z__ <- as.numeric(ctdata$z__)
  }
  merged <- .merge_preserve_order(data, ctdata, by = colnames(ctdata)[-1])

  if (is.null(weights)) {
    wtfn <- switch(target,
                   "ate" = ate,
                   "ett" = ett,
                   stop("invalid target"))
    weights <- wtfn(design, data = merged)
  } else {
    if (length(weights) != nrow(merged)) {
      stop("weights must be same length as rows of `data`")
    }
    if (!is.numeric(weights)) {
      stop("weights must be numeric")
    }
  }


  if (!is.null(cov_adj_model)) {
    model <- lm(merged[, outcome] ~ z__ + offset(cov_adj),
                data = merged, weights = weights, ...)
  } else {
    model <- lm(merged[, outcome] ~ z__,
                data = merged, weights = weights, ...)
  }

  return(new("DirectAdjusted", model, Design = design, target = target))
}
