##' When included as a predictor in a model formula where weights (via
##' \code{ate()} or \code{ett()} being passed to the \code{weights=} argument)
##' or a covariance adjustment model (via \code{cov_adj()} being passed to the
##' \code{offset=} argument) are included, this will extract the treatment
##' variable from the \code{Design}.
##'
##' While \code{adopters()} can be used in any situation, it is most useful for
##' scenarios where the treatment variable is non-binary and the \code{Design}
##' contains a \code{Dichotomy}. For example, say \code{q} is a 3-level ordinal
##' treatment variable, and the binary comparison of interest is captured in
##' \code{dichotomy = q == 3 ~ q < 3}. If you were to fit a model including
##' \code{q} as a predictor, e.g. \code{lm(y ~ q, ...)}, \code{lm} would treat
##' \code{q} as the full ordinal variable. On the other hand, by calling
##' \code{lm(y ~ adopters(), weights = ate(des), ...)}, \code{adopters()} will
##' generate the appropriate binary variable to allow estimation of treatment
##' effects.
##'
##' @title Obtain Treatment from Design
##' @param design Optional \code{Design}. If the \code{Design} can't be
##'   identified in the model (usually because neither weights (\code{ate()} or
##'   \code{ett()}) nor a covariate adjustment model (\code{cov_adj()}) are
##'   found), the \code{Design} can be passed diretly.
##' @return The treatment variable to be placed in the regression formula.
##' @export
adopters <- function(design = NULL) {
  if (is.null(design)) {
    design <- .get_design()
  }

  # Only thing we need from the data is cluster info to enable later merge
  form <- as.formula(paste("~", paste(var_names(design, "u"),
                                      collapse = "+")))

  data <- .get_data_from_model("adopters", form)


  # Extract treatment and unitofassignment variables from the Design
  treatment_uoa <- cbind(treatment(design),
                         design@structure[, var_names(design, "u"),
                                             drop = FALSE])

  # Merge the extracted pieces from the Design with the data being used to build
  # the model.
  treatment_data <- .merge_preserve_order(data, treatment_uoa,
                                          by = var_names(design, "u"))

  treatment <- tryCatch(treatment_data[, var_names(design, "t")],
                        error = function(e) {
                          # if treatment variable already exists in data, there
                          # will be a .x and .y version; e.g. z.x and z.y, so
                          # we'll extract the ".y" version (the second one)
                          # since the merge above has the treatment from the
                          # Design second.
                          treatment_data[, paste0(var_names(design, "t"), ".y")]
                        })

  return(treatment)
}
