##' @title Obtain Treatment from Design
##' @param design Optional design. If \code{NULL}, the Design will be identified
##'   from its inclusion in a weight (\code{ate()} or \code{ett()}) or from a
##'   covariate adjustment model (\code{cov_adj()}).
##' @return The treatment variable to be placed in the regression formula.
##' @export
adopters <- function(design = NULL) {
  if (is.null(design)) {
    design <- .get_design()
  }

  data <- .get_data_from_model(design@call$formula)


  treatment_uoa <- cbind(treatment(design),
                         design@structure[, var_names(design, "u"),
                                             drop = FALSE])

  treatment_data <- .merge_preserve_order(data, treatment_uoa,
                                          by = var_names(design, "u"))

  treatment <- tryCatch(treatment_data[, var_names(design, "t")],
                        error = function(e) {
                          treatment_data[, paste0(var_names(design, "t"), ".y")]
                          # if treatment variable already exists in data, there
                          # will be a .x and .y version; e.g. z.x and z.y
                        })

  return(treatment)
}
