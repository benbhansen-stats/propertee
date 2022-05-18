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
