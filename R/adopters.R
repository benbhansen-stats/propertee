##' @title Adopters
##' @param design Optional design
##' @return Treatment
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
                        })

  return(treatment)
}
