##' When included as a predictor in a model formula where weights (via
##' \code{ate()} or \code{ett()} being passed to the \code{weights=} argument)
##' or a covariance adjustment model (via \code{cov_adj()} being passed to the
##' \code{offset=} argument) are included, this will extract the treatment
##' variable from the \code{Design}.
##'
##' While \code{assigned()} can be used in any situation, it is most useful for
##' scenarios where the treatment variable is non-binary and the \code{Design}
##' contains a \code{Dichotomy}. For example, say \code{q} is a 3-level ordinal
##' treatment variable, and the binary comparison of interest is captured in
##' \code{dichotomy = q == 3 ~ q < 3}. If you were to fit a model including
##' \code{q} as a predictor, e.g. \code{lm(y ~ q, ...)}, \code{lm} would treat
##' \code{q} as the full ordinal variable. On the other hand, by calling
##' \code{lm(y ~ assigned(), weights = ate(des), ...)}, \code{assigned()} will
##' generate the appropriate binary variable to allow estimation of treatment
##' effects.
##'
##' If called outside of a model call and without a \code{data} argument, this
##' will extract the treatment from the \code{design}.
##'
##' @title Obtain Treatment from Design
##' @param design Optional \code{Design}. If the \code{Design} can't be
##'   identified in the model (usually because neither weights (\code{ate()} or
##'   \code{ett()}) nor a covariate adjustment model (\code{cov_adj()}) are
##'   found), the \code{Design} can be passed diretly.
##' @param data Optional data set. By default `assigned()` will attempt to
##'   identify the appropriate data, if this fails (or you want to overwrite
##'   it), you can pass the data here.
##' @return The treatment variable to be placed in the regression formula.
##' @rdname AssignedAliases
##' @export
assigned <- function(design = NULL, data = NULL) {
  if (is.null(design)) {
    design <- .get_design()
  }

  # Only thing we need from the data is uoa info to enable later merge
  form <- as.formula(paste("~", paste(var_names(design, "u"),
                                      collapse = "+")))

  if (is.null(data)) {
    suppressWarnings(data <- try(.get_data_from_model("assigned", form),
                                 silent = TRUE))
    if (is(data, "try-error")) {
      warning(paste("`data` cannot be found. Extracting treatment",
                    "from `design` instead."))
      return(treatment(design)[, 1])
    }
  }

  # Extract treatment and unitofassignment variables from the Design
  treatment_uoa <- cbind(treatment(design, binary = "ifany"),
                         design@structure[, var_names(design, "u"),
                                             drop = FALSE])

  # Merge the extracted pieces from the Design with the data being used to build
  # the model.
  treatment_data <- .merge_preserve_order(data, treatment_uoa,
                                          by = var_names(design, "u"),
                                          all.x = TRUE)

  if (is_dichotomized(design)) {
    txtname <- "z__"
  } else {
    txtname <- var_names(design, "t")
  }
  treatment <- tryCatch(treatment_data[, txtname],
                        error = function(e) {
                          # if treatment variable already exists in data, there
                          # will be a .x and .y version; e.g. z.x and z.y, so
                          # we'll extract the ".y" version (the second one)
                          # since the merge above has the treatment from the
                          # Design second.
                          treatment_data[, paste0(txtname, ".y")]
                        })

  return(treatment)
}

##' @rdname AssignedAliases
##' @export
adopters <- assigned

##' @rdname AssignedAliases
##' @export
a. <- assigned

##' @rdname AssignedAliases
##' @export
z. <- assigned
