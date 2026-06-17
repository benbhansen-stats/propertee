##' @title Obtain Treatment from StudySpecification
##'
##' @description When passing a \code{lm} object to [lmitt()], extract and use
##'   the treatment variable specified in the \code{StudySpecification}.
##'
##' @details When passing a \code{lm} object to [lmitt()], the treatment
##'   variable in the \code{formula} passed to [lm()] needs to be identifiable.
##'   Rather than placing the treatment variable directly in the \code{formula},
##'   use one of these functions, to allow [lmitt()] to identify the treatment
##'   variable.
##'
##'   To keep the formula in the [lm()] call concise, instead of passing
##'   \code{specification} and \code{data} arguments to these functions, one can
##'   pass a \code{WeightedStudySpecification} object to the \code{weights}
##'   argument of the [lm()] call or a \code{SandwichLayer} object to the
##'   \code{offset} argument.
##'
##'   Alternatively, you can pass the \code{specification} and \code{data}
##'   arguments.
##'
##'   While \code{assigned()} can be used in any situation, it is most useful
##'   for scenarios where the treatment variable is non-binary and the
##'   \code{StudySpecification} contains a \code{Dichotomy}. For example, say
##'   \code{q} is a 3-level ordinal treatment variable, and the binary
##'   comparison of interest is captured in \code{dichotomy = q == 3 ~ q < 3}.
##'   If you were to fit a model including \code{q} as a predictor, e.g.
##'   \code{lm(y ~ q, ...)}, \code{lm} would treat \code{q} as the full ordinal
##'   variable. On the other hand, by calling \code{lm(y ~ assigned(), weights =
##'   ate(spec), ...)}, \code{assigned()} will generate the appropriate binary
##'   variable to allow estimation of treatment effects.
##'
##'   If called outside of a model call and without a \code{data} argument, this
##'   will extract the treatment from the \code{specification}. If this is the
##'   goal, the [treatment()] function is better suited for this purpose.
##'
##' @param specification Optional \code{StudySpecification}. If the
##'   \code{StudySpecification} can't be identified in the model (usually
##'   because neither weights (\code{ate()} or \code{ett()}) nor a covariate
##'   adjustment model (\code{cov_adj()}) are found), the
##'   \code{StudySpecification} can be passed diretly.
##' @param data Optional data set. By default [assigned()] will attempt to
##'   identify the appropriate data, if this fails (or you want to overwrite
##'   it), you can pass the data here.
##' @inheritParams ett
##' @return The treatment variable to be placed in the regression formula.
##' @rdname AssignedAliases
##' @export
##' @examples
##' data(simdata)
##' spec <- obs_spec(z ~ uoa(uoa1, uoa2), data = simdata)
##' mod <- lm(y ~ assigned(), data = simdata, weights = ate(spec))
##' lmittmod <- lmitt(mod)
##' summary(lmittmod, vcov.type = "CR0")
assigned <- function(specification = NULL, data = NULL, dichotomy = NULL) {
  if (is.null(specification)) {
    specification <- .get_spec()
  }

  if (is.null(data)) {
    # Get uoa info to enable merge in .expand_txt, and get non-treatment variables
    # in the dichotomy for .apply_dichotomy
    divars <- setdiff(all.vars(dichotomy), c(var_names(specification, "t"), "."))
    form <- as.formula(paste("~", paste(c(divars, var_names(specification, "u")),
                                        collapse = "+")))
    suppressWarnings(data <- try(.get_data_from_model("assigned", form),
                                 silent = TRUE))
    if (is(data, "try-error")) {
      warning(paste("`data` cannot be found. Extracting treatment",
                    "from `specification` instead."))
      return(treatment(specification)[, 1])
    }
  }

  if (specification@unit_of_assignment_type == "none") {
    data[["..uoa.."]] <- rownames(data)
  }


  tt <- treatment(specification, binary = FALSE)
  tt <- .expand_txt(tt, data, specification)
  if (is.null(dichotomy)) {
    possible_dichotomies <- .find_dichotomies()
    dichotomy <- .validate_dichotomy(possible_dichotomies)
  } else {
    dichotomy <- .validate_dichotomy(dichotomy)
  }


  if (!is.null(dichotomy)) {
    return(.apply_dichotomy(cbind(tt, data), dichotomy))
  } else {
    return(tt[,1])
  }
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
