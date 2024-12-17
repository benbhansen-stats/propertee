##' @title Generate Direct Adjusted Weights for Treatment Effect Estimation
##'
##' @description These should primarily be used inside models. See Details.
##'   [ate()] creates weights to estimate the Average Treatment Effect and
##'   [ett()] creates weights to estimate Effect of Treatment on the Treated.
##'
##' @details These functions should primarily be used in the \code{weight}
##'   argument of [lmitt()] or[lm()]. All arguments are optional if used within
##'   those functions. If used on their own, \code{specification} and
##'   \code{data} must be provided.
##'
##'   In a \code{StudySpecification} with \code{block}s, the weights are
##'   generated as a function of the ratio of the number of treated units in a
##'   block versus the total number of units in a block.
##'
##'   In any blocks where that ratio is 0 or 1 (that is, all units in the block
##'   have the same treatment status), the weights will be 0. In effect this
##'   removes from the target population any block in which there is no basis
##'   for estimating either means under treatment or means under control.
##'
##'   If block is missing for a given observation, a weight of 0 is applied.
##'
##'   A \code{dichotomy} is specified by a \code{formula} consisting of a
##'   conditional statement on both the left-hand side (identifying treatment
##'   levels associated with "treatment") and the right hand side (identifying
##'   treatment levels associated with "control"). For example, if your
##'   treatment variable was called \code{dose} and doses above 250 are
##'   considered treatment, you might write:
##'
##'   \code{dichotomy(spec) <- dose > 250 ~ dose <= 250}
##'
##'   The period (\code{.}) can be used to assign all other units of assignment.
##'   For example, we could have written the same treatment regime as either
##'
##'   \code{dichotomy(spec) <- dose > 250 ~ .}
##'
##'   or
##'
##'   \code{dichotomy(spec) <- . ~ dose <= 250}
##'
##'   The \code{dichotomy} formula supports Relational Operators (see
##'   [Comparison]), Logical Operators (see [Logic]), and \code{%in%} (see
##'   [match()]).
##'
##'   The conditionals need not assign all values of treatment to control or
##'   treatment, for example, \code{dose > 300 ~ dose < 200} does not assign
##'   \code{200 <= dose <= 300} to either treatment or control. This would be
##'   equivalent to manually generating a binary variable with \code{NA}
##'   whenever \code{dose} is between 200 and 300. Standard errors will reflect
##'   the sizes of the comparison groups specified by the \code{dichotomy}.
##'
##'   Code for the computation of the weights was contributed by Tim Lycurgus.
##'
##' @param specification optional; a \code{StudySpecification} object created by
##'   one of \code{rct_spec()}, \code{rd_spec()}, or \code{obs_spec()}.
##' @param dichotomy optional; a formula defining the dichotomy of the treatment
##'   variable if it isn't already \code{0}/\code{1}. See details for more
##'   information. If \code{ett()} or \code{ate()} is called within a
##'   \code{lmitt()} call that specifies a \code{dichotomy} argument, that
##'   \code{dichotomy} will be used if the argument here has not been specified.
##' @param by optional; named vector or list connecting names of unit of
##'   assignment/ variables in \code{specification} to unit of
##'   assignment/unitid/cluster variables in \code{data}. Names represent
##'   variables in the StudySpecification; values represent variables in the
##'   data. Only needed if variable names differ.
##' @param data optional; the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the \code{weights} argument of
##'   a model.
##' @return a \code{WeightedStudySpecification} object, which is a vector of
##'   numeric weights
##' @export
##' @rdname WeightCreators
##' @examples
##' data(simdata)
##' spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
##' summary(lmitt(y ~ 1, data = simdata, specification = spec, weights = ate()))
ett <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ett",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @export
##' @rdname WeightCreators
ate <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ate",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' Called from \code{ate()} or \code{ett()}.
##' @title (Internal) Worker function for weight calculation
##' @param specification a \code{StudySpecification} object created by one of
##'   \code{rct_spec()}, \code{rd_spec()}, or \code{obs_spec()}.
##' @param dichotomy optional; a formula defining the dichotomy of the treatment
##'   variable if it isn't already \code{0}/\code{1}. See details of help for
##'   \code{ate()} or \code{ett()} e.g. for details.
##' @param target One of "ate" or "ett"; \code{ate()} and \code{ett()} chooses
##'   these automatically.
##' @param by optional; named vector or list connecting names of unit of
##'   assignment/ variables in \code{specification} to unit of
##'   assignment/cluster variables in \code{data}. Names represent variables in
##'   the StudySpecification; values represent variables in the data. Only
##'   needed if variable names differ.
##' @param data optionally the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the \code{weights} argument of
##'   a model.
##' @return a \code{WeightedStudySpecification} object
##' @keywords internal
.weights_calc <- function(specification, target, dichotomy, by, data) {
  if (!(target %in% c("ate", "ett"))) {
    stop("Invalid weight target")
  }

  if (is.null(specification)) {
    specification <- .get_spec()
  }

  # get `dichotomy` argument and validate against any up the call stack
  if (is.null(dichotomy)) {
    possible_dichotomies <- .find_dichotomies()
    dichotomy <- .validate_dichotomy(possible_dichotomies)
  } else {
    dichotomy <- .validate_dichotomy(dichotomy)
  }

  if (is.null(data)) {
    # Only thing we need from the data is unit of assignment info to enable
    # later merge
    form <- as.formula(paste("~", paste(var_names(specification, "u"),
                                        collapse = "+")))

    data <- .get_data_from_model("weights", form, by)
  } else {
    # #174 handle tibbles
    data <- .as_data_frame(data)
  }

  if (!is.null(by)) {
    # .update_by handles checking input
    specification <- .update_by(specification, data, by)
  }

  # getting the unique rows of `data` will ensure we calculate weights at the
  # level of assignment
  txt <- .bin_txt(specification,
                  unique(data[, var_names(specification, "u"), drop = FALSE]),
                  dichotomy)

  #### generate weights

  if (length(var_names(specification, "b")) == 0) {
    # If no block is specified, then e_z is the proportion of units of
    # assignment who receive the treatment.
    e_z <- mean(txt, na.rm = TRUE)

    if (target == "ate") {
      weights <- txt / e_z + (1 - txt) / (1 - e_z)
    } else if (target == "ett") {
      weights <- txt + (1 - txt) * e_z / (1 - e_z)
    }
  } else {
    # If a block is specified, then e_z varies by block and is the proportion
    # of units of assignments within the block that receive the treatment.

    # Identify number of units per block, and number of treated units per block
    # NOTE 5/21/24: since .bin_txt() returns NA elements, need to replace
    # `blocks(specification)` with a matrix that has NA elements at the same indices
    # (and is aligned to the order of the ID's in `data`)
    blks <- as.data.frame(matrix(nrow = length(txt),
                                 ncol = length(var_names(specification, "b")),
                                 dimnames = list(rownames = NULL,
                                                 colnames = var_names(specification, "b"))))

    blks[!is.na(txt),] <- .merge_preserve_order(
      unique(data[, var_names(specification, "u"), drop=FALSE]),
      cbind(specification@structure[, var_names(specification, "b"), drop=FALSE],
            specification@structure[, var_names(specification, "u"), drop=FALSE]),
      all = FALSE, all.x = TRUE
    )[!is.na(txt), var_names(specification, "b")]

    # Generate e_z and merge back into blks
    e_z <- aggregate(txt ~ ., data = blks, FUN=function(x) sum(x)/length(x))
    blks <- .merge_preserve_order(blks, e_z, all.x = TRUE)
    # Rename for clarity
    names(blks)[ncol(blks)] <- "e_z"

    to_reset_to_0 <- blks$e_z == 1 | blks$e_z == 0

    if (target == "ate") {
      weights <- txt / blks$e_z + (1 - txt) / (1 - blks$e_z)
    } else if (target == "ett") {
      weights <- txt + (1 - txt) * blks$e_z / (1 - blks$e_z)
    }

    if (any(to_reset_to_0, na.rm = TRUE)) {
      weights[to_reset_to_0] <- 0
    }

  }

  return(.join_spec_weights(weights,
                            specification,
                            target = target,
                            data = data,
                            dichotomy = dichotomy))
}

##' Helper function called during creation of the weights via \code{ate()} or
##' \code{ett()}
##' @title (Internal) Expand unit of assignment level weights to the level of
##'   the data
##' @param weights a vector of weights sorted according to the
##'   \code{StudySpecification}
##' @param specification a \code{StudySpecification}
##' @param target One of "ate" or "ett"
##' @param data New data
##' @param dichotomy formula used to specify a dichotomy of a non-binary
##'   treatment variable. The output \code{WeightedStudySpecification} object
##'   will store this as its \code{dichotomy} slot, unless it is NULL, in which
##'   case it will be translated to an empty \code{formula}.
##' @return a \code{WeightedStudySpecification}
##' @keywords internal
.join_spec_weights <- function(weights, specification, target, data, dichotomy) {
  uoanames <- var_names(specification, "u")

  # Merge uoa data with weights at uoa level
  # NOTE 5/21/24: changed this from taking `specification@structure` because the order
  # of the weights now reflects the ordering of the data rather than `specification@structure`
  # due to changes in .bin_txt()
  uoadata <- unique(data[, var_names(specification, "u"), drop = FALSE])
  uoadata$specification_weights <- weights

  # Merge with data to expand weights to unit of analysis level
  weights <- .merge_preserve_order(data, uoadata,
                                   by = uoanames,
                                   all.x = TRUE)$specification_weights

  # Replace NA weights with 0 so they don't contribute to the model, but aren't
  # droppde
  weights[is.na(weights)] <- 0

  return(new("WeightedStudySpecification",
             weights,
             StudySpecification = specification,
             target = target,
             dichotomy = as.formula(dichotomy)))
}
