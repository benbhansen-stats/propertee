#' @include validWeights.R
NULL


##' Called from \code{ate()} or \code{ett()}.
##' @title (Internal) Worker function for weight calculation
##' @param specification a \code{StudySpecification} object created by one of
##'   \code{rct_spec()}, \code{rd_spec()}, or \code{obs_spec()}.
##' @param target One of "ate" or "ett"; \code{ate()} and \code{ett()} chooses
##'   these automatically.
##' @param weightAlias An alias for the weight target, currently one of "ate",
##'   "ett", "att". Chosen by \code{ate()} and \code{ett()} automatically.
##' @param dichotomy optional; a formula defining the dichotomy of the treatment
##'   variable if it isn't already \code{0}/\code{1}. See details of help for
##'   \code{ate()} or \code{ett()} e.g. for details.
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
.weights_calc <- function(specification,
                          target,
                          weightAlias,
                          dichotomy,
                          by,
                          data) {
  if (!(.isValidWeightTarget(target))) {
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
  
  # get the units of assignment from the StudySpecification to calculate weights
  txt <- .bin_txt(specification, dichotomy = dichotomy)
  
  #### generate weights
  
  if (!has_blocks(specification)) {
    # If no block is specified, then e_z is the proportion of units of
    # assignment who receive the treatment.
    e_z <- mean(txt, na.rm = TRUE)
  } else {
    # If a block is specified, then e_z varies by block and is the proportion
    # of units of assignments within the block that receive the treatment.
    blks <- specification@structure[, c(var_names(specification, "u"),
                                        var_names(specification, "b"))]
    
    # Generate e_z and merge back into blks
    e_z_ <- aggregate(txt ~ .,
                      data = blks[,var_names(specification, "b"), drop=FALSE],
                      FUN=function(x) sum(x)/length(x))
    blks <- .merge_preserve_order(blks, e_z_, all.x = TRUE)
    # Rename for clarity
    names(blks)[ncol(blks)] <- "e_z"

    e_z <- blks[["e_z"]]
    
  }
  
  weights  <-
    switch(target,
           ate = txt / e_z + (1 - txt) / (1 - e_z),
           ett = txt + (1 - txt) * e_z / (1 - e_z),
           etc= (1-txt) + txt * (1-e_z)/e_z,
           ato= txt*(1-e_z) + (1-txt)*e_z
    )
  
  to_reset_to_0 <- e_z == 1 | e_z == 0
  if (any(to_reset_to_0, na.rm = TRUE)) {
    weights[to_reset_to_0] <- 0
  }
  
  return(.join_spec_weights(weights,
                            specification,
                            target = target,
                            weightAlias = weightAlias,
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
##' @param weightAlias Any currently supported alias
##' @param data New data
##' @param dichotomy formula used to specify a dichotomy of a non-binary
##'   treatment variable. The output \code{WeightedStudySpecification} object
##'   will store this as its \code{dichotomy} slot, unless it is NULL, in which
##'   case it will be translated to an empty \code{formula}.
##' @return a \code{WeightedStudySpecification}
##' @keywords internal
.join_spec_weights <- function(weights,
                               specification,
                               target,
                               weightAlias,
                               data,
                               dichotomy) {
  uoanames <- var_names(specification, "u")
  
  # Merge uoa data with weights at uoa level
  uoadata <- specification@structure[, var_names(specification, "u"), drop = FALSE]
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
             weightAlias = weightAlias,
             dichotomy = as.formula(dichotomy)))
}
