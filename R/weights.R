##' @title Generate Direct Adjusted Weights
##' @param design a \code{Design} object created by one of \code{rct_design()},
##'   \code{rd_design()}, or \code{obs_design()}.
##' @param dichotomy optional; a formula defining the dichotomy of the treatment
##'   variable if it isn't already \code{0}/\code{1}. See details of help for
##'   \code{rct_design()} e.g. for details.
##' @param by optional; named vector or list connecting names of cluster/unit of
##'   assignment variables in \code{design} to cluster/unit of assignment
##'   variables in \code{data}. Names represent variables in the Design; values
##'   represent variables in the data. Only needed if variable names differ.
##' @param data optionally the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the \code{weights} argument of
##'   a model.
##' @return a \code{WeightedDesign} object
##' @export
##' @rdname WeightCreators
ett <- function(design = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(design = design,
                       target = "ett",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @export
##' @rdname WeightCreators
ate <- function(design = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(design = design,
                       target = "ate",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' Called from \code{ate()} or \code{ett()}.
##' @title (Internal) Worker function for weight calculation
##' @param design a \code{Design} object created by one of \code{rct_design()},
##'   \code{rd_design()}, or \code{obs_design()}.
##' @param dichotomy optional; a formula defining the dichotomy of the treatment
##'   variable if it isn't already \code{0}/\code{1}. See details of help for
##'   \code{rct_design()} e.g. for details.
##' @param target One of "ate" or "ett"; \code{ate()} and \code{ett()} chooses
##'   these automatically.
##' @param by optional; named vector or list connecting names of cluster/unit of
##'   assignment variables in \code{design} to cluster/unit of assignment
##'   variables in \code{data}. Names represent variables in the Design; values
##'   represent variables in the data. Only needed if variable names differ.
##' @param data optionally the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the \code{weights} argument of
##'   a model.
##' @return a \code{WeightedDesign} object
.weights_calc <- function(design, target, dichotomy, by, data) {
  if (!(target %in% c("ate", "ett"))) {
    stop("Invalid weight target")
  }

  if (!is.null(dichotomy)) {
    if (!inherits(dichotomy, "formula")) {
      stop("`dichotomy` must be a `formula`")
    }
    if (is_dichotomized(design)) {
      warning(paste("design is already dichotomized; over-writing",
                    "with new `dichotomy`"))
    }
    dichotomy(design) <- dichotomy
  }

  if (is.null(design)) {
    design <- .get_design()
  }

  if (is.null(data)) {
    # Only thing we need from the data is cluster info to enable later merge
    form <- as.formula(paste("~", paste(var_names(design, "u"),
                                        collapse = "+")))

    data <- .get_data_from_model("weights", form, by)
  } else if (!is.data.frame(data)) {
    stop("`data` must be `data.frame`")
  }

  if (!is.null(by)) {
    # .update_by handles checking input
    design <- .update_by(design, data, by)
  }

  # Ensure treatment is binary
  if (!is_binary_or_dichotomized(design)) {
    stop(paste("To calculate weights, treatment must either be 0/1 binary,\n",
               "or the Design must have a dichotomy."))
  }


  #### generate weights

  txt <- .bin_txt(design)

  if (length(var_names(design, "b")) == 0) {
    # If no block is specified, then e_z is the proportion of clusters who
    # receive the treatment.
    e_z <- mean(txt, na.rm = TRUE)

    if (target == "ate") {
      weights <- txt / e_z + (1 - txt) / (1 - e_z)
    } else if (target == "ett") {
      weights <- txt + (1 - txt) * e_z / (1 - e_z)
    }
  } else {
    # If a block is specified, then e_z varies by block and is the proportion
    # of clusters within the block that receive the treatment.

    # Identify number of units per block, and number of treated units per block
    block_units <- table(blocks(design)[!is.na(txt), ])
    block_tx_units <- tapply(txt,
                             blocks(design),
                             FUN = sum,
                             na.rm = TRUE)
    e_z <- block_tx_units / block_units

    to_reset_to_0 <- e_z == 1 | e_z == 0
    to_reset_to_0 <- to_reset_to_0[blocks(design)[, 1]]

    # Expand e_z to cluster-level
    e_z <- as.numeric(e_z[blocks(design)[, 1]])

    if (target == "ate") {
      weights <- txt / e_z + (1 - txt) / (1 - e_z)
    } else if (target == "ett") {
      weights <- txt + (1 - txt) * e_z / (1 - e_z)
    }

    if (any(to_reset_to_0, na.rm = TRUE)) {
      weights[to_reset_to_0] <- 0
    }

  }

  return(.join_design_weights(weights, design, target = target, data = data))
}

##' Helper function called during creation of the weights via \code{ate()} or
##' \code{ett()}
##' @title (Internal) Expand unit of assignment level weights to the level of
##'   the data
##' @param weights a vector of weights sorted according to the \code{Design}
##' @param design a \code{Design}
##' @param target One of "ate" or "ett"
##' @param data New data
##' @return a \code{WeightedDesign}
##' @keywords internal
.join_design_weights <- function(weights, design, target, data) {

  uoanames <- var_names(design, "u")

  # Merge uoa data with weights at uoa level
  uoadata <- design@structure[, uoanames, drop = FALSE]
  uoadata$design_weights <- weights

  # Merge with data to expand weights to unit of analysis level
  weights <- .merge_preserve_order(data, uoadata,
                                   by = uoanames,
                                   all.x = TRUE)$design_weights

  return(new("WeightedDesign",
             weights,
             Design = design,
             target = target))
}
