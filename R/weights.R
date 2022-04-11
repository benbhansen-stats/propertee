##' @title Generate Direct Adjusted Weights
##' @param design a Design object created by one of `rct_design`, `rd_design`,
##'   or `obs_design`.
##' @param data optionally the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the `weights` argument of a
##'   model.
##' @param by optional; named vector or list connecting names of cluster/unit of
##'   assignment variables in `design` to cluster/unit of assignment variables
##'   in `data`. Names represent variables in the Design; values represent
##'   variables in the data. Only needed if variable names differ.
##' @param dichotomize optionally, a formula defining the dichotomization of the
##'   treatment variable if it isn't already \code{0}/\code{1}. See details of
##'   help for \code{rct_design()} e.g. for details.
##' @return a WeightedDesign object
##' @export
##' @rdname WeightCreators
ett <- function(design, data = NULL, by = NULL, dichotomize = NULL) {
  .weights_calc(design, data, by, "ett", dichotomize = dichotomize)
}

##' @export
##' @rdname WeightCreators
ate <- function(design, data = NULL, by = NULL, dichotomize = NULL) {
  .weights_calc(design, data, by, "ate", dichotomize = dichotomize)
}

# (Internal) Calculates weights
.weights_calc <- function(design, data, by, target, dichotomize) {
  if (!(target %in% c("ate", "ett"))) {
    stop("Invalid weight target")
  }

  if (is.null(data)) {
    data <- .get_data_from_model(design@call$formula, by)
  }

  if (!is.null(dichotomize)) {
    if (is_dichotomized(design)) {
      warning("design is already dichotomized; over-writing with new `dichotomize`")
    }
    dichotomization(design) <- dichotomize
  }

  if (!is.null(by)) {
    design <- .update_by(design, data, by)
  }

  # Ensure treatment is binary
  if (!has_binary_treatment(design)) {
    stop(paste("To calculate weights, treatment must either be 0/1 binary,\n",
               "or the Design must have a dichotomization."))
  }


  #### generate weights

  txt <- .bin_txt(design)

  if (length(var_names(design, "b")) == 0) {
    # If no block is specified, then e_z is the proportion of clusters who receive
    # the treatment.
    e_z <- mean(txt, na.rm = TRUE)

    if (target == "ate") {
      weights <- txt / e_z + (1 - txt) / (1 - e_z)
    } else if (target == "ett") {
      weights <- txt + (1 - txt) * e_z / (1 - e_z)
    }
  } else {
    # If a block is specified, then e_z varies by block and is the proportion
    # of clusters within the block that receive the treatment.

    # Identify number of units per block, and number of treated units per blodk
    block_units <- table(blocks(design)[!is.na(txt), ])
    block_tx_units <- tapply(txt,
                             blocks(design),
                             FUN = sum,
                             na.rm = TRUE)
    e_z <- block_tx_units / block_units

    # Expand e_z to cluster-level
    e_z <- as.numeric(e_z[blocks(design)[, 1]])

    if (target == "ate") {
      weights <- txt / e_z + (1 - txt) / (1 - e_z)
    } else if (target == "ett") {
      weights <- txt + (1 - txt) * e_z / (1 - e_z)
    }
  }

  .join_design_weights(weights, design, target = target, data = data)
}
