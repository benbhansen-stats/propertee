##' @title Generate Direct Adjusted Weights
##' @param design a Design object created by one of \code{rct_design()},
##'   \code{rd_design()}, or \code{obs_design()}.
##' @param dichotomy optionally, a formula defining the dichotomy of the
##'   treatment variable if it isn't already \code{0}/\code{1}. See details of
##'   help for \code{rct_design()} e.g. for details.
##' @param by optional; named vector or list connecting names of cluster/unit of
##'   assignment variables in \code{design} to cluster/unit of assignment
##'   variables in \code{data}. Names represent variables in the Design; values
##'   represent variables in the data. Only needed if variable names differ.
##' @param data optionally the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the \code{weights} argument of
##'   a model.
##' @return a WeightedDesign object
##' @export
##' @rdname WeightCreators
ett <- function(design = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  .weights_calc(design = design,
                target = "ett",
                dichotomy = dichotomy,
                by = by,
                data = data)
}

##' @export
##' @rdname WeightCreators
ate <- function(design = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  .weights_calc(design = design,
                target = "ate",
                dichotomy = dichotomy,
                by = by,
                data = data)
}

# (Internal) Calculates weights
.weights_calc <- function(design, target, dichotomy, by, data) {
  if (!(target %in% c("ate", "ett"))) {
    stop("Invalid weight target")
  }

  if (!is.null(dichotomy)) {
    if (!is(dichotomy, "formula")) {
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
    data <- .get_data_from_model("weights", design@call$formula, by)
  } else if (!is.data.frame(data)) {
    stop("`data` must be `data.frame`")
  }

  if (!is.null(by)) {
    # .update_by handles checking input
    design <- .update_by(design, data, by)
  }

  # Ensure treatment is binary
  if (!has_binary_treatment(design)) {
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
