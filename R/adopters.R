adopters <- function(design = NULL) {
  if (is.null(design)) {
    # Identify all frames with a weights argument
    weights_args <- lapply(sys.calls(), `[[`, "weights")
    # Loop over each frame which has a weight argument.
    # Its most likely the first frame, but perhaps not.
    for (i in which(!vapply(weights_args, is.null, logical(1)))) {
      possible_design <- get("weights", sys.frame(i))
      if (is(possible_design, "WeightedDesign")) {
        # If we have a WeightedDesign, save it and break
        wd <- possible_design
        break()
      }
    }
  }

  data <- .get_data_from_model_2(wd@Design@call$formula)


  treatment_uoa <- cbind(treatment(wd@Design),
                         wd@Design@structure[, var_names(wd@Design, "u"),
                                             drop = FALSE])

  treatment_data <- .merge_preserve_order(data, treatment_uoa,
                                          by = var_names(wd@Design, "u"))

  treatment <- tryCatch(treatment_data[, var_names(wd@Design, "t")],
                        error = function(e) {
                          treatment_data[, paste0(var_names(wd@Design, "t"), ".y")]
                        })

  return(treatment)
}



# Internal function to try and retrieve the data from the model when `ate` or
# `ett` are called without a data argument
.get_data_from_model_2 <- function(form = NULL, search_function = NULL, by = NULL) {
  design <- get("design", sys.frame(-1))


  # Below we try to be intelligent about finding the appropriate data. If this
  # fails, we may have need for a brute force method that just loops through
  # frames and looks for a `data` object.
  .fallback <- function() {
    for (i in seq_len(sys.nframe())) {
      try(data <- get("data", envir = sys.frame(i)),
          silent = TRUE)
      if (!is.null(data) && is.data.frame(data)) {
        break()
      }
    }
    if (is.null(data) || !is.data.frame(data)) {
      stop("Could not determine appropriate data")
    }
    return(data)
  }

  data <- NULL

  # Obtain the names of all functions in the callstack
  fns_called <- as.character(lapply(sys.calls(), `[[`, 1))

  # Identify all frames which have `model.frame.default` called (or `ittestimate`)
  model_frame_pos <- which(fns_called %in% c("model.frame.default", "ittestimate"))

  # Identify the frames in which `ate` or `ett` are passed as weights
  adopter_args <- lapply(sys.calls(), `[[`, "formula")
  adopters_pos <- which(grepl("adopters([a-zA-Z0-9]*)", lapply(adopter_args, deparse)))

  # At this point, we know all frames in which the weights are passed,
  # and all frames which are calls to `model.frame.default`. This will identify
  # which frames are both.
  model_frame_and_adopters_pos <- intersect(model_frame_pos, adopters_pos)

  if (length(model_frame_pos) < 1) {
    warning(paste0("No call to `model.frame` with Design weights in the call stack found.\n",
                   "Trying fallback method to obtain data"))
    try(data <- .fallback(),
        silent = TRUE)
  } else if (length(model_frame_and_adopters_pos) > 1) {
    stop(paste0("Multiple models with weights found on the call stack.\n",
                "Try splitting the call into multiple lines rather than a single nested model."))
  } else {

    environment(form) <- environment(get("formula", sys.frame(model_frame_and_adopters_pos)))
    try(data <- do.call(data.frame,
                        c(model.frame(form,
                                      data = get("data",
                                                 sys.frame(model_frame_and_weights_pos)),
                                      na.action = na.pass),
                                     check.names = FALSE)),
        silent = TRUE)
  }

  # Ensure data is actually data at this point - if not, try fallback method
  if (is.null(data) || !is.data.frame(data)) {
    warning("Unable to detect data by normal means, trying fallback method to obtain data")
    data <- .fallback()
  }

  # Finally if even fallback fails, exit
  if (is.null(data) || !is.data.frame(data)) {
    stop("Could not determine appropriate data")
  }

  data <- .rename_model_frame_columns(data)[["renamedModelFrame"]]

  return(data)
}
