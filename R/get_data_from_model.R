# Internal function to try and retrieve the data from the model when `ate` or
# `ett` are called without a data argument
.get_data_from_model <- function(form) {
  stopifnot(is(as.formula(form), "formula"))

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
  fnsCalled <- as.character(lapply(sys.calls(), `[[`, 1))

  # Identify all frames which have `model.frame.default` called (or `ittestimate`)
  modelFramePos <- which(fnsCalled %in% c("model.frame.default", "ittestimate"))

  # Identify the frames in which `ate` or `ett` are passed as weights
  weightArgs <- sapply(sys.calls(), `[[`, "weights")

  # If the call is `weights = ate(des)`, then we'll be able to find "ate"
  # directly. This corresponds to `weightsPos`. However, if the call is more
  # complex, like `weights = 3*ate(des)` or `weights = sqrt(ate(des))`, then
  # we'll need to find `ate(` instead, leading to `atePos` and `ettPos`.
  weightsPos <- sapply(weightArgs, `[[`, 1) %in% c("ate", "ett")
  atePos <- sapply(lapply(weightArgs, as.character), function(x) any(grepl("^ate\\(", x)))
  ettPos <- sapply(lapply(weightArgs, as.character), function(x) any(grepl("^ett\\(", x)))

  weightsPos <- which(weightsPos | atePos | ettPos)

  # At this point, we know all frames in which the weights are passed,
  # and all frames which are calls to `model.frame.default`. This will identify
  # which frames are both.
  modelFrameandWeightsPos <- intersect(modelFramePos, weightsPos)

  if (length(modelFramePos) < 1) {
    warning(paste0("No call to `model.frame` with Design wights in the call stack found.\n",
                   "Trying fallback method to obtain data"))
    try(data <- .fallback(),
        silent = TRUE)
  } else {
    try(data <- do.call(data.frame,
                        c(model.frame(form,
                                      data = get("data",
                                                 sys.frame(modelFrameandWeightsPos)),
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
