# Internal function to try and retrieve the data from the model when `ate` or
# `ett` are called without a data argument
.get_data_from_model <- function(form, unitOfAssignmentIds = NULL) {

  if (!(is.call(form) | is(form, "formula") | is.name(form))) {
    stop("internal error: form must be a formula or name")
  }

  # Evaluate as needed.
  if (is.name(form)) {
    form <- dynGet(form) # dynGet searches iteratively backwards in the
                         # callstack. Per the documentation, it is experimental
                         # and should be used with caution. We should only be hitting
                         # this if the user passes a predefined formula, e.g.:
                         # f <- y ~ x, RCT_Design(f)
  } else if (is.call(form)) {
    form <- as.formula(form)
  }
  # By this point we should have hit an earlier error
  if (!is(form, "formula")) {
    stop("internal error: unable to convert form to formula class")
  }

  if (!is.null(unitOfAssignmentIds)) {
    # if we're passed unitOfAssignmentIds, need to update the formula.

    # First, convert passed-as-string new variable names to `name` objects
    # to avoid quotation (y ~ "x")` in formula
    unitOfAssignmentIds <- lapply(unitOfAssignmentIds, as.name)

    # Next use substitute to do the replacing.
    # Cannot just call `substitute(form, unitOfAssignmentIds)` since we want to pass
    # the actual formula object, not just the name "form"
    form <- as.formula(do.call("substitute", list(form, unitOfAssignmentIds)))
  }

  # update formula to always use unitOfAssignment, since if this is the original
  # call to *_Design, user may have used cluster/uoa/unitid
  form <- .updateFormToUnitOfAssignment(form)

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
  weightArgs <- lapply(sys.calls(), `[[`, "weights")

  # If the call is `weights = ate(des)`, then we'll be able to find "ate"
  # directly. This corresponds to `weightsPos`. However, if the call is more
  # complex, like `weights = 3*ate(des)` or `weights = sqrt(ate(des))`, then
  # we'll need to find `ate(` instead, leading to `atePos` and `ettPos`.
  weightsPos <- lapply(weightArgs, `[[`, 1) %in% c("ate", "ett")
  atePos <- vapply(lapply(weightArgs, as.character), function(x) any(grepl("^ate\\(", x)), TRUE)
  ettPos <- vapply(lapply(weightArgs, as.character), function(x) any(grepl("^ett\\(", x)), TRUE)

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
  } else if (length(modelFrameandWeightsPos) > 1) {
    stop(paste0("Multiple models with weights found on the call stack.\n",
                "Try splitting the call into multiple lines rather than a single nested model."))
  } else {
    environment(form) <- environment(get("formula", sys.frame(modelFrameandWeightsPos)))
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
