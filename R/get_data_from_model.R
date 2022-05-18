# (Internal) Whenever a function in a model (ate/ett/cov_adj/adopters) is called
# without an explicit `data` argument, this will attempt to extract the data
# from the model itself.
# which_fn = "weights" or "adopters", helps separate logic for the two
# functions.
# form = formula passed to model.frame.
# by = translation of cluster ID names, passed down from weights.
.get_data_from_model <- function(which_fn,
                                 form = NULL,
                                 by = NULL) {

  if (!which_fn %in% c("weights", "adopters")) {
    stop(paste("Internal error: which_fn is invalid,", which_fn))
  }

  # Formula should be from a Design; it is used inside `model.frame` below
  if (!(is.call(form) | is(form, "formula") | is.name(form))) {
    stop("internal error: form must be a formula or name")
  }

  # Evaluate as needed.
  if (is.name(form)) {
    # dynGet searches iteratively backwards in the
    # callstack. Per the documentation, it is experimental
    # and should be used with caution. We should only be hitting
    # this if the user passes a predefined formula, e.g.:
    # f <- y ~ x; rct_design(f, ...)
    form <- dynGet(form)
  } else if (is.call(form)) {
    form <- as.formula(form)
  }
  # By this point we should have hit an earlier error
  if (!is(form, "formula")) {
    stop("internal error: unable to convert form to formula class")
  }

  if (!is.null(by)) {
    .check_by(by)
    # if we're passed by, need to update the formula.

    # First, convert passed-as-string new variable names to `name` objects
    # to avoid quotation (y ~ "x")` in formula
    by <- lapply(by, as.name)

    # Next use substitute to do the replacing.
    # Cannot just call `substitute(form, by)` since we want to pass
    # the actual formula object, not just the name "form"
    form <- as.formula(do.call("substitute", list(form, by)))
  }

  # update formula to always use unit_of_assignment, since if this is the
  # original call to *_Design, user may have used cluster/uoa/unitid
  form <- .update_form_to_unit_of_assignment(form)

  data <- NULL

  # Obtain the names of all functions in the callstack
  fns_called <- as.character(lapply(sys.calls(), `[[`, 1))

  # Identify all frames which have `model.frame.default` called (or
  # `ittestimate`)
  mf_pos <- which(fns_called %in% c("model.frame.default", "ittestimate"))

  # identify whether we're looking inside weights or adopters
  if (which_fn == "weights") {
    # Find all frames with `weights` argument
    weights_args <- lapply(sys.calls(), `[[`, "weights")

    # Within each `weights` argument, we're looking for the string "ate(" or
    # "ett(", ensuring that the previous character (if there is one: ?) is not a
    # letter, thus ensuring we don't have some oddity like `translate(x)`
    # triggering.
    ate_pos <- vapply(lapply(weights_args, deparse),
                      function(x) any(grepl("[^a-zA-Z]?ate\\(", x)), TRUE)
    ett_pos <- vapply(lapply(weights_args, deparse),
                      function(x) any(grepl("[^a-zA-Z]?ett\\(", x)), TRUE)

    fn_pos <- which(ate_pos | ett_pos)
  } else if (which_fn == "adopters") {
    # Find all frames with `formula`
    adopter_args <- lapply(sys.calls(), `[[`, "formula")

    # Search these formula for "adopters(", ensuring that the previous character
    # is not a letter, in case there are functions like customadopters().
    fn_pos <- which(grepl("[^a-zA-Z]?adopters\\(",
                                lapply(adopter_args, deparse)))
  }


  if (length(mf_pos) == 0) {
    # If no model.frames were identified
    warning(paste0("No call to `model.frame` with Design weights in the ",
                   "call stack found."))
  } else {
    # We've identified at least one model.frame.default
    if (length(mf_pos) > 1) {
      # If we have more than one, pick the one that has weights/adopters.
      mf_pos <- mf_pos[mf_pos %in% fn_pos]
      if (length(mf_pos) > 1) {
        # still multiple; too confusing!
        stop("Multiple calls to model.frame.default contain flexida elements")
      }
    }

    # Regenerate the model.frame with the appropriate data
    environment(form) <-
      environment(get("formula", sys.frame(mf_pos)))
    try(data <- do.call(data.frame,
                        c(model.frame(form,
                                      data = get("data", sys.frame(mf_pos)),
                                      na.action = na.pass),
                          check.names = FALSE)),
        silent = TRUE)
  }

  # Ensure data is actually data at this point - if not, try fallback method
  if (is.null(data) || !is.data.frame(data)) {
    warning(paste("Unable to detect data by normal means,",
                  "trying fallback method to obtain data"))
    data <- .fallback_data_search()
  }

  # Finally if even fallback fails, exit
  if (is.null(data) || !is.data.frame(data)) {
    stop("Could not determine appropriate data")
  }

  data <- .rename_model_frame_columns(data)[["renamedModelFrame"]]

  return(data)
}


# (Internal) We try to be intelligent about finding the appropriate data. If
# this fails, we may have need for a brute force method that just loops through
# frames and looks for a `data` object.
.fallback_data_search <- function() {
  for (i in seq_len(sys.nframe())) {
    try({
      data <- get("data", envir = sys.frame(i))
    }, silent = TRUE)
    if (!is.null(data) && is.data.frame(data)) {
      break()
    }
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("Could not determine appropriate data")
  }
  return(data)
}
