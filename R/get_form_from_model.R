.get_form_from_model <- function(which_fn) {
  form <- NULL

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

    try(form <- formula(get("formula", sys.frame(mf_pos))),
        silent = TRUE)
  }

    # Finally if even fallback fails, exit
  if (is.null(form) || !is(form, "formula")) {
    stop("Could not determine appropriate formula")
  }

  return(form)

}
