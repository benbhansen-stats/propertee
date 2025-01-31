##' Whenever a function in a model
##' (\code{ate()}/\code{ett()}/\code{cov_adj()}/\code{assigned()}) is called
##' without an explicit \code{data=} argument, this will attempt to extract the
##' data from the model itself.
##'
##' The \code{form} specifies what columns of the data are needed. For current
##' use cases (\code{ate()}/\code{ett()} and \code{assigned()}), this will be
##' only the unit of assignment variables, so e.g. \code{form = ~ uoavar}, to
##' enable merging of UOA level variables to the model data. However, this can
##' easily be expanded if other variables are needed.
##' @title (Internal) Locate data in call stack
##' @param which_fn Identify calling function, "weights" or "assigned", helps
##'   separate logic for the two functions.
##' @param form Formula on which to apply \code{model.frame()}. See details
##' @param by translation of unit of assignment/unitid/cluster ID names, passed
##'   down from weights.
##' @return \code{data.frame}
##' @keywords internal
.get_data_from_model <- function(which_fn,
                                 form = NULL,
                                 by = NULL) {

  if (!which_fn %in% c("weights", "cov_adj", "assigned")) {
    stop(paste("Internal error: which_fn is invalid,", which_fn))
  }

  # Formula should be from a StudySpecification; it is used inside `model.frame` below
  if (!(is.call(form) | inherits(form, "formula") | is.name(form))) {
    stop("internal error: form must be a formula or name")
  }

  # Evaluate as needed.
  if (is.name(form)) {
    # dynGet searches iteratively backwards in the
    # callstack. Per the documentation, it is experimental
    # and should be used with caution. We should only be hitting
    # this if the user passes a predefined formula, e.g.:
    # f <- y ~ x; rct_spec(f, ...)
    form <- dynGet(form)
  } else if (is.call(form)) {
    form <- as.formula(form)
  }
  # By this point we should have hit an earlier error
  if (!inherits(form, "formula")) {
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

  data <- NULL

  # Obtain the names of all functions in the callstack
  fns_called <- as.character(lapply(sys.calls(), `[[`, 1))

  # Identify all frames which have `model.frame.default` called (or
  # `lmitt`)
  mf_pos <- which(fns_called %in% c("model.frame.default", "lmitt.formula"))

  # identify whether we're looking inside weights or assigned
  if (which_fn == "weights") {
    # Find all frames with `weights` argument

    weights_args <- lapply(sys.calls(), `[[`, "weights")

    # Within each `weights` argument, we're looking for the string "ate(" or
    # "ett(", ensuring that the previous character (if there is one: ?) is not a
    # letter, thus ensuring we don't have some oddity like `translate(x)`
    # triggering.
    pos <- vapply(weights_args, function(aa) {
      is.call(aa) && (aa[[1]] == quote(ate) || aa[[1]] == quote(ett) ||
                        any(vapply(aa, function(x) {
                         x == quote(ate()) || x == quote(ett())
                        }, logical(1))))
    },   logical(1))

    fn_pos <- which(pos)
  } else if (which_fn == "cov_adj") {
    offset_args <- lapply(sys.calls(), `[[`, "offset")

    pos <- vapply(offset_args, function(aa) {
      is.call(aa) && (aa[[1]] == quote(cov_adj) ||
                        any(vapply(aa, function(x) {
                         x == quote(cov_adj)
                        }, logical(1))))
    },  logical(1))

    fn_pos <- which(pos)
  } else if (which_fn == "assigned") {
    # Find all frames with `formula`
    adopter_args <- lapply(sys.frames(), function(x)
      tryCatch(get("formula", envir = x, inherits = FALSE),
               error = function(e) NULL))

    # Ensure only terms or formula are found #124
    adopter_args <- lapply(adopter_args, function(x) {
      if (!(is(x, "formula") | is(x, "terms"))) {
        return(NULL)
      } else {
        return(x)
      }
    })

    # Search these formula for "assigned(", ensuring that the previous character
    # is not a letter, in case there are functions like customassigned().
    pos <- vapply(adopter_args, function(x) {
      any(c("assigned", "z.", "a.", "adopters") %in%
            all.names(as.formula(x)))
    },
    logical(1))

    fn_pos <- which(pos)
  }


  if (length(mf_pos) == 0) {
    # If no model.frames were identified
    warning(paste0("No call to `model.frame` with StudySpecification weights in the ",
                   "call stack found."))
  } else {
    # We've identified at least one model.frame.default
    if (length(mf_pos) > 1) {
      # If we have more than one, pick the one that has weights/assigned that is
      # lowest in the call-stack. (Pick lowest to avoid any issues where there
      # is pre-processing done at an in-bteween step.)
      mf_pos <- max(mf_pos[mf_pos %in% fn_pos])
    }

    # Regenerate the model.frame with the appropriate data
    environment(form) <-
      environment(get("formula", sys.frame(mf_pos)))

    # Try to get the data from the appropriate frame
    passed_data <- tryCatch(get("data", sys.frame(mf_pos)))

    ### #193 - There's no harm adding this every time, as if it's not needed,
    ### it'll be dropped in the model.frame call below
    passed_data$..uoa.. <- rownames(passed_data)

    if (is.null(passed_data)) {
      # If we can't get the data directly, try getting it via name
      passed_data <- tryCatch(eval(sys.call(mf_pos)$data,
                                   sys.frame(mf_pos)),
                              error = function(e) NULL)
    }
    try(data <- model.frame(form,
                            data = passed_data,
                            na.action = na.pass),
        silent = TRUE)
  }

  # search for a teeMod object in the call stack
  if (is.null(data) || !is.data.frame(data)) {
    for (i in seq_len(sys.nframe())) {
      # recursion stops at the rlang_trace_top_env option if it's been set (i.e.
      # in test_that), otherwise where the call stack ends
      if (identical(parent.frame(i), getOption("rlang_trace_top_env")) |
          is.null(sys.call(-i))) break

      objs <- sapply(ls(parent.frame(i)), get, envir = parent.frame(i))
      lmitts <- vapply(objs, inherits, logical(1), "teeMod")
      if (any(lmitts)) {
        found_lmitt <- objs[[which(lmitts)]]
        try(data <- get("data", envir = environment(formula(found_lmitt))),
            silent = TRUE)
        break
      }
    }
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

  return(data)
}


##' We try to be intelligent about finding the appropriate data. If this fails,
##' we may have need for a brute force method that just loops through frames and
##' looks for a `data` object.
##' @title (Internal) Fallback brute force method to locate \code{data} in the
##'   call stack.
##' @return If found, the data.
##' @keywords internal
.fallback_data_search <- function() {
  for (i in seq(1, sys.nframe())) {
    data <- mget("data", envir = sys.frame(i), ifnotfound = list(NULL))$data
    if (!is.null(data) && is.data.frame(data)) {
      break()
    }
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("Could not determine appropriate data")
  }
  return(data)
}
