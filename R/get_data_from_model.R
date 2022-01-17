# Internal function to try and retrieve the data from the model when `ate` or
# `ett` are called without a data argument
.get_data_from_model <- function() {

  # Below we try to be intelligent about finding the appropriate data. If this
  # fails, we may have need for a brute force method that just loops through
  # frames and looks for a `data` object.
  .fallback <- function() {
    for (i in seq_len(sys.nframe())) {
      tryCatch({
        data <- get("data", envir = sys.frame(i))
      })
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

  # A list of supported interactive calls. Probably want to move this somewhere
  # more generic eventually.
  modelsSupported <- c("lm", "ittestimate")

  # Obtain the names of all functions in the callstack
  fnsCalled <- as.character(lapply(sys.calls(), "[[", 1))

  # Identify all frames which contains a function calling a support model
  modelFramePos <- which(fnsCalled %in% modelsSupported)
  if (length(modelFramePos) > 1) {
    stop("Multiple model calls found in call stack")
  }
  if (length(modelFramePos) < 1) {
    warning("No supported models found in stack. Trying fallback method to obtain data")
    return(.fallback())
  }
  ###*** The model call won't be in the first frame if user calls something like
  # `summary(lm(...`. Also, doing it this way allows us to easily change from
  # looking for a supported model, to looking for a call to `model.frame`.

  modelcall <- fnsCalled[modelFramePos]

  ##########################
  ##### Model type: lm #####
  ##########################
  if (modelcall == "lm") {
    try(data <- get("data", envir = sys.frame(modelFramePos)),
        silent = TRUE)
    # if the above line fails, `data` will remain NULL and an informative
    # error will be printed below
  }

  ###################################
  ##### Model type: ittestimate #####
  ###################################
  if (modelcall == "ittestimate") {
    try(data <- get("data", envir = sys.frame(modelFramePos)),
        silent = TRUE)
    # if the above line fails, `data` will remain NULL and an informative
    # error will be printed below
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
  return(data)
}
