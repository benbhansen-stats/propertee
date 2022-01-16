# Internal function to try and retrieve the data from the model when `ate` or
# `ett` are called without a data argument
.get_data_from_model <- function() {
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
    stop(paste0("No supported model type found. Supported models: ",
                paste(modelsSupported, sep = ", ")))
  }
  ###*** The model call won't be in the first frame if user calls something like
  # `summary(lm(...`. Also, doing it this way allows us to easily change from
  # looking for a supported model, to looking for a call to `model.frame`.

  modelcall <- fnsCalled[modelFramePos]

  ##########################
  ##### Model type: lm #####
  ##########################
  if (modelcall == "lm") {
   data <- get("data", envir = sys.frame(modelFramePos))
  }

  ###################################
  ##### Model type: ittestimate #####
  ###################################
  if (modelcall == "ittestimate") {
   data <- get("data", envir = sys.frame(modelFramePos))
  }

  # Ensure data is actually data at this point
  if (is.null(data) || !is.data.frame(data)) {
    stop("Could not determine appropriate data")
  }
  return(data)
}
