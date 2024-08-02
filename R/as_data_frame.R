##' @title Convert object to \code{data.frame} or produce meaningful error
##' @param x An object
##' @return \code{x} as a \code{data.frame}
##' @keywords internal
.as_data_frame <- function(x) {
  tryCatch( {
    result <- as.data.frame(x)
    return(result)
  },
  error = function(e) {
    msg <- paste("Failed to convert to data frame:",
                 conditionMessage(e))
    stop(msg, call. = FALSE)
  },
  warning = function(w) {
    msg <- paste("Warning occurred while converting to data frame:",
                 conditionMessage(w))
    warning(msg, call. = FALSE)
    return(as.data.frame(x))
  }
  )
}
