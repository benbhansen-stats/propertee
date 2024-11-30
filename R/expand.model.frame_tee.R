##' A variation of expand.model.frame which works for \code{teeMod}
##' objects
##'
##' When building a \code{teeMod} object inside \code{lmitt()}, we do a lot of
##' manipulation of the variables involved in the model such that by the time
##' the \code{teeMod} is produced, neither the outcome nor predictors actually
##' fit in the model exist in the \code{data} passed into the call.
##'
##' (E.g. to be specific, if a user calls \code{myda <- lmitt(y ~ 1, data =
##' mydata)}, then \code{model.frame(myda)} would contain column names not found
##' in \code{mydata}.)
##'
##' This is a clone of \code{stats::expand.model.frame()} which has one addition
##' - after extracting the \code{model$call$data} from \code{model}, it adds
##' columns from \code{model.frame(model)} to the object. This ensures that the
##' additional variables created during \code{lmitt()} can be found.
##'
##' Trivial modifications from \code{stats::expand.model.frame()} include
##' ensuring \code{model} is a \code{teeMod} object, and using the \code{::}
##' syntax as appropriate.
##' @title Add new variables to a model frame from a \code{teeMod} object
##' @param model A \code{teeMod} object
##' @param extras one-sided formula or vector of character strings describing
##'   new variables to be added
##' @param envir an environment to evaluate things in
##' @param na.expand logical; see \code{stats::expand.model.frame} for details
##' @return A \code{data.frame}
##' @keywords internal
.expand.model.frame_teeMod <- function (model,
                                    extras,
                                    envir = environment(formula(model)),
                                    na.expand = FALSE) {
  # R4.2.3 or earlier
  if (as.numeric(version$major) < 4 |
        (as.numeric(version$major) == 4 & as.numeric(version$minor) < 3)) {
    stopifnot(is(model, "teeMod")) # JE addition
    f <- stats::formula(model) # JE modification
    data <- eval(model$call$data, envir)
    data <- cbind(data, stats::model.frame(model)) # JE addition
    ff <- foo ~ bar + baz
    gg <- if (is.call(extras))
            extras
    else str2lang(paste("~", paste(extras, collapse = "+")))
    ff[[2L]] <- f[[2L]]
    ff[[3L]][[2L]] <- f[[3L]]
    ff[[3L]][[3L]] <- gg[[2L]]
    if (!na.expand) {
      naa <- model$call$na.action
      subset <- model$call$subset
      rval <- eval(call("model.frame", ff, data = data, subset = subset,
                        na.action = naa), envir)
    }
    else {
      subset <- model$call$subset
      rval <- eval(call("model.frame", ff, data = data, subset = subset,
                        na.action = I), envir)
      oldmf <- stats::model.frame(model) # JE modification
      keep <- match(rownames(oldmf), rownames(rval))
      rval <- rval[keep, ]
      class(rval) <- "data.frame"
    }
    return(rval)
  } else {
    stopifnot(is(model, "teeMod")) # JE addition
    f <- stats::formula(model) # JE modification
    cl <- getCall(model)
    data <- cl$data
    data <- cbind(data, stats::model.frame(model)) # JE addition
    f <- formula(model)
    ff <- foo ~ bar + baz
    gg <- if (is.call(extras))
        extras
    else str2lang(paste("~", paste(extras, collapse = "+")))
    ff[[2L]] <- f[[2L]]
    ff[[3L]][[2L]] <- f[[3L]]
    ff[[3L]][[3L]] <- gg[[2L]]
    environment(ff) <- envir
    if (!na.expand) {
        rval <- eval(call("model.frame", ff, data = data,
            subset = cl$subset, na.action = cl$na.action), envir)
    }
    else {
        rval <- eval(call("model.frame", ff, data = data,
            subset = cl$subset, na.action = I), envir)
        oldmf <- stats::model.frame(model) # JE modification
        keep <- match(rownames(oldmf), rownames(rval))
        rval <- rval[keep, ]
        class(rval) <- "data.frame"
    }
    return(rval)
}

}
