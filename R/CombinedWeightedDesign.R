#' @include WeightedDesign.R
NULL
# The above ensures that `WeightedDesign` is defined prior to
# `CombinedWeightedDesign`

CombinedWeightedDesign <- setClass("CombinedWeightedDesign",
                                   contains = "WeightedDesign",
                                   slots = c(dichotomies = "list",
                                             keys = "list"))

setValidity("CombinedWeightedDesign", function(object) {
  if (length(object@dichotomies) != length(object@keys)) {
    return("discrepancy between length of @dichotomies and @keys")
  }
  TRUE
})

##' @rdname WeightedDesignOps
##' @param force_dichotomy_equal if \code{FALSE} (default), Designs are
##'   considered equivalent even if their \code{dichotomy} differs. If
##'   \code{TRUE}, \code{@dichotomy} must also be equal.
##' @export
##' @importFrom methods slot
##' @importFrom stats formula
setMethod("c", signature(x = "WeightedDesign"),
          function(x, ..., force_dichotomy_equal = FALSE) {
            dots <- list(...)
            # x must be a WeightedDesign to get here; ensure all other elements
            # are as well
            if (any(vapply(dots, function(k) !is(k, "WeightedDesign"), TRUE))) {
              stop("WeightedDesigns can only be combined with other WeightedDesigns")
            }

            # Ensure all WeightedDesigns have the same target
            targets <- c(x@target, lapply(dots, methods::slot, "target"))
            if (length(unique(targets)) > 1) {
              stop(paste("WeightedDesigns can only be concatenated with",
                         "the same target (ate or ett)"))
            }

            designs <- c(x@Design, lapply(dots, methods::slot, "Design"))
            # temporarily remove calls from the Designs to avoid weird discrepancies here
            destmp <- lapply(designs, function(x) {
              x@call <- call("c") # placeholder for unique check
              x
            })

            if (length(unique(destmp)) == 1) {
              # if all Designs are identical (sans @call), we can c together
              # just the numeric weights, and pull out Design/target from the
              # first one
              return(new("WeightedDesign",
                         Reduce(c, lapply(dots, methods::slot, ".Data"), init = x@.Data),
                         Design = x@Design,
                         target = x@target))
            } else {
              if (force_dichotomy_equal) {
                stop("When `force_dichotomy_equal` is `TRUE`, Designs must be identical")
              }
            }

            # If we've made it this far, at least one Design is different from
            # another.

            # Store and remove dichotomys
            dichotomies <- lapply(designs, dichotomy)
            destmp <- lapply(destmp, `dichotomy<-`, stats::formula())

            if (length(unique(destmp)) > 1) {
              # If we're here, we're checking Designs that are forced to have
              # equal calls and dichotomy, so any remaining difference is
              # non-neglible and we cannot proceed.
              stop(paste("Cannot combine WeightedDesigns from Designs which",
                         "differ on elements other than `dichotomy`"))
            }

            # Now that we've made it this far, we've ensured that all designs
            # are identical except dichotomies.

            lengths <- c(length(x), vapply(dots, length, 1))
            keys <- lapply(lengths, seq_len)

            for (i in 2:length(lengths)) {
              keys[[i]] <- keys[[i]] + sum(lengths[1:(i - 1)])
            }

            # The `Reduce` call c's together all the .Data in ... . The `init`
            # argument adds the .Data to it.
            wd <- new("WeightedDesign",
                      Reduce(c, lapply(dots, methods::slot , ".Data"), init = x@.Data),
                      Design = x@Design,
                      target = x@target)
            # Note that the Design is `x@Design` which will have some particular
            # dichotomization. I think its better to force this in general in
            # WeightedDesign (that wd@Design passes has_binary_treatment())
            # since that'll be important in non-`c`'d scenarios. In
            # CombinedWeightedDesigns, we'll have to make sure to ignore that.

            return(new("CombinedWeightedDesign",
                       wd,
                       dichotomies = dichotomies,
                       keys = keys))
          })
