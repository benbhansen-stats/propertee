# Supported weights (targets) refer to the actual calculations. All targets
# should have at minimum one alias, themselves.

##' @title Valid Weights
##'
##' @description These track and access the weights we support as well as any aliases.
##'
##' @details "target" refers to weight calculations we support.
##'
##' "alias" refers to all possible names. Every "target" has at least one
##' "alias" (itself) and may have more.
##'
##' \code{.isValidWeightTarget()} and \code{.isValidWeightAlias()} identify
##' whether a given input (from a user) is a value weighting name.
##'
##' \code{.listValidWeightTargets()} and \code{.listValidWeightAliases()} are
##' for returning nicely formatted strings for messages to users.
##'
##' IMPORTANT: Adding new aliases MUST correspond to new functions defined in
##' weights_exported.R, with possible adjustments to calculations in
##' weights_internal.R.
##' @param target String
##' @param alias String
##' @return Logical for \code{.isValid*}, and a string for \code{.listValid*}.
##' @keywords internal
##' @rdname validWeights
.validWeights <- list(
  targets = c("ate", "ett", "etc", "ato"),
  aliases = c("ate", "ett", "att", "etc", "atc", "ato", "olw",
              "owt", "pwt")
)

##' @keywords internal
##' @rdname validWeights
.isValidWeightTarget <- function(target) {
  stopifnot(is.character(target))
  return(target %in% .validWeights$targets)
}

##' @keywords internal
##' @rdname validWeights
.isValidWeightAlias <- function(alias) {
  stopifnot(is.character(alias))
  return(alias %in% .validWeights$aliases)
}

##' @keywords internal
##' @rdname validWeights
.listValidWeightTargets <- function() {
  return(paste(.validWeights$targets, collapse = ", "))
}


##' @keywords internal
##' @rdname validWeights
.listValidWeightAliases <- function() {
  return(paste(.validWeights$aliases, collapse = ", "))
}
