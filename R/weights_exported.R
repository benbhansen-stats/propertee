#' @include weights_internal.R
NULL

############## ***IMPORTANT*** ######################
# If adding new weights here (either direct targets or aliases, YOU MUST UPDATE
# the `.validWeights` object inside validWeights.R

##' @title Generate Direct Adjusted Weights for Treatment Effect Estimation
##'
##' @description These should primarily be used inside models. See Details.
##'
##' @details These functions should primarily be used in the \code{weight}
##'   argument of [lmitt()] or[lm()]. All arguments are optional if used within
##'   those functions. If used on their own, \code{specification} and
##'   \code{data} must be provided.
##'
##'   - \code{ate} - Average treatment effect. Aliases: \code{ate()}.
##'   - \code{ett} - Effect of treatment on the treated. Aliases: \code{ett()},
##'                        \code{att()}.
##'   - \code{etc} - Effect of treatment on controls. Aliases: \code{etc()},
##'                        \code{atc()}.
##'   - \code{ato} - Overlap-weighted average effect. Aliases: \code{ato()},
##'                        \code{olw}, \code{owt}, \code{pwt}.
##'
##'   In a \code{StudySpecification} with \code{block}s, the weights are
##'   generated as a function of the ratio of the number of treated units in a
##'   block versus the total number of units in a block.
##'
##'   In any blocks where that ratio is 0 or 1 (that is, all units in the block
##'   have the same treatment status), the weights will be 0. In effect this
##'   removes from the target population any block in which there is no basis
##'   for estimating either means under treatment or means under control.
##'
##'   If block is missing for a given observation, a weight of 0 is applied.
##'
##'   A \code{dichotomy} is specified by a \code{formula} consisting of a
##'   conditional statement on both the left-hand side (identifying treatment
##'   levels associated with "treatment") and the right hand side (identifying
##'   treatment levels associated with "control"). For example, if your
##'   treatment variable was called \code{dose} and doses above 250 are
##'   considered treatment, you might write:
##'
##'   \code{ate(..., dichotomy = dose > 250 ~ dose <= 250}
##'
##'   The period (\code{.}) can be used to assign all other units of assignment.
##'   For example, we could have written the same treatment regime as either
##'
##'   \code{etc(..., dichotomy = dose > 250 ~ .}
##'
##'   or
##'
##'   \code{olw(..., dichotomy = . ~ dose <= 250}
##'
##'   The \code{dichotomy} formula supports Relational Operators (see
##'   [Comparison]), Logical Operators (see [Logic]), and \code{%in%} (see
##'   [match()]).
##'
##'   The conditionals need not assign all values of treatment to control or
##'   treatment, for example, \code{dose > 300 ~ dose < 200} does not assign
##'   \code{200 <= dose <= 300} to either treatment or control. This would be
##'   equivalent to manually generating a binary variable with \code{NA}
##'   whenever \code{dose} is between 200 and 300. Standard errors will reflect
##'   the sizes of the comparison groups specified by the \code{dichotomy}.
##'
##'   Tim Lycurgus contributed code for the computation of weights. The
##'   \sQuote{overlap weight} concept is due to Li, Morgan and Zaslavsky
##'   (2018), although the current implementation differs from that
##'   discussed in their paper in that it avoids estimated propensity scores.
##'
##' @param specification optional; a \code{StudySpecification} object created by
##'   one of \code{rct_spec()}, \code{rd_spec()}, or \code{obs_spec()}.
##' @param dichotomy optional; a formula defining the dichotomy of the treatment
##'   variable if it isn't already \code{0}/\code{1}. See details for more
##'   information. If \code{ett()} or \code{ate()} is called within a
##'   \code{lmitt()} call that specifies a \code{dichotomy} argument, that
##'   \code{dichotomy} will be used if the argument here has not been specified.
##' @param by optional; named vector or list connecting names of unit of
##'   assignment/ variables in \code{specification} to unit of
##'   assignment/unitid/cluster variables in \code{data}. Names represent
##'   variables in the StudySpecification; values represent variables in the
##'   data. Only needed if variable names differ.
##' @param data optional; the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the \code{weights} argument of
##'   a model.
##' @return a \code{WeightedStudySpecification} object, which is a vector of
##'   numeric weights
##' @export
##' @rdname WeightCreators
##' @references Li, Fan, Kari Lock Morgan, and Alan M. Zaslavsky.
##' "Balancing covariates via propensity score weighting." Journal of the American Statistical Association 113, no. 521 (2018): 390-400.
##' @examples
##' data(simdata)
##' spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
##' summary(lmitt(y ~ 1, data = simdata, specification = spec, weights = ate()), vcov.type = "CR0")
ett <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ett",
                       weightAlias = "ett",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @rdname WeightCreators
##' @export
att <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ett",
                       weightAlias = "att",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @rdname WeightCreators
##' @export
ate <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ate",
                       weightAlias = "ate",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @rdname WeightCreators
##' @export
etc <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "etc",
                       weightAlias = "etc",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @rdname WeightCreators
##' @export
atc <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "etc",
                       weightAlias = "atc",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}


##' @rdname WeightCreators
##' @export
ato <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ato",
                       weightAlias = "ato",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @rdname WeightCreators
##' @export
olw <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ato",
                       weightAlias = "olw",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @rdname WeightCreators
##' @export
owt <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ato",
                       weightAlias = "owt",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}

##' @rdname WeightCreators
##' @export
pwt <- function(specification = NULL, dichotomy = NULL, by = NULL, data = NULL) {
  return(.weights_calc(specification = specification,
                       target = "ato",
                       weightAlias = "pwt",
                       dichotomy = dichotomy,
                       by = by,
                       data = data))
}
