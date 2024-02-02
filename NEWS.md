# v0.3.5

## Improvements
* Diagonal elements of `vcovDA()` matrices lacking sufficient degrees of freedom for estimation are returned as NA's rather than numeric zeros. This is a deviation from the `sandwich` package that aims to provide clarity to results that may otherwise appear as negative diagonal elements of the vcov matrix

# v0.3.4

## Improvements
* Computational performance for `estfun.DirectAdjusted` has been improved

## Bug Fixes
* No more errors due to under-the-hood duplication of a moderator variable
* `absorb=TRUE` estimates have been corrected in the case when all observations in a stratum have 0 weights due to only treated or control units of assignment existing in the stratum

# v0.3.3

## Added Features
* `vcovDA()` can accept user-created variance estimation functions that start with the prefix `.vcov_`; the `type` argument should take the rest of the function name as an input
* Variance estimation for robust GLM's (models fit using `robustbase::glmrob`) is now accommodated
* HC1 variance estimates are now accommodated

# v0.3.2

## Added Features
* Effect estimation for continuous moderator variables is now supported

## Non-Breaking Changes
* `vcovDA()` will return NA's for the entries of the covariance matrix that lack sufficient degrees of freedom for an estimate. Informative warnings will accompany the matrix, further indicating which standard errors have been NA'd out.

## Bug Fixes
* Functions for generating weights, `ate()` and `ett()`, return weights of 0 rather than infinity for blocks that contain treated units but no control units.
* Prior covariate adjustment fits were previously incorporated into variance estimation differently depending on whether one created a `SandwichLayer` object before calling `lmitt()` or called `cov_adj()` in the `offset` argument of the `lmitt()` call. This has been corrected, and both ways return the same variance estimates.
* Covariate adjustment models that admit rectangular bread matrices, such as those
produced by `robustbase::lmrob`, are now accommodated given the reformulated estimating
equations in versions `v0.1.1` and later.
* A contrasts error raised by `model.matrix()` in certain `cov_adj()` calls has been resolved.

# v0.3.1

## Breaking Changes
* We now order `DirectAdjusted` objects' matrix of estimating equations based on user-specified ID columns or unit of assignment ID's.
* The `stats::update` function can no longer be called on `DirectAdjusted` objects.

## Non-Breaking Changes
* `DirectAdjusted` objects now have `lmitt_call` slots.
* `summary` calls on `DirectAdjusted` objects accept `vcov.type` arguments to specify the desired standard error calculation shown in the output. Acceptable types follow the documentation for `vcovDA`.
* Shown or printed `DirectAdjusted` objects return more comprehensible labels for ITT effect outputs.

## R Version Compatibility
* Now compatible with R 4.3. Particularly, we advise users working with R 4.3 to avoid `expand.model.frame` calls on `DirectAdjusted` objects and instead use the internal function `.expand.model.frame.DA` when necessary.

# v0.2.1

## Breaking Changes
* Stratum fixed effects and subgroup moderating effects can now be accounted for via the `absorb` argument. Previous versions did not properly support this functionality. Valid standard errors under absorption, however, have not been confirmed.

# v0.1.1

## Breaking Changes
* We have reformulated the estimating equations used to derive standard errors. In estimation settings we accommodate, testing has not revealed any differences in standard error estimates between the previous and current estimating equations, but we do not assure this is the case for all possible situations.

# v0.0.1

* Compatible with R 4.2.3
* Introduces functionality for direct adjusted and design-informed standard errors accommodating covariance adjustment in the model-based setting
* Cluster-robust standard errors can only be estimated using the HC0 estimator
