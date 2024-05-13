# propertee 0.3.9

## Bug Fixes
* Scaling constants have been updated in `estfun.teeMod()` to account for a previously missing factor of sqrt(n / n_C) applied to contributions to the covariance adjustment model estimating equations

# propertee 0.3.8

## Breaking Changes
* When model-based standard errors clustered at the level of assignment are called for in a blocked design, `vcovDA()` clusters units of assignment in small blocks, blocks with only one treated or control unit, together.

# propertee 0.3.7

## Breaking Changes
* `vcovDA()` scales estimating equations using different constants than it did before

# propertee 0.3.6

## Bug Fixes
* Previous procedure for aligning contributions to estimating equations from first-stage and second-stage models failed when column(s) used for alignment had NA's. Outputs of `vcovDA()` were liable to change from call to call as a result. This has been fixed.

# propertee 0.3.5

## Improvements
* Diagonal elements of `vcovDA()` matrices lacking sufficient degrees of freedom for estimation are returned as NA's rather than numeric zeros. This is a deviation from the `sandwich` package that aims to provide clarity to results that may otherwise appear as negative diagonal elements of the vcov matrix

## Bug Fixes
* When `lmitt()` is called with a blocked design and `absorb=TRUE`, the block-centered assignment and, if applicable, moderator and assignment:moderator interaction columns, are no longer centered on the grand mean of the column. This ensures blocks that do not satisfy positivity of the assignment variable (or positivity within a factor level) do not contribute to effect estimation
* `lmitt()` now accepts references to formula objects

# propertee 0.3.4

## Improvements
* Computational performance for `estfun.teeMod` has been improved

## Bug Fixes
* No more errors due to under-the-hood duplication of a moderator variable
* `absorb=TRUE` estimates have been corrected in the case when all observations in a stratum have 0 weights due to only treated or control units of assignment existing in the stratum

# propertee 0.3.3

## Added Features
* `vcovDA()` can accept user-created variance estimation functions that start with the prefix `.vcov_`; the `type` argument should take the rest of the function name as an input
* Variance estimation for robust GLM's (models fit using `robustbase::glmrob`) is now accommodated
* HC1 variance estimates are now accommodated

# propertee 0.3.2

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

# propertee 0.3.1

## Breaking Changes
* We now order `teeMod` objects' matrix of estimating equations based on user-specified ID columns or unit of assignment ID's.
* The `stats::update` function can no longer be called on `teeMod` objects.

## Non-Breaking Changes
* `teeMod` objects now have `lmitt_call` slots.
* `summary` calls on `teeMod` objects accept `vcov.type` arguments to specify the desired standard error calculation shown in the output. Acceptable types follow the documentation for `vcovDA`.
* Shown or printed `teeMod` objects return more comprehensible labels for ITT effect outputs.

## R Version Compatibility
* Now compatible with R 4.3. Particularly, we advise users working with R 4.3 to avoid `expand.model.frame` calls on `teeMod` objects and instead use the internal function `.expand.model.frame.DA` when necessary.

# propertee 0.2.1

## Breaking Changes
* Stratum fixed effects and subgroup moderating effects can now be accounted for via the `absorb` argument. Previous versions did not properly support this functionality. Valid standard errors under absorption, however, have not been confirmed.

# propertee 0.1.1

## Breaking Changes
* We have reformulated the estimating equations used to derive standard errors. In estimation settings we accommodate, testing has not revealed any differences in standard error estimates between the previous and current estimating equations, but we do not assure this is the case for all possible situations.

# propertee 0.0.1

* Compatible with R 4.2.3
* Introduces functionality for direct adjusted and design-informed standard errors accommodating covariance adjustment in the model-based setting
* Cluster-robust standard errors can only be estimated using the HC0 estimator
