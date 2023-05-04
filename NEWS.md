# v0.3.2

## Bug Fixes
* Covariate adjustment models that admit rectangular bread matrices, such as those
produced by `robustbase::lmrob`, are now accommodated given the reformulated estimating
equations in versions `v0.1.1` and later

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
