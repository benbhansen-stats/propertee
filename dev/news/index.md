# Changelog

## **propertee** 1.0.4

- Fix routines for CR2/CV2 standard errors and associated degrees of
  freedom that previously errored when the provided `teeMod` had a
  non-empty `na.action`, or the `teeMod` or its `weights` argument had a
  `dichotomy` slot
- Return NaN for the degrees of freedom when CR2/CV2 standard errors
  have deficient degrees of freedom
- [`as.lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/as_lmitt.md)
  finds `StudySpecification` objects passed to
  [`assigned()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)/[`a.()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)/[`z.()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)/[`adopters()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)
  when the object isn’t explicitly passed to
  [`as.lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/as_lmitt.md)
- Two routines for estimating degrees of freedom are now available: the
  default in
  [`summary.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/teeMod_summary.md)
  is `dof.type=stata`, which uses the number of clusters (or
  observations, in the absence of clustering) less one. We suggest
  specifying `dof.type="IK"` in
  [`summary.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/teeMod_summary.md)
  calls, which uses a routine adapted from Imbens and Kolesár (2016),
  which adapts a routine from Bell and McCaffrey (2002). This routine
  has been available in previous versions for
  [`summary.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/teeMod_summary.md)
  calls with `type=CR2`, but it has been updated so it better matches
  the outputs of the `dfadjust` package maintained by Kolesár. See the
  documentation for the internal function
  [`.compute_IK_dof()`](https://benbhansen-stats.github.io/propertee/dev/reference/dot-compute_IK_dof.md)
  for references and further details.

## **propertee** 1.0.3

CRAN release: 2025-10-28

- Fix bug in
  [`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
  when covariance adjustment model formula includes transformations of
  variables

## **propertee** 1.0.2

- Fix bug in variance estimation code when the mean response/offset in a
  subpopulation of the control condition is not identified

## **propertee** 1.0.1

CRAN release: 2025-08-26

- Fix small test suite bugs

## **propertee** 1.0.0

CRAN release: 2025-08-21

- First CRAN submission

## **propertee** 0.7.0

### New Features

- HC2 and CR2 standard errors are now available for `teeMod` objects.
  CR2 corrections use a new, fast computation that obviates the need for
  obtaining the spectral decomposition of orthocomplements of
  cluster-specific projection matrices.
- Standard errors for `teeMod` objects now optionally include
  methodological developments from a forthcoming paper from Wasserman
  and Hansen:
  1.  When computing the meat of the sandwich standard errors, the bias
      correction to the residuals from estimating intent-to-treat
      effects is allowed to differ from the bias correction to residuals
      from a covariance adjustment model fit. The `type` argument of
      [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
      determines the correction for the former, while the
      `cov_adj_rcorrect` argument determines the latter. The findings
      from the simulation study in the paper inform the default
      arguments of
      [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
      (and thus,
      [`summary.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/teeMod_summary.md)):
      an HC2/CR2 correction for the residuals of the intent-to-treat
      effect estimates, and an HC1/CR1 correction or the residuals from
      the covariance adjustment model fit.
  2.  When units used for estimating intent-to-treat effects also
      contribute to fitting the covariance adjustment model, the boolean
      `loco_residuals` argument of
      [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
      indicates the residuals from effect estimation associated with
      individuals in the overlapping units should be replaced with
      residuals that use a leave-one-out estimate of the covariance
      adjustment model parameters.
- Degrees of freedom provided in
  [`summary.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/teeMod_summary.md)
  now reflect clustering: for CR0 and CR1 standard errors, the
  associated degrees of freedom are one less than the number of clusters
  used for estimation; for CR2 standard errors, degrees of freedom are
  computed using the approach of Imbens and Kolesár (2016) (see
  documentation of
  [`.compute_IK_dof()`](https://benbhansen-stats.github.io/propertee/dev/reference/dot-compute_IK_dof.md)
  for citation) and leverage the fast computational routine mentioned
  above
- Units of assignment in `StudySpecification` are now optional, though
  recommended.

### Bug Fixes

- Fix minor bugs with the `dichotomy` argument of
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)

## **propertee** 0.6.1

### Bug Fixes

- [`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md),
  [`unitid()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md),
  [`cluster()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md),
  [`uoa()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md),
  [`block()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md),
  and
  [`forcing()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  no longer fail automatically when passed a non-numeric or
  non-character variable. Now, they will first attempt to convert the
  variable to a character variable.
- An internal routine in
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  no longer fails when the environment in which
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  is called differs from the environment in which the
  `StudySpecification` associated with the `teeMod` is created.

## **propertee** 0.6.0

### New Features

- In their `coefficients` element, `teeMod` objects now report estimates
  of mean quantities in the control condition (response and, if
  applicable, predictions of response). See the
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  man page for further details.
- Introduces
  [`etc()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  (effect of the treatment on controls) and
  [`ato()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  (overlap-weighted effect) weighting functions. `atc` is an alias for
  [`etc()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md),
  while `olw`, `owt`, and `pwt` are aliases for
  [`ato()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md).
- Bump minimum R version to 4.1.0 to allow internal usage of pipes and
  anonymous functions.

### Bug Fixes

- When
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)-type
  functions are called with `data` arguments that do not include rows
  associated with all units of assignment specified in the
  `StudySpecification` object, the resulting weights reflect assignment
  probabilities across all units of assignment in the
  `StudySpecification`, not only those represented in `data`.

## **propertee** 0.5.2

### New Features

- Users fitting multiple `teeMod` objects to test multiple outcome
  variables or different levels of a factor treatment variable can pass
  those models to an `mmm` object from the `multcomp` package, then pass
  the `mmm` object to `multcomp`’s `glht()` function to obtain standard
  errors estimated using
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md).
  This requires passing `vcov=vcov_tee` to `glht()`. On a technical
  note, `**propertee**` now “suggests” installing `multcomp`.

### Bug Fixes

- Weights produced by
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  and
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  no longer produce 0’s or NA’s for `StudySpecification` objects created
  with multiple columns specifying the blocking scheme

## **propertee** 0.5.1

### Bug Fixes

- Variance estimation routine fixes miscalculations in the case when
  there is covariance adjustment and rows are omitted from the
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  fit due to `NA`’s in the covariates or treatment assignment

## **propertee** 0.5.0

### Major Changes

- All references to “design” have been changed to “specification”.
  - `*_design` is now `*_spec` (e.g. `rct_design` is now `rct_spec`)
  - `Design` objects are now `StudySpecification` objects
  - The `design=` argument to
    [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
    is now `specification=`.

## **propertee** 0.4.1

### New Features

- Passing `absorb = TRUE` to `lmitt` without specifying a block proceeds
  as if the entire sample is a single block.

### Bug Fixes

- Fixed bug where use of `dichotomy` and moderator variables in
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  could lead to errors due to too long of a formula.

## **propertee** 0.4.0

### New Features

- [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md),
  weights calculation functions
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  and
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md),
  and assignment vector generation function
  [`assigned()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)
  now accept a `dichotomy` argument that can be used for studies with
  time-varying treatment assignment. The `Design` object, unlike before,
  will not carry information about this dichotomization. Instead, the
  information stored there reflecting when units were assigned to
  treatment (if they were assigned to treatment) will be leveraged to
  create inverse probability of assignment weights and assignment
  indicators for datasets that have longitudinal data for the study
  units.

### Bug Fixes

- Standard error calculations no longer error when a `by` column is used
  to uniquely identify rows in the covariance adjustment or effect
  estimation sample that cannot be distinguished with information in the
  `Design` alone

## **propertee** 0.3.10

### Bug Fixes

- Linking unit of assignments to clusters for variance estimation no
  longer errors when `Design` objects are created with a `tibble`
- [`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
  does not error with covariance adjustment models fit with
  [`robustbase::glmrob()`](https://rdrr.io/pkg/robustbase/man/glmrob.html)

## **propertee** 0.3.9

### Bug Fixes

- Scaling constants have been updated in
  [`estfun.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/estfun.teeMod.md)
  to account for a previously missing factor of sqrt(n / n_C) applied to
  contributions to the covariance adjustment model estimating equations

## **propertee** 0.3.8

### Breaking Changes

- When model-based standard errors clustered at the level of assignment
  are called for in a blocked design,
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  clusters units of assignment in small blocks, blocks with only one
  treated or control unit, together.

## **propertee** 0.3.7

### Breaking Changes

- [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  scales estimating equations using different constants than it did
  before

## **propertee** 0.3.6

### Bug Fixes

- Previous procedure for aligning contributions to estimating equations
  from first-stage and second-stage models failed when column(s) used
  for alignment had NA’s. Outputs of
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  were liable to change from call to call as a result. This has been
  fixed.

## **propertee** 0.3.5

### Improvements

- Diagonal elements of
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  matrices lacking sufficient degrees of freedom for estimation are
  returned as NA’s rather than numeric zeros. This is a deviation from
  the `sandwich` package that aims to provide clarity to results that
  may otherwise appear as negative diagonal elements of the vcov matrix

### Bug Fixes

- When
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  is called with a blocked design and `absorb=TRUE`, the block-centered
  assignment and, if applicable, moderator and assignment:moderator
  interaction columns, are no longer centered on the grand mean of the
  column. This ensures blocks that do not satisfy positivity of the
  assignment variable (or positivity within a factor level) do not
  contribute to effect estimation
- [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  now accepts references to formula objects

## **propertee** 0.3.4

### Improvements

- Computational performance for `estfun.teeMod` has been improved

### Bug Fixes

- No more errors due to under-the-hood duplication of a moderator
  variable
- `absorb=TRUE` estimates have been corrected in the case when all
  observations in a stratum have 0 weights due to only treated or
  control units of assignment existing in the stratum

## **propertee** 0.3.3

### Added Features

- [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  can accept user-created variance estimation functions that start with
  the prefix `.vcov_`; the `type` argument should take the rest of the
  function name as an input
- Variance estimation for robust GLM’s (models fit using
  [`robustbase::glmrob`](https://rdrr.io/pkg/robustbase/man/glmrob.html))
  is now accommodated
- HC1 variance estimates are now accommodated

## **propertee** 0.3.2

### Added Features

- Effect estimation for continuous moderator variables is now supported

### Non-Breaking Changes

- [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
  will return NA’s for the entries of the covariance matrix that lack
  sufficient degrees of freedom for an estimate. Informative warnings
  will accompany the matrix, further indicating which standard errors
  have been NA’d out.

### Bug Fixes

- Functions for generating weights,
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  and
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md),
  return weights of 0 rather than infinity for blocks that contain
  treated units but no control units.
- Prior covariate adjustment fits were previously incorporated into
  variance estimation differently depending on whether one created a
  `SandwichLayer` object before calling
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  or called
  [`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
  in the `offset` argument of the
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  call. This has been corrected, and both ways return the same variance
  estimates.
- Covariate adjustment models that admit rectangular bread matrices,
  such as those produced by
  [`robustbase::lmrob`](https://rdrr.io/pkg/robustbase/man/lmrob.html),
  are now accommodated given the reformulated estimating equations in
  versions `v0.1.1` and later.
- A contrasts error raised by
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) in
  certain
  [`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
  calls has been resolved.

## **propertee** 0.3.1

### Breaking Changes

- We now order `teeMod` objects’ matrix of estimating equations based on
  user-specified ID columns or unit of assignment ID’s.
- The [`stats::update`](https://rdrr.io/r/stats/update.html) function
  can no longer be called on `teeMod` objects.

### Non-Breaking Changes

- `teeMod` objects now have `lmitt_call` slots.
- `summary` calls on `teeMod` objects accept `vcov.type` arguments to
  specify the desired standard error calculation shown in the output.
  Acceptable types follow the documentation for `vcov_tee`.
- Shown or printed `teeMod` objects return more comprehensible labels
  for ITT effect outputs.

### R Version Compatibility

- Now compatible with R 4.3. Particularly, we advise users working with
  R 4.3 to avoid `expand.model.frame` calls on `teeMod` objects and
  instead use the internal function `.expand.model.frame_teeMod` when
  necessary.

## **propertee** 0.2.1

### Breaking Changes

- Stratum fixed effects and subgroup moderating effects can now be
  accounted for via the `absorb` argument. Previous versions did not
  properly support this functionality. Valid standard errors under
  absorption, however, have not been confirmed.

## **propertee** 0.1.1

### Breaking Changes

- We have reformulated the estimating equations used to derive standard
  errors. In estimation settings we accommodate, testing has not
  revealed any differences in standard error estimates between the
  previous and current estimating equations, but we do not assure this
  is the case for all possible situations.

## **propertee** 0.0.1

- Compatible with R 4.2.3
- Introduces functionality for direct adjusted and design-informed
  standard errors accommodating covariance adjustment in the model-based
  setting
- Cluster-robust standard errors can only be estimated using the HC0
  estimator
