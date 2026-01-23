# Make ID's to pass to the `cluster` argument of `vcov_tee()`

`.make_uoa_ids()` returns a factor vector of cluster ID's that align
with the order of the units of observations' contributions in
[`estfun.teeMod()`](https://benbhansen-stats.github.io/propertee/reference/estfun.teeMod.md).
This is to ensure that when
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/reference/var_estimators.md)
calls
[`sandwich::meatCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html),
the `cluster` argument aggregates the correct contributions to
estimating equations within clusters.

## Usage

``` r
.make_uoa_ids(x, vcov_type, cluster = NULL, ...)
```

## Arguments

- x:

  a fitted `teeMod` object

- vcov_type:

  a string indicating model-based or design-based covariance estimation.
  Currently, "MB", "CR", and "HC" are the only strings registered as
  indicating model-based estimation.

- cluster:

  character vector or list; optional. Specifies column names that appear
  in both the covariance adjustment and direct adjustment model
  dataframes. Defaults to NULL, in which case unit of assignment columns
  indicated in the StudySpecification will be used for clustering. If
  there are multiple clustering columns, they are concatenated together
  for each row and separated by "\_".

- ...:

  arguments passed to methods

## Value

A vector with length equal to the number of unique units of observation
used to fit the two models. See Details of
[`estfun.teeMod()`](https://benbhansen-stats.github.io/propertee/reference/estfun.teeMod.md)
for the method for determining uniqueness.
