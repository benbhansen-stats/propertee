# Convert `lm` object into `teeMod`

Converts the output of [`lm()`](https://rdrr.io/r/stats/lm.html) into a
`teeMod` object, for standard errors that account for block and cluster
information carried with the `lm`'s weights, and/or an offset
incorporating predictions of the outcome from a separate model.

## Usage

``` r
as.lmitt(x, specification = NULL)

as.teeMod(x, specification = NULL)
```

## Arguments

- x:

  `lm` object with weights containing a `WeightedStudySpecification`, or
  an offset from
  [`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md).

- specification:

  Optional, explicitly specify the `StudySpecification` to be used. If
  the `StudySpecification` is specified elsewhere in `x` (e.g. passed as
  an argument to any of
  [`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md),
  [`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md),
  [`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md)
  or
  [`assigned()`](https://benbhansen-stats.github.io/propertee/reference/AssignedAliases.md))
  it will be found automatically and does not need to be passed here as
  well. If different `StudySpecification` objects are passed (either
  through the `lm` in weights or covariance adjustment, or through this
  argument), an error will be produced.

## Value

`teeMod` object

## Details

The formula with which `x` was created must include a treatment
identifier (e.g.
[`assigned()`](https://benbhansen-stats.github.io/propertee/reference/AssignedAliases.md)).
If a model-based offset is incorporated, the model's predictions would
have to have been extracted using
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md)
as opposed to `predict{}` in order for `teeMod` standard error
calculations to reflect propagation of error from these predictions.
This mechanism only supports treatment main effects: to estimate
interactions of treatment assignment with a moderator variable, use
[`lmitt()`](https://benbhansen-stats.github.io/propertee/reference/lmitt.md)
instead of [`lm()`](https://rdrr.io/r/stats/lm.html) followed by
`as.lmitt()`.
