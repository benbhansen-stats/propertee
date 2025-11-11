# (Internal) Worker function for weight calculation

Called from
[`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
or
[`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md).

## Usage

``` r
.weights_calc(specification, target, weightAlias, dichotomy, by, data)
```

## Arguments

- specification:

  a `StudySpecification` object created by one of
  [`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
  [`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
  or
  [`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md).

- target:

  One of "ate" or "ett";
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  and
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  chooses these automatically.

- weightAlias:

  An alias for the weight target, currently one of "ate", "ett", "att".
  Chosen by
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  and
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  automatically.

- dichotomy:

  optional; a formula defining the dichotomy of the treatment variable
  if it isn't already `0`/`1`. See details of help for
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  or
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  e.g. for details.

- by:

  optional; named vector or list connecting names of unit of assignment/
  variables in `specification` to unit of assignment/cluster variables
  in `data`. Names represent variables in the StudySpecification; values
  represent variables in the data. Only needed if variable names differ.

- data:

  optionally the data for the analysis to be performed on. May be
  excluded if these functions are included as the `weights` argument of
  a model.

## Value

a `WeightedStudySpecification` object
