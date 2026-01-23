# (Internal) Expand unit of assignment level weights to the level of the data

Helper function called during creation of the weights via
[`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
or
[`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)

## Usage

``` r
.join_spec_weights(
  weights,
  specification,
  target,
  weightAlias,
  data,
  dichotomy
)
```

## Arguments

- weights:

  a vector of weights sorted according to the `StudySpecification`

- specification:

  a `StudySpecification`

- target:

  One of "ate" or "ett"

- weightAlias:

  Any currently supported alias

- data:

  New data

- dichotomy:

  formula used to specify a dichotomy of a non-binary treatment
  variable. The output `WeightedStudySpecification` object will store
  this as its `dichotomy` slot, unless it is NULL, in which case it will
  be translated to an empty `formula`.

## Value

a `WeightedStudySpecification`
