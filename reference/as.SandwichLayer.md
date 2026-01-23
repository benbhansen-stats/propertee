# Convert a `PreSandwichLayer` to a `SandwichLayer` with a `StudySpecification` object

`as.SandwichLayer()` uses the `StudySpecification` object passed to the
`specification` argument to populate the slots in a `SandwichLayer`
object that a `PreSandwichLayer` does not have sufficient information
for.

## Usage

``` r
as.SandwichLayer(x, specification, by = NULL, Q_data = NULL)
```

## Arguments

- x:

  a `PreSandwichLayer` object

- specification:

  a `StudySpecification` object

- by:

  optional; a string or named vector of unique identifier columns in the
  data used to create `specification` and the data used to fit the
  covariance adjustment model. Default is NULL, in which case unit of
  assignment columns are used for identification (even if they do not
  uniquely identify units of observation). If a named vector is
  provided, names should represent variables in the data used to create
  `specification`, while values should represent variables in the
  covariance adjustment data.

- Q_data:

  dataframe of direct adjustment sample, which is needed to generate the
  `keys` slot of the `SandwichLayer` object. Defaults to NULL, in which
  case if `by` is NULL, the data used to create `specification` is used,
  and if `by` is not NULL, appropriate data further up the call stack
  (passed as arguments to
  [`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md)
  or
  [`lmitt.formula()`](https://benbhansen-stats.github.io/propertee/reference/lmitt.md),
  for example) is used.

## Value

a `SandwichLayer` object
