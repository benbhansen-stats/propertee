# (Internal) Add columns for merging covariance adjustment and direct adjustment samples to model formula

(Internal) Add columns for merging covariance adjustment and direct
adjustment samples to model formula

## Usage

``` r
.update_ca_model_formula(model, by = NULL, specification = NULL)
```

## Arguments

- model:

  any model that inherits from a `glm`, `lm`, or
  [`robustbase::lmrob`](https://rdrr.io/pkg/robustbase/man/lmrob.html)
  object

- by:

  optional; a string or named vector of unique identifier columns in the
  data used to create `specification` and the data used to fit the
  covariance adjustment model. Default is NULL, in which case unit of
  assignment columns are used for identification (even if they do not
  uniquely identify units of observation). If a named vector is
  provided, names should represent variables in the data used to create
  `specification`, while values should represent variables in the
  covariance adjustment data.

- specification:

  a `StudySpecification` object. Default is NULL, in which case a
  `StudySpecification` object is sought from higher up the call stack.

## Value

formula

## Details

This function is typically used prior to
[`.get_data_from_model()`](https://benbhansen-stats.github.io/propertee/dev/reference/dot-get_data_from_model.md)
and incorporates information provided in a `by` argument to ensure the
necessary columns for merging the two samples are included in any
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) calls.
