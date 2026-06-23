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
  object.

- by:

  named vector; optional.

- specification:

  a `StudySpecification` object; optional. If not provided, it can be
  retrieved from the call stack if passed to suitable calls.

## Value

formula

## Details

This function is typically used prior to
[`.get_data_from_model()`](https://benbhansen-stats.github.io/propertee/reference/dot-get_data_from_model.md)
and incorporates information provided in a `by` argument to ensure the
necessary columns for merging the two samples are included in any
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) calls.
