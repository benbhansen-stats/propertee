# (Internal) Return ID's used to order observations in the covariance adjustment sample

(Internal) Return ID's used to order observations in the covariance
adjustment sample

## Usage

``` r
.sanitize_C_ids(x, id_col = NULL, sorted = FALSE, ...)
```

## Arguments

- x:

  a `SandwichLayer` object

- id_col:

  character vector or list; optional. Specifies column names that appear
  in botn the covariance adjustment and direct adjustmet samples.
  Defaults to NULL, in which case unit of assignment columns in the
  `SandwichLayer`'s `StudySpecification` slot will be used to generate
  ID's.

- ...:

  arguments passed to methods

## Value

A vector of length equal to the number of units of observation used to
fit the covariance adjustment model
