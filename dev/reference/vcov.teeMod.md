# Compute variance-covariance matrix for fitted `teeMod` model

An S3method for [`stats::vcov`](https://rdrr.io/r/stats/vcov.html) that
computes standard errors for `teeMod` models using
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md).

## Usage

``` r
# S3 method for class 'teeMod'
vcov(object, ...)
```

## Arguments

- object:

  a fitted `teeMod` model

- ...:

  additional arguments to
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md).

## Value

A variance-covariance matrix with row and column entries for the
estimated coefficients in `x`, the marginal mean outcome in the control
condition, the marginal mean `offset` in the control condition (if an
`offset` is provided), and if a moderator variable is specified in the
formula for `x`, the mean interaction in the control condition of the
outcome and `offset` with the moderator variable

## Details

`vcov.teeMod()` wraps around
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md),
so additional arguments passed to `...` will be passed to the
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
call. See documentation for
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
for information about necessary arguments.
