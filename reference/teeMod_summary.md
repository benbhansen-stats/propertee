# Summarizing `teeMod` objects

[`summary()`](https://rdrr.io/r/base/summary.html) method for class
`teeMod`

## Usage

``` r
# S3 method for class 'teeMod'
summary(object, vcov.type = "HC0", ...)

# S3 method for class 'summary.teeMod'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- object:

  `teeMod` object

- vcov.type:

  A string indicating the desired variance estimator. See
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/reference/var_estimators.md)
  for details on accepted types.

- ...:

  Additional arguments to
  [`vcov_tee()`](https://benbhansen-stats.github.io/propertee/reference/var_estimators.md),
  such as the desired finite sample heteroskedasticity-robust standard
  error adjustment.

- x:

  `summary.teeMod` object

- digits:

  the number of significant digits to use when printing.

- signif.stars:

  logical. If ‘TRUE’, ‘significance stars’ are printed for each
  coefficient.

## Value

object of class `summary.teeMod`

## Details

If a `teeMod` object is fit with a `SandwichLayer` offset, then the
usual [`stats::summary.lm()`](https://rdrr.io/r/stats/summary.lm.html)
output is enhanced by the use of covariance-adjusted sandwich standard
errors, with t-test values recalculated to reflect the new standard
errors.
