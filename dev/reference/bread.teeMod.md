# Extract bread matrix from a `teeMod` model fit

An S3method for
[`sandwich::bread`](https://sandwich.R-Forge.R-project.org/reference/bread.html)
that extracts the bread of the direct adjustment model sandwich
covariance matrix.

## Usage

``` r
# S3 method for class 'teeMod'
bread(x, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

- ...:

  arguments passed to methods

## Value

A variance-covariance matrix with row and column entries for the
estimated coefficients in `x`, the marginal mean outcome in the control
condition, the marginal mean `offset` in the control condition (if an
`offset` is provided), and if a moderator variable is specified in the
formula for `x`, the mean interaction in the control condition of the
outcome and `offset` with the moderator variable

## Details

This function is a thin wrapper around
[`.get_tilde_a22_inverse()`](https://benbhansen-stats.github.io/propertee/dev/reference/sandwich_elements_calc.md).
