# Confidence intervals with standard errors provided by `vcov.teeMod()`

An S3method for [`stats::confint`](https://rdrr.io/r/stats/confint.html)
that uses standard errors computed using
[`vcov.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/vcov.teeMod.md).
Additional arguments passed to this function, such as `cluster` and
`type`, specify the arguments of the
[`vcov.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/vcov.teeMod.md)
call.

## Usage

``` r
# S3 method for class 'teeMod'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  a fitted `teeMod` model

- parm:

  a specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  the confidence level required.

- ...:

  additional arguments to pass to
  [`vcov.teeMod()`](https://benbhansen-stats.github.io/propertee/dev/reference/vcov.teeMod.md)

## Value

A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter. These will be labelled as (1-level)/2 and 1 -
(1-level)/2 in % (by default 2.5% and 97.5%)

## Details

Rather than call
[`stats::confint.lm()`](https://rdrr.io/r/stats/confint.html),
`confint.teeMod()` calls
[`.confint_lm()`](https://benbhansen-stats.github.io/propertee/dev/reference/dot-confint_lm.md),
a function internal to the `propertee` package that ensures additional
arguments in the `...` of the `confint.teeMod()` call are passed to the
internal [`vcov()`](https://rdrr.io/r/stats/vcov.html) call.
