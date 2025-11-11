# Produce confidence intervals for linear models

Produce confidence intervals for linear models

## Usage

``` r
.confint_lm(object, parm, level = 0.95, ...)
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

`.confint_lm()` is a copy of
[`stats::confint.lm`](https://rdrr.io/r/stats/confint.html) but passes
arguments in `...` to the [`vcov()`](https://rdrr.io/r/stats/vcov.html)
call. When called on a `teeMod` model, this produces confidence
intervals where standard errors are computed based on the desired
formulation of the
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md)
call.
