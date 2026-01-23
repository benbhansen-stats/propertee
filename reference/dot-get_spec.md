# (Internal) Locate a `StudySpecification` in the call stack

[`assigned()`](https://benbhansen-stats.github.io/propertee/reference/AssignedAliases.md)/[`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)/[`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)/[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md)
all need the `StudySpecification` to operate. If any are called in the
model without a `specification=` argument, this function sees if it can
find the `StudySpecification` in another of these functions.

## Usage

``` r
.get_spec(NULL_on_error = FALSE)
```

## Arguments

- NULL_on_error:

  if `TRUE`, returns `NULL` if a `StudySpecification` object is not
  found rather than throwing an error.

## Value

A `StudySpecification`, or `NULL` if `NULL_on_error` is `TRUE` and the
`StudySpecification` can't be found.

## Details

Note that it will never look inside
[`assigned()`](https://benbhansen-stats.github.io/propertee/reference/AssignedAliases.md)
(gets complicated in formulas), only in weights or
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md).
E.g.

`lm(y ~ assigned(), weights = ate(spec), offest = cov_adj(mod1))`

`lm(y ~ assigned(), weights = ate(), offest = cov_adj(mod1, specification = spec))`

will both work, but

`lm(y ~ assigned(spec), weights = ate(), offest = cov_adj(mod1))`

will fail.
