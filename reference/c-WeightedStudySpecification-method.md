# Concatenate weights

Given several variations of weights generated from a single
`StudySpecification`, combine into a single weight.

## Usage

``` r
# S4 method for class 'WeightedStudySpecification'
c(x, ..., warn_dichotomy_not_equal = FALSE)
```

## Arguments

- x, :

  .. a `WeightedStudySpecification` object, typically created from
  [`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
  or
  [`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)

- ...:

  any number of additional `WeightedStudySpecification` objects with
  equivalent `StudySpecification` to `x` and eachother

- warn_dichotomy_not_equal:

  if `FALSE` (default), `WeightedStudySpecification`s are considered
  equivalent even if their `dichotomy` differs. If `TRUE`, a warning is
  produced.

## Value

A numeric `vector` with the weights concatenated in the input order.

## Details

Concatenating `WeightedStudySpecification` objects with
[`c()`](https://rdrr.io/r/base/c.html) requires both individual
`WeightedStudySpecification` objects to come from the same
`StudySpecification` and have the same target (e.g all created with
[`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
or all created with
[`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md),
no mixing-and-matching). All arguments to
[`c()`](https://rdrr.io/r/base/c.html) must be
`WeightedStudySpecification`.

`WeightedStudySpecification` objects may be concatenated together even
without having the same `@dichotomy` slot. This procedure only prompts a
warning for differing dichotomies if the argument
`warn_dichotomy_not_equal` is set to `TRUE`.

## Examples

``` r
data(simdata)
spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
w1 <- ate(spec, data = simdata[1:30,])
w2 <- ate(spec, data = simdata[31:40,])
w3 <- ate(spec, data = simdata[41:50,])
c_w <- c(w1, w2, w3)
c(length(w1), length(w2), length(w3), length(c_w))
#> [1] 30 10 10 50

spec <- rct_spec(dose ~ unit_of_assignment(uoa1, uoa2), data = simdata)
w1 <- ate(spec, data = simdata[1:10, ], dichotomy = dose >= 300 ~ .)
w2 <- ate(spec, data = simdata[11:30, ], dichotomy = dose >= 200 ~ .)
w3 <- ate(spec, data = simdata[31:50, ], dichotomy = dose >= 100 ~ .)
c_w <- c(w1, w2, w3)
```
