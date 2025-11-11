# (Internal) Applies dichotomy to treatment

Given a dichotomy formula and a `data.frame` with a treatment variable
and any variables in the formula, returns a `vector` containing only
`0`, `1`, or `NA`.

## Usage

``` r
.apply_dichotomy(txt, dichotomy)
```

## Arguments

- txt:

  A named `data.frame` containing a column of the treatment, such as
  that produed by `treatment(myspecification)`, and any variables
  specified in `dichotomy`.

- dichotomy:

  A formula specifying how to dichotomize the non-binary treatment
  column in `txt` (or a call that evaluates to a formula). See the
  Details section of the
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  or
  [`att()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  help pages for information on specifying this formula

## Value

A `vector` of binary treatments
