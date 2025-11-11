# (Internal) Extracts treatment as binary `vector`

(Internal) Extracts treatment as binary `vector`

## Usage

``` r
.bin_txt(spec, data = NULL, dichotomy = NULL)
```

## Arguments

- spec:

  A `StudySpecification`, used to get treatment assignment information

- data:

  A dataframe with unit of assignment information and, if a dichotomy is
  provided, columns specified therein

- dichotomy:

  Optional, a formula. See the Details section of the
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  or
  [`att()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  help pages for information on specifying the formula

## Value

A `vector` of binary treatments

## Details

If a `dichotomy` is specified or the `StudySpecification` has a
treatment variable consisting only of 0/1 or `NA`, then returns the
binary treatment. Otherwise (it has a non-binary treatment and lacks a
dichotomy) it errors.
