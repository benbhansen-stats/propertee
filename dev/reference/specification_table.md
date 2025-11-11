# Table of elements from a `StudySpecification`

Produces a table (1-dimensional, or 2-dimensional if `y` is specified)
of the elements of the `StudySpecification`.

## Usage

``` r
specification_table(
  specification,
  x,
  y = NULL,
  sort = FALSE,
  decreasing = TRUE,
  use_var_names = FALSE,
  ...
)

stable(
  specification,
  x,
  y = NULL,
  sort = FALSE,
  decreasing = TRUE,
  use_var_names = FALSE,
  ...
)
```

## Arguments

- specification:

  A `StudySpecification` object

- x:

  One of "treatment", "unit of assignment", (synonym "uoa"), "block".
  Abbreviations are accepted. "unit of assignment" can be replaced by
  "unitid" or "cluster" if the `StudySpecification` was created with
  that element.

- y:

  Optionally, another string similar to `x`. A 1-dimensional table is
  produced if `y` is left at its default, `NULL`.

- sort:

  Ignored if `y` is not `NULL`. If `FALSE` (default), one-way table is
  sorted according to "names" of levels. If set to `TRUE`, one-way table
  is sorted according to values.

- decreasing:

  If `sort` is `TRUE`, choose whether to sort descending (`TRUE`,
  default) or ascending (`FALSE`).

- use_var_names:

  If `TRUE`, name dimensions of table returned by variable names. If
  `FALSE` (default), name by their function (e.g. "treatment" or
  "blocks"). Passing the `dnn` argument in `...` (an argument of
  [`table()`](https://rdrr.io/r/base/table.html)) overrides whatever is
  requested here.

- ...:

  additional arguments [`table()`](https://rdrr.io/r/base/table.html)

## Value

A table of the requested variables.

## Examples

``` r
data(simdata)
spec <- obs_spec(z ~ unit_of_assignment(uoa1, uoa2) + block(bid),
                  data = simdata)
specification_table(spec, "treatment")
#> treatment
#> 0 1 
#> 6 4 
specification_table(spec, "treatment", "block", sort = TRUE, use_var_names = TRUE)
#>    z
#> bid 0 1
#>   1 3 1
#>   2 2 1
#>   3 1 2
```
