# (Internal) Locate data in call stack

Whenever a function in a model
([`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)/[`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)/[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md)/[`assigned()`](https://benbhansen-stats.github.io/propertee/reference/AssignedAliases.md))
is called without an explicit `data=` argument, this will attempt to
extract the data from the model itself.

## Usage

``` r
.get_data_from_model(which_fn, form = NULL, by = NULL)
```

## Arguments

- which_fn:

  Identify calling function, "weights" or "assigned", helps separate
  logic for the two functions.

- form:

  Formula on which to apply
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html). See
  details

- by:

  translation of unit of assignment/unitid/cluster ID names, passed down
  from weights.

## Value

`data.frame`

## Details

The `form` specifies what columns of the data are needed. For current
use cases
([`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)/[`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
and
[`assigned()`](https://benbhansen-stats.github.io/propertee/reference/AssignedAliases.md)),
this will be only the unit of assignment variables, so e.g.
`form = ~ uoavar`, to enable merging of UOA level variables to the model
data. However, this can easily be expanded if other variables are
needed.
