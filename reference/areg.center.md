# Group-center akin to Stata's `areg`

From the Stata documentation: `areg` begins by recalculating Y and X and
to have mean 0 within the groups specified by `absorb()`. The overall
mean of each variable is then added back in.

## Usage

``` r
areg.center(mm, grp, wts = NULL, grand_mean_center = FALSE)
```

## Arguments

- mm:

  Matrix of variables to center

- grp:

  Group to center on

- wts:

  Optional weights

- grand_mean_center:

  Optional center output at `mean(var)`

## Value

Vector of group-centered values
