# (Internal) Order observations used to fit a `teeMod` model and a prior covariance adjustment model

(Internal) Order observations used to fit a `teeMod` model and a prior
covariance adjustment model

## Usage

``` r
.order_samples(x, by = NULL, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

- by:

  character vector of columns to get ID's for ordering from. Default is
  NULL, in which case unit of assignment ID's are used for ordering.

- ...:

  arguments passed to methods

## Value

A list of four named vectors. The `Q_not_C` element holds the ordering
for units of observation in the direct adjustment sample but not the
covariance adjustment samples; `Q_in_C` and `C_in_Q`, the ordering for
units in both; and `C_not_Q`, the ordering for units in the covariance
adjustment sample only. `Q_in_C` and `C_in_Q` differ in that the names
of the `Q_in_C` vector correspond to row indices of the original matrix
of estimating equations for the direct adjustment model, while the names
of `C_in_Q` correspond to row indices of the matrix of estimating
equations for the covariance adjustment model. Similarly, the names of
`Q_not_C` and `C_not_Q` correspond to row indices of the direct
adjustment and covariance adjustment samples, respectively. Ultimately,
the order of
[`.make_uoa_ids()`](https://benbhansen-stats.github.io/propertee/reference/dot-make_uoa_ids.md)
and
[`estfun.teeMod()`](https://benbhansen-stats.github.io/propertee/reference/estfun.teeMod.md)
is given by concatenating the vectors stored in `Q_not_C`, `Q_in_C`, and
`C_not_q`.

## Details

`.order_samples()` underpins the ordering for
[`.make_uoa_ids()`](https://benbhansen-stats.github.io/propertee/reference/dot-make_uoa_ids.md)
and
[`estfun.teeMod()`](https://benbhansen-stats.github.io/propertee/reference/estfun.teeMod.md).
This function orders the outputs of those functions, but also informs
how the original matrices of contributions to estimating equations need
to be indexed to align units of observations' contributions to both sets
of estimating equations.  
  
When a `by` argument is provided to
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md),
it is used to construct the order of `.order_samples()`.
