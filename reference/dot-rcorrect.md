# (Internal) Bias correct residuals contributing to standard errors of a `teeMod`

(Internal) Bias correct residuals contributing to standard errors of a
`teeMod`

## Usage

``` r
.rcorrect(resids, x, model, type, ...)
```

## Arguments

- resids:

  numeric vector of residuals to correct

- x:

  teeMod object

- model:

  string indicating which model the residuals are from. `"itt"`
  indicates correction to the residuals of `x`, and `"cov_adj"`
  indicates correction to the residuals of the covariance adjustment
  model. This informs whether corrections should use information from
  `x` or the `fitted_covariance_model` slot of the `SandwichLayer`
  object in the `offset` for corrections.

- type:

  string indicating the desired bias correction. Can be one of
  `"(HC/CR/MB)0"`, `"(HC/CR/MB)1"`, or `"(HC/CR/MB)2"`

- ...:

  additional arguments passed from up the call stack; in particular, the
  `cluster_cols` argument, which informs whether to cluster and provide
  CR2 corrections instead of HC2 corrections, as well as the correction
  for the number of clusters in the CR1 correction. This may also
  include a `by` argument.
