# Compute residuals for a `teeMod` object with leave-one-out estimates of the `offset`

Compute residuals for a `teeMod` object with leave-one-out estimates of
the `offset`

## Usage

``` r
.compute_loo_resids(x, cluster, ...)
```

## Arguments

- x:

  a `teeMod` object

- cluster:

  vector of column names that identify clusters

## Details

The residual for any observation also used for fitting the
`fitted_covariance_model` stored in the `offset` of `x` is replaced by
an estimated residual that uses a cluster leave-one-out estimate of the
`fitted_covariance_model` for generating a value of the `offset`.
