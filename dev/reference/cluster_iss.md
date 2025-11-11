# Use properties of idempotent matrices to cheaply compute inverse symmetric square roots of cluster-specific subsets of projection matrices

Use properties of idempotent matrices to cheaply compute inverse
symmetric square roots of cluster-specific subsets of projection
matrices

## Usage

``` r
cluster_iss(
  tm,
  cluster_unit,
  cluster_ids = NULL,
  cluster_var = NULL,
  exclude = na.action(tm),
  tol = 1e-09,
  ...
)
```

## Arguments

- exclude:

  index of units to exclude from computing the correction; for example,
  if they're NA's
