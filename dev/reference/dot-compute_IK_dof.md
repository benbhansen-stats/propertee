# Compute the degrees of freedom of a sandwich standard error with HC2 correction

Compute the degrees of freedom of a sandwich standard error with HC2
correction

## Usage

``` r
.compute_IK_dof(
  tm,
  ell,
  cluster = NULL,
  bin_y = FALSE,
  exclude = na.action(tm),
  tol = 1e-09
)
```

## References

Guido W. Imbens and Michael Kolesár. "Robust Standard Errors in Small
Samples: Some Practical Advice". In: *The Review of Economics and
Statistics* 98.4 (Oct. 2016), pp. 701-712.
