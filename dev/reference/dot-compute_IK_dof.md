# (Internal) Compute the degrees of freedom of a contrast of a sandwich variance estimate as proposed in Imbens and Kolesár (2016)

(Internal) Compute the degrees of freedom of a contrast of a sandwich
variance estimate as proposed in Imbens and Kolesár (2016)

## Usage

``` r
.compute_IK_dof(
  tm,
  ell,
  vcov.type,
  cluster_ids = NULL,
  cluster = NULL,
  tol = 1e-09
)
```

## Arguments

- tm:

  `teeMod` object.

- ell:

  numeric vector.

- vcov.type:

  character.

- cluster_ids:

  optional, vector of ID's for clustering degrees of freedom estimate.
  If not provided, default is the ID's associated with `cluster`.

- cluster:

  optional, character identifiying the clustering variable if
  `cluster_ids` is not provided. If not provided, defaults to the unit
  of assignment columns specified in the `StudySpecification`.

- tol:

  optional, numeric. Should not be changed.

## Details

`ell` should be a vector of length equal to the number of estimated
coefficients in the `teeMod` object. This excludes coefficients printed
in `show.teeMod` with `:(Intercept)` suffixes. The degrees of freedom
for a single standard error will specify for `ell` a vector of all zeros
except one element, which will have a 1 in the location corresponding to
the coefficient of interest.

`cluster_ids` should be ordered in alignment with the dataframe passed
to
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md).
It should not exclude NA's because the function will exclude them where
necessary.

`vcov.type` takes the same arguments as the `type` argument in
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/dev/reference/var_estimators.md).

## References

Guido W. Imbens and Michael Kolesár. Robust Standard Errors in Small
Samples: Some Practical Advice". *The Review of Economics and
Statistics*, 98(4):701-712, October 2016.
