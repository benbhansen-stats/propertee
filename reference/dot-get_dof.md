# (Internal) Compute the degrees of freedom of a contrast of a sandwich variance estimate associated with a `teeMod`

(Internal) Compute the degrees of freedom of a contrast of a sandwich
variance estimate associated with a `teeMod`

## Usage

``` r
.get_dof(
  x,
  ell,
  vcov.type,
  dof.type = c("stata", "IK"),
  cluster_ids = NULL,
  cluster = NULL,
  ...
)
```

## Arguments

- x:

  `teeMod` object.

- ell:

  numeric.

- vcov.type:

  character.

- dof.type:

  character, either `IK` or `stata`.

- cluster_ids:

  optional, vector of ID's for clustering degrees of freedom estimate.
  If not provided, default is the ID's associated with `cluster`.

- cluster:

  optional, character identifiying the clustering variable if
  `cluster_ids` is not provided. If not provided, defaults to the unit
  of assignment columns specified in the `StudySpecification`.

- ...:

  Additional arguments passed from calls higher up the stack. These
  arguments are not used within this function.

## Details

With `dof_type="IK"`, one obtains a degrees of freedom estimate as given
by the procedure propose in Imbens and Kolesár (2016). For more details
see the documentation of
[`.compute_IK_dof()`](https://benbhansen-stats.github.io/propertee/reference/dot-compute_IK_dof.md).
With `dof_type="stata"`, one obtains degrees of freedom equal to the
number of clusters less one, the default provided by STATA software (see
Cameron and Miller, 2015 or Bell and McCaffrey, 2002).

`ell` should be an integer specifying a location in the vector of
estimated coefficients (ignoring the coefficients suffixed by
`:(Intercept)`) or a vector of the same length as the vector of
estimated coefficients (ignoring the coefficients suffixed by
`:(Intercept)`).

`vcov.type` takes the same arguments as the `type` argument in
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/reference/var_estimators.md).

`cluster_ids` should be ordered in alignment with the dataframe passed
to
[`lmitt()`](https://benbhansen-stats.github.io/propertee/reference/lmitt.md).
It should not exclude NA's because the function will exclude them where
necessary.

## References

Guido W. Imbens and Michael Kolesár. Robust Standard Errors in Small
Samples: Some Practical Advice". *The Review of Economics and
Statistics*, 98(4):701-712, October 2016.

A. Colin Cameron and Douglas L. Miller. A Practitioner's Guide to
Cluster-Robust Inference. *The Journal of Human Resources*,
50(2):317-372, 2015.

Robert M. Bell and Daniel F. McCaffrey. Bias Reduction in Standard
Errors for Linear Regression with Multi-Stage Samples. *Survey
Methodology*, 28(2):169-181, December 2002.
