# (Internal) Compute \\(I\_{i} - H\_{ii})^{-1/2}\\ as part of CR2 variance estimates

(Internal) Compute \\(I\_{i} - H\_{ii})^{-1/2}\\ as part of CR2 variance
estimates

## Usage

``` r
cluster_iss(tm, cluster_unit, cluster_ids = NULL, ...)
```

## Arguments

- tm:

  `teeMod` object.

- cluster_unit:

  cluster to subset observations to. Must be found in `cluster_ids`.

- cluster_ids:

  optional, ID's for clustering standard errors. If not provided,
  default is the ID's associated with `cluster` if provided, otherwise
  the unit of assignment ID's.

- ...:

  Additional arguments passed from calls higher in the stack. One that
  may be used is `cluster`, which identifies the clustering variable if
  `cluster_ids` is not provided. If `cluster` is not provided, the
  default is the unit of assignment variable specified in the
  `StudySpecification`.

## Details

The notation \\I\_{i}\\ and \\H\_{ii}\\ comes from Bell and McCaffrey
(2002). The matrix \\I\_{i}\\ is an identity matrix with number of rows
equal to the number of observations in cluster \\i\\, and \\H\_{ii}\\
subsets the hat matrix associated with the regression fit stored in `tm`
to the rows associated with observations in cluster \\i\\.

When possible, the function uses the method in Wasserman (2026) to
cheaply compute the inverse symmetric square root.

`cluster_ids` should be ordered in alignment with the dataframe passed
to
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md).
It should not exclude NA's because the function will exclude them.

## References

Robert M. Bell and Daniel F. McCaffrey. Bias Reduction in Standard
Errors for Linear Regression with Multi-Stage Samples. *Survey
Methodology*, 28(2):169-181, December 2002.

Joshua Wasserman. Methods for Causal Inference in Settings with
Clustered Data Subject to Missingness and Measurement Error. Unpublished
thesis, June 2026.
