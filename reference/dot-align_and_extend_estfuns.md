# (Internal) Align the dimensions and rows of direct adjustment and covariance adjustment model estimating equations matrices

(Internal) Align the dimensions and rows of direct adjustment and
covariance adjustment model estimating equations matrices

## Usage

``` r
.align_and_extend_estfuns(x, ctrl_means_ef_mat = NULL, by = NULL, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

- ctrl_means_ef_mat:

  optional, a matrix of estimating equations corresponding to the
  estimates of the marginal (and possibly conditional) means of the
  outcome and `offset` in the control condition. These are aligned and
  extended in the same way as the matrix of estimating equations for `x`
  and `cbind`ed to them

- by:

  optional, a character vector indicating columns that uniquely identify
  rows in the dataframe used for fitting `x` and the dataframe passed to
  the `data` argument of the covariance adjustment model fit. The
  default is `NULL`, in which case the unit of assignment columns
  specified in the `StudySpecification` slot of `x` are used.

- ...:

  mostly arguments passed to methods, but the special case is the
  argument `loco_residuals`, which indicates the offsets in the
  residuals of `x` should be replaced by versions that use
  leave-one-cluster-out estimates of the covariance model

## Value

A list of two matrices, one being the aligned contributions to the
estimating equations for the direct adjustment model, and the other
being the aligned contributions to the covariance adjustment model.

## Details

`.align_and_extend_estfuns()` first extracts the matrices of
contributions to the empirical estimating equations for the direct
adjustment and covariance adjustment models; then, it pads the matrices
with zeros to account for units of observation that appear in one
model-fitting sample but not the other; finally it orders the matrices
so units of observation (or if unit of observation-level ordering is
impossible, units of assignment) are aligned.
