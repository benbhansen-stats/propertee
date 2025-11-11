# Variance/Covariance for `teeMod` objects

Compute robust sandwich variance estimates with optional covariance
adjustment

## Usage

``` r
vcov_tee(x, type = NULL, cluster = NULL, ...)

.vcov_DB0(x, ...)

.vcov_DB(x, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

- type:

  a string indicating the desired bias correction for the residuals of
  `x`. Default makes no bias correction. See Details for supported types

- cluster:

  a vector indicating the columns that define clusters. The default is
  the unit of assignment columns in the `StudySpecification` stored in
  `x`. These columns should appear in the dataframe used for fitting `x`
  as well as the dataframe passed to the covariance model fit in the
  case of prior covariance adjustment. See Details

- ...:

  arguments to be passed to the internal variance estimation function,
  such as `cov_adj_rcorrect` and `loco_residuals`. If `x` has a
  `SandwichLayer` object in its offset, The former specifies the bias
  correction to the residuals of the covariance model, and the latter
  indicates whether the offset should be replaced with predictions from
  leave-one-cluster-out fits of the covariance adjustment model. See
  Details

## Value

A variance-covariance matrix with row and column entries for the
estimated coefficients in `x`, the marginal mean outcome in the control
condition, the marginal mean `offset` in the control condition (if an
`offset` is provided), and if a moderator variable is specified in the
formula for `x`, the mean interaction in the control condition of the
outcome and `offset` with the moderator variable

## Details

Variance estimates will be clustered on the basis of the columns
provided to `cluster` (or obtained by the default behavior). As a
result, providing `"HCx"` or `"CRx"` to `type` will produce the same
variance estimate given that `cluster` remains the same.

With prior covariance adjustment, unless the `data` argument of the
covariance model fit is the same as the `data` argument for fitting `x`
and the `StudySpecification` of `x` has been created with a formula of
the form `trt_col ~ 1`, the column(s) provided to `cluster` must appear
in the dataframes in both `data` arguments, even if the clustering
structure does not exist, per se, in the covariance adjustment sample.
For instance, in a finely stratified randomized trial, one might desire
standard errors clustered at the block level, but the covariance
adjustment model may include auxiliary units that did not participate in
the trial. In this case, in the `data` argument of the fitted covariance
model, the column(s) passed to `cluster` should have the block ID's for
rows overlapping with the `data` argument used for fitting `x`, and NA's
for any auxiliary units. `vcov_tee()` will treat each row with an NA as
its own cluster.

For ITT effect estimates without covariance adjustment, `type`
corresponds to the variance estimate desired. Supported options include:

- `"MB0"`, `"HC0"`, and `"CR0"` for model-based HC/CR0 standard errors

- `"MB1"`, `"HC1"`, and `"CR1"` for model-based HC/CR1 standard errors
  (for `"MB1"` and `"HC1"`, this is \\n/(n - 2)\\, and for `"CR1"`, this
  is \\g\cdot(n-1)/((g-1)\cdot(n-2))\\, where \\g\\ is the number of
  clusters in the sample used for fitting `x`)

- `"MB2"`, `"HC2"`, and `"CR2"` for model-based HC/CR2 standard errors

- `"DB0"` for design-based HC0 variance estimates

The `type` argument does not correspond to existing variance estimators
in the literature in the case of prior covariance adjustment. It
specifies the bias correction to the residuals of `x`, but the residuals
of the covariance model are corrected separately based on the
`cov_adj_rcorrect` argument. The `cov_adj_rcorrect` argument takes the
same options as `type` except `"DB0"`. When the covariance model
includes rows in the treatment condition for fitting, the residuals of
`x` are further corrected by having the values of `offset` replaced by
predictions that use coefficient estimates that leave out rows in the
same cluster (as defined by the `cluster` argument).

Design-based variance estimation is implemented for `teeMod` models
satisfying requirements including:

- The model uses `rct_spec` as `StudySpecification`

- The model only estimates a main treatment effect
