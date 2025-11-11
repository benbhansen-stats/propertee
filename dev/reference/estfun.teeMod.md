# Extract empirical estimating equations from a `teeMod` model fit

An S3method for
[`sandwich::estfun`](https://sandwich.R-Forge.R-project.org/reference/estfun.html)
for producing a matrix of contributions to the direct adjustment
estimating equations.

## Usage

``` r
# S3 method for class 'teeMod'
estfun(x, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

- ...:

  arguments passed to methods, most importantly those that define the
  bias corrections for the residuals of `x` and, if applicable, a
  `fitted_covariance_model` stored in its offset

## Value

An \\n\times k\\ matrix of empirical estimating equations for `x`. `k`
includes the model intercept, main effects of treatment and moderator
variables, any moderator effects, and marginal and conditional means of
the outcome (and `offset`, if provided) in the control condition. See
Details for definition of \\n\\.

## Details

If a prior covariance adjustment model has been passed to the `offset`
argument of the `teeMod` model using
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md),
`estfun.teeMod()` incorporates contributions to the estimating equations
of the covariance adjustment model.  
  
The covariance adjustment sample may not fully overlap with the direct
adjustment sample, in which case `estfun.teeMod()` returns a matrix with
the same number of rows as the number of unique units of observation
used to fit the two models. Uniqueness is determined by matching units
of assignment used to fit the covariance adjustment model to units of
assignment in the `teeMod` model's `StudySpecification` slot; units of
observation within units of assignment that do not match are additional
units that add to the row count.  
  
The`by` argument in
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
can provide a column or a pair of columns (a named vector where the name
specifies a column in the direct adjustment sample and the value a
column in the covariance adjustment sample) that uniquely specifies
units of observation in each sample. This information can be used to
align each unit of observation's contributions to the two sets of
estimating equations. If no `by` argument is provided and units of
observation cannot be uniquely specified, contributions are aligned up
to the unit of assignment level. If standard errors are clustered no
finer than that, they will provide the same result as if each unit of
observation's contributions were aligned exactly.  
  
This method incorporates bias corrections made to the residuals of `x`
and, if applicable, the covariance model stored in its `offset`. When
its crossproduct is taken (perhaps after suitable summing across rows
within clusters), it provides a heteroskedasticity- (or cluster-) robust
estimate of the meat matrix of the variance-covariance of the parameter
estimates in `x`.
