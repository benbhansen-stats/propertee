# Generate matrix of estimating equations for `lmrob()` fit

Generate matrix of estimating equations for `lmrob()` fit

Extract bread matrix from an `lmrob()` fit

## Usage

``` r
# S3 method for class 'lmrob'
estfun(x, ...)

# S3 method for class 'lmrob'
bread(x, ...)
```

## Arguments

- x:

  An `lmrob` object produced using an MM/SM estimator chain

- ...:

  Additional arguments to be passed to `bread`

## Value

A \\n\times \\(p+1) matrix where the first column corresponds to the
scale estimate and the remaining \\p\\ colums correspond to the
coefficients

A \\p\times \\(p+1) matrix where the first column corresponds to the
scale estimate and the remaining \\p\\ colums correspond to the
coefficients

## Details

This is part of a workaround for an issue in the robustbase code
affecting sandwich covariance estimation. The issue in question is issue
\#6471, robustbase project on R-Forge. This function contributes to
providing sandwich estimates of covariance-adjusted standard errors for
robust linear covariance adjustment models.

This is part of a workaround for an issue in the robustbase code
affecting sandwich covariance estimation. The issue in question is issue
\#6471, robustbase project on R-Forge. This function contributes to
providing sandwich estimates of covariance-adjusted standard errors for
robust linear covariance adjustment models.

## Author

Ben B. Hansen
