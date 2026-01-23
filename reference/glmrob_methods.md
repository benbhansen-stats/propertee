# Extract empirical estimating equations from a `glmbrob` model fit

Extract empirical estimating equations from a `glmbrob` model fit

Extract bread matrix from an `lmrob()` fit

## Usage

``` r
# S3 method for class 'glmrob'
estfun(x, ...)

# S3 method for class 'glmrob'
bread(x, ...)
```

## Arguments

- x:

  a fitted `lmrob` object

- ...:

  arguments passed to methods

## Value

matrix, estimating functions evaluated at data points and fitted
parameters

matrix, inverse Hessian of loss as evaluated at fitted parameters
