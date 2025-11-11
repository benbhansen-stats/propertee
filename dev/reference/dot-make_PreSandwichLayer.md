# (Internal) Get covariance adjustments and their gradient with respect to covariance adjustment model coefficients

`.make_PreSandwichLayer()` takes a fitted covariance adjustment model
passed to the `model` argument and generates adjustments to outcomes for
observations in the `newdata` argument. It also evaluates the gradient
of the adjustments taken with respect to the coefficients at the
coefficient estimates.

## Usage

``` r
.make_PreSandwichLayer(model, newdata = NULL, ...)
```

## Arguments

- model:

  a fitted model to use for generating covariance adjustment values

- newdata:

  a dataframe with columns called for in `model`

- ...:

  additional arguments to pass on to `model.frame` and `model.matrix`.
  These cannot include `na.action`, `xlev`, or `contrasts.arg`: the
  former is fixed to be `na.pass`, while the latter two are provided by
  elements of the `model` argument.

## Value

A `PreSandwichLayer` object
