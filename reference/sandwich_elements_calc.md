# (Internal) Estimate components of the sandwich covariance matrix returned by `vcov_tee()`

(Internal) Estimate components of the sandwich covariance matrix
returned by
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/reference/var_estimators.md)

## Usage

``` r
.get_a22_inverse(x, ...)

.get_a11_inverse(x)

.get_a21(x, ...)

.get_tilde_a22_inverse(x, ...)

.get_tilde_a21(x)
```

## Arguments

- x:

  a fitted `teeMod` model

- ...:

  arguments passed to `bread` method

## Value

`.get_a22_inverse()`/`.get_tilde_a22_inverse()`: A \\2\times 2\\ matrix
corresponding to an intercept and the treatment variable in the direct
adjustment model

`.get_a11_inverse()`: A \\p\times p\\ matrix where \\p\\ is the
dimension of the covariance adjustment model, including an intercept

`.get_a21()`/`.get_tilde_a21()`: A \\2\times p\\ matrix where the number
of rows are given by the intercept and the treatment variable in the
direct adjustment model, and the number of columns are given by the
dimension of the covariance adjustment model

## Details

`.get_a22_inverse()`/`.get_tilde_a22_inverse()`: \\A\_{22}^{-1}\\ is the
"bread" of the sandwich covariance matrix returned by
[`vcov_tee()`](https://benbhansen-stats.github.io/propertee/reference/var_estimators.md)
whether one has fit a prior covariance adjustment model or not.

`.get_a11_inverse()`: \\A\_{11}^{-1}\\ is the "bread" of the sandwich
covariance matrix for the covariance adjustment model. This matrix
contributes to the meat matrix of the direct adjustment sandwich
covariance matrix.

`.get_a21()`/`.get_tilde_a21()`: \\A\_{21}\\ is the gradient of the
estimating equations for the direct adjustment model taken with respect
to the covariance adjustment model parameters. This matrix is the
crossproduct of the prediction gradient for the units of observation in
\\\mathcal{Q}\\ and the model matrix of the direct adjustment model.
