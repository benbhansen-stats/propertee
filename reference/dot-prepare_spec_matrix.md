# (Internal) Helper function for design-based meat matrix calculation

(Internal) Helper function for design-based meat matrix calculation

## Usage

``` r
.prepare_spec_matrix(x, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

## Value

a \\m \times (p+2)\\ matrix of cluster sums of design-based estimating
equations scaled by \\\sqrt{m\_{b0}m\_{b1}}/m\_{b}\\. Here \\m\\ is the
number of clusters, \\p\\ is the number of covariates used in the prior
covariance adjustment (excluding intercept)
