# (Internal) Helper function for design-based meat matrix calculation

(Internal) Helper function for design-based meat matrix calculation

## Usage

``` r
.cov_mat_est(XXz, bidz)
```

## Value

estimated upper and lower bounds of covariance matrix of estimating
function vectors under either treatment or control

## Details

Diagonal elements are estimated by sample variances Off-diagonal
elements are estimated using the Young's elementary inequality
