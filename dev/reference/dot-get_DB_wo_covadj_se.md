# (Internal) Design-based variance for models without covariance adjustment

(Internal) Design-based variance for models without covariance
adjustment

## Usage

``` r
.get_DB_wo_covadj_se(x, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

## Value

design-based variance estimate of the main treatment effect estimate

## Details

Calculate bread matrix for design-based variance estimate for `teeMod`
models without covariance adjustment and without absorbed effects
