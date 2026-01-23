# Adjust residuals for both-sides absorption

Adjust residuals for both-sides absorption

## Usage

``` r
block_center_residuals(x)
```

## Arguments

- x:

  a fitted `teeMod` model

## Value

the fitted `teeMod` with updated block center residuals.

## Details

This function subtracts off the block residual mean function \\\hat
\alpha(v_b, \theta)\\ for each observation from model residuals
