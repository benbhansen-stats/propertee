# Instruct `cov_adj()` to find a default reference value for columns in a covariance adjustment model

Instruct
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md)
to find a default reference value for columns in a covariance adjustment
model

## Usage

``` r
default()
```

## Value

empty list that inherits from class `set_to_reference_default`

## Details

When set as the entry of a list passed to the `set_to_reference`
argument of
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md),
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/reference/cov_adj.md)
uses the following replacement values depending on the column type:

- factor/character: minimum level

- numeric: minimum value

- logical: FALSE
