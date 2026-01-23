# (Internal) Replace standard errors for moderator effect estimates with insufficient degrees of freedom with `NA`

(Internal) Replace standard errors for moderator effect estimates with
insufficient degrees of freedom with `NA`

## Usage

``` r
.check_df_moderator_estimates(
  vmat,
  model,
  cluster_cols,
  model_data = quote(data),
  envir = environment(formula(model))
)
```

## Arguments

- vmat:

  output of `.vcov_XXX()` called with input to `model` argument below as
  the first argument

- model:

  a fitted `teeMod` model

- cluster_cols:

  a character vector indicating the column(s) defining cluster ID's

- model_data:

  dataframe or name of dataframe used to fit `model`

- envir:

  environment to get `model_data` from if `model_data` has class `name`

## Value

A variance-covariance matrix with NA's for estimates lacking sufficient
degrees of freedom
