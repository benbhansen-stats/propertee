# (Internal) Aggregate weights and outcomes to cluster level

(Internal) Aggregate weights and outcomes to cluster level

## Usage

``` r
.aggregate_to_cluster(x, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

## Value

a list of

- a data frame of cluster weights, outcomes, treatments, and block ids;

- treatment id column name;

- block id column name

## Details

aggregate individual weights and outcomes to cluster weighted sums
