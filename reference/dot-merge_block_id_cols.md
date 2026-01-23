# (Internal) Merge multiple block IDs

(Internal) Merge multiple block IDs

## Usage

``` r
.merge_block_id_cols(df, ids)
```

## Arguments

- df:

  a data frame

- ids:

  a vector of block IDs, column names of df

## Value

a data frame with a column that contains unique block number IDs

## Details

merge multiple block ID columns by the value combinations and store the
new block ID in the column `ids[1]`
