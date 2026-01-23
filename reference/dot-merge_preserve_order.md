# (Internal) Merge `data.frame`s ensuring order of first `data.frame` is maintained

(Internal) Merge `data.frame`s ensuring order of first `data.frame` is
maintained

## Usage

``` r
.merge_preserve_order(x, ...)
```

## Arguments

- x:

  `data.frame` whose ordering is to be maintained

- ...:

  Additional arguments to
  [`merge()`](https://rdrr.io/r/base/merge.html), particularly a second
  `data.frame` and a `by=` argument.

## Value

Merged `data.frame` with the same ordering as `x`.
