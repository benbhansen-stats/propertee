# (Internal) Fallback brute force method to locate `data` in the call stack.

We try to be intelligent about finding the appropriate data. If this
fails, we may have need for a brute force method that just loops through
frames and looks for a `data` object.

## Usage

``` r
.fallback_data_search()
```

## Value

If found, the data.
