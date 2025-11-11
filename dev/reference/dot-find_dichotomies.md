# (Internal) Find `dichotomy` formulas in the call stack

(Internal) Find `dichotomy` formulas in the call stack

## Usage

``` r
.find_dichotomies()
```

## Value

A list where elements are either formulas or NULL, depending on whether
a `dichotomy` argument was found in
[`lmitt.formula()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
or its `weights` argument

## Details

`.find_dichotomies()` searches for
[`lmitt.formula()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
calls and their `weights` arguments for any `dichotomy` arguments.
