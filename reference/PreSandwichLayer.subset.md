# `PreSandwichLayer` and `SandwichLayer` subsetting

Return subset of a `PreSandwichLayer` or `SandwichLayer` which meets
conditions.

## Usage

``` r
# S4 method for class 'PreSandwichLayer'
subset(x, subset)

# S4 method for class 'PreSandwichLayer'
x[i]
```

## Arguments

- x:

  `PreSandwichLayer` or `SandwichLayer` object

- subset:

  Logical vector identifying values to keep or drop

- i:

  indices specifying elements to extract or replace. See
  [`help("[")`](https://rdrr.io/r/base/Extract.html) for further
  details.

## Value

`x` subset by `subset` or `i`
