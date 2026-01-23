# `WeightedStudySpecification` subsetting

Provides functionality to subset the weights of a
`WeightedStudySpecification` object.

## Usage

``` r
# S4 method for class 'WeightedStudySpecification'
subset(x, subset)

# S4 method for class 'WeightedStudySpecification'
x[i]
```

## Arguments

- x:

  `WeightedStudySpecification` object

- subset:

  Logical vector identifying values to keep or drop

- i:

  indices specifying elements to extract or replace. See
  [`help("[")`](https://rdrr.io/r/base/Extract.html) for further
  details.

## Value

A `WeightedStudySpecification` object which is a subsetted version of
`x`.
