# `StudySpecification` Structure Information

Obtaining a `data.frame` which encodes the specification information.

## Usage

``` r
get_structure(specification)

# S4 method for class 'StudySpecificationStructure'
show(object)
```

## Arguments

- specification:

  a `StudySpecification` object

- object:

  a `StudySpecificationStructure` object, typically the output of
  `get_structure`

## Value

A `StudySpecificationStructure` object containing the structure of the
`specification` as a `data.frame`.

## Examples

``` r
data(simdata)
spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
get_structure(spec)
#>    z uoa1 uoa2 bid
#> 1  0    1    1   1
#> 2  0    1    2   1
#> 3  0    2    1   1
#> 4  1    2    2   1
#> 5  0    3    1   2
#> 6  0    3    2   2
#> 7  1    4    1   2
#> 8  0    4    2   3
#> 9  1    5    1   3
#> 10 1    5    2   3
```
