# Identify fine strata

Identify blocks in a `StudySpecification` with exactly one treated or
one control unit of assignment.

## Usage

``` r
identify_small_blocks(spec)
```

## Arguments

- spec:

  A `StudySpecification` object.

## Value

Logical vector with length given by the number of blocks in
`StudySpecification`
