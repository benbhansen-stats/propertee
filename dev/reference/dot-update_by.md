# (Internal) Use `by` to update `StudySpecification` with new variable names

Helper function used to update the variable names of a
`StudySpecification` when user passes a `by=` argument to align variable
names between data sets.

## Usage

``` r
.update_by(specification, data, by)
```

## Arguments

- specification:

  A `StudySpecification`

- data:

  `Data set`

- by:

  named vector or list connecting names of unit of
  assignment/unitid/cluster variables in `specification` to unit of
  assignment/unitid/cluster variables in `data`. Names represent
  variables in the StudySpecification; values represent variables in the
  data.

## Value

A `StudySpecification` with updated variable names
