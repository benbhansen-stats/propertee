# (Internal) A few checks to ensure `by=` is valid

Thie ensures that the `by=` argument is of the proper type, is named,
and consists of only unique entries.

## Usage

``` r
.check_by(by)
```

## Arguments

- by:

  named vector or list connecting names of unit of
  assignment/unitid/cluster variables in `specification` to unit of
  assignment/unitid/cluster variables in `data`. Names represent
  variables in the StudySpecification; values represent variables in the
  data.

## Value

`NULL` if no errors are found
