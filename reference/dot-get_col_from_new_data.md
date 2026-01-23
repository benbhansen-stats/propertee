# (Internal) Extract specified `type` from new data set

(Internal) Extract specified `type` from new data set

## Usage

``` r
.get_col_from_new_data(
  specification,
  newdata,
  type,
  by = NULL,
  implicitBlock = FALSE,
  ...
)
```

## Arguments

- specification:

  A `StudySpecification`

- newdata:

  A `data.frame`, which may or may not be the one which was used to
  create `specification`. It must have the units of assignment
  variable(s) (though `by=` argument can be used if the name differ),
  and will appropriately merge with the `specification` the blocks,
  treatment or forcings.

- type:

  One of "t", "f", or "b".

- by:

  optional; named vector or list connecting names of unit of
  assignment/unitid/cluster variables in `specification` to unit of
  assignment/unitid/cluster variables in `data`. Names represent
  variables in the StudySpecification; values represent variables in the
  data. Only needed if variable names differ.

- implicitBlock:

  If the `StudySpecification` does not include a block, `TRUE` will
  return a constant 1 for the blocks if `type` requests it.

- ...:

  Additional arguments to
  [`merge()`](https://rdrr.io/r/base/merge.html).

## Value

The column(s) belonging to the requested `type` in
