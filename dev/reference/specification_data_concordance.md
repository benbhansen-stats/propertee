# Check for variable agreement within units of assignment

Useful for debugging purposes to ensure that there is concordance
between variables in the `StudySpecification` and data.

## Usage

``` r
specification_data_concordance(
  specification,
  data,
  by = NULL,
  warn_on_nonexistence = TRUE
)
```

## Arguments

- specification:

  a `StudySpecification` object

- data:

  a new data set, presumably not the same used to create
  `specification`.

- by:

  optional; named vector or list connecting names of variables in
  `specification` to variables in `data`. Names represent variables in
  `specification`; values represent variables in `data`. Only needed if
  variable names differ.

- warn_on_nonexistence:

  default `TRUE`. If a variable does not exist in `data`, should this be
  flagged? If `FALSE`, silently move on if a variable doesn't exist in
  `data`.

## Value

invisibly `TRUE` if no warnings are produced, `FALSE` if any warnings
are produced.

## Details

Consider the following scenario: A `StudySpecification` is generated
from some dataset, "data1", which includes a block variable "b1". Within
each unique unit of assignment/unitid/cluster of "data1", it must be the
case that "b1" is constant. (Otherwise the creation of the
`StudySpecification` will fail.)

Next, a model is fit which includes weights generated from the
`StudySpecification`, but on dataset "data2". In "data2", the block
variable "b1" also exists, but due to some issue with data cleaning,
does not agree with "b1" in "data1".

This could cause errors, either directly (via actual error messages) or
simply produce nonsense results. `specification_data_concordance()` is
specificationed to help debug these scenarios by providing information
on whether variables in both the data used in the creation of
`specification` ("data1" in the above example) and some new dataset,
`data`, ("data2" in the above example) have any inconsistencies.
