# (Internal) Ensures replacement column for `StudySpecification` is a `data.frame`.

Helper function for `StudySpecification` replacers to ensure replacement
is a properly named `data.frame`

## Usage

``` r
.convert_to_data.frame(value, specification, type)
```

## Arguments

- value:

  A `vector` or `data.frame` containing a replacement.

- specification:

  A `StudySpecification`

- type:

  One of "t", "f", "u" or "b"

## Value

`data.frame` containing named column(s)

## Details

When given a replacement set of values (e.g `vector` or `matrix`), this
ensures that the replacement is a named `data.frame`.

Input `vector`: Since it cannot be named, a vector can only be used to
replace an existing component. If the existing component has more than 1
column, uses the name of the first column.

Input `matrix` or `data.frame`: If unnamed and replacing existing
component, must have no more columns than original component. (If less
columns, uses the name of the first few columns.) If named, can replace
any number of columns.
