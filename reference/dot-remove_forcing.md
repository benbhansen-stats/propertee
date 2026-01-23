# (Internal) Removes the forcing column entirely from a `StudySpecification`

In preparation for converting an RD `StudySpecification` to another
`StudySpecification`, this will strip the forcing variable entirely. It
is removed from the data (both `@structure` and `@column_index`), as
well as from the formula stored in `@call`.

## Usage

``` r
.remove_forcing(spec)
```

## Arguments

- spec:

  A `StudySpecification`

## Value

The `StudySpecification` without any forcing variable

## Details

Note that the output `StudySpecification` will fail a validity check
(with `validObject()`) due to an RD `StudySpecification` requiring a
forcing variable, so change the `@type` immediately.
