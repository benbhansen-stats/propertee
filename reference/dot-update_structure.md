# (Internal) Replaces `type` columns in `specification` with `new`

Assumes
[`.convert_to_data.frame()`](https://benbhansen-stats.github.io/propertee/reference/dot-convert_to_data.frame.md)
has already been called on `new`

## Usage

``` r
.update_structure(specification, new, type)
```

## Arguments

- specification:

  A `StudySpecification`

- new:

  A named `data.frame` with the replacement, should be the output of
  [`.convert_to_data.frame()`](https://benbhansen-stats.github.io/propertee/reference/dot-convert_to_data.frame.md).

- type:

  One of "t", "f", "u" or "b\\.

## Value

The updated `StudySpecification`
