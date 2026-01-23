# (Internal) Checks newdata/by argument for specification accessors

(Internal) Checks newdata/by argument for specification accessors

## Usage

``` r
.specification_accessors_newdata_validate(newdata, by)
```

## Arguments

- newdata:

  newdata argument from e.g.
  [`treatment()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_extractreplace.md),
  [`blocks()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_extractreplace.md),
  etc

- by:

  from e.g.
  [`treatment()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_extractreplace.md),
  [`blocks()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_extractreplace.md),
  etc. See
  [`.check_by()`](https://benbhansen-stats.github.io/propertee/reference/dot-check_by.md)

## Value

Invisibly `TRUE`. Warns or errors as appropriate.
