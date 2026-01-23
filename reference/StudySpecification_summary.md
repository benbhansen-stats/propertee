# Summarizing `StudySpecification` objects

[`summary()`](https://rdrr.io/r/base/summary.html) method for class
`StudySpecification`.

## Usage

``` r
# S3 method for class 'StudySpecification'
summary(object, ..., treatment_binary = TRUE)

# S3 method for class 'summary.StudySpecification'
print(x, ..., max_unit_print = 3)
```

## Arguments

- object:

  `StudySpecification` object, usually a result of a call to
  [`rct_spec()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_objects.md),
  [`obs_spec()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_objects.md),
  or
  [`rd_spec()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_objects.md).

- ...:

  Ignored

- treatment_binary:

  Should the treatment be dichotomized if `object` contains a
  `dichotomy`? Ignored if `object` does not contain a `dichotomy`.

- x:

  `summary.StudySpecification` object, usually as a result of a call to
  `summary.StudySpecification()`

- max_unit_print:

  Maximum number of treatment levels to print in treatment table

## Value

The `StudySpecification` or `summary.StudySpecification`object,
invisibly
