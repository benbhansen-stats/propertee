# (Internal) Perform checks on formula for creation of StudySpecification.

Checks performed:

- Ensure presence of no more than one of
  [`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md),
  [`cluster()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  or
  [`unitid()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md).

- Disallow multiple
  [`block()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  or multiple
  [`forcing()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  terms.

- Disallow
  [`forcing()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  unless in RDD.

## Usage

``` r
.check_spec_formula(form, allow_forcing = FALSE)
```

## Arguments

- form:

  A formula passed to
  [`.new_StudySpecification()`](https://benbhansen-stats.github.io/propertee/dev/reference/dot-new_StudySpecification.md)

- allow_forcing:

  Binary whether
  [`forcing()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  is allowed (`TRUE` for RDD, `FALSE` for RCT and Obs).

## Value

`TRUE` if all checks pass, otherwise errors.
