# (Internal) Create a new `StudySpecification` object.

Helper function to create a new `StudySpecification`. Called internally
from
[`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
[`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
or
[`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md).

## Usage

``` r
.new_StudySpecification(
  form,
  data,
  type,
  subset = NULL,
  call = NULL,
  na.fail = TRUE,
  called_from_lmitt = FALSE
)
```

## Arguments

- form:

  Formula to create StudySpecification, see help for `rcr_spec()`,
  [`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  or
  [`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
  for details.

- data:

  The data set

- type:

  One of "RCT", "RD", or "Obs"

- subset:

  Any subset information

- call:

  The call generating the `StudySpecification`.

- na.fail:

  Should it error on NA's (`TRUE`) or remove them (`FALSE`)?

- called_from_lmitt:

  Logical; was this called inside
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md),
  or was it called from `*_spec()` (default).

## Value

A new StudySpecification object
