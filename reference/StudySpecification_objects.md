# Generates a `StudySpecification` object with the given specifications.

Generate a randomized control treatment StudySpecification
(`rct_spec()`), or an observational StudySpecification (`obs_spec()`),
or a regression discontinuity StudySpecification (`rd_spec()`).

## Usage

``` r
rct_spec(formula, data, subset = NULL, na.fail = TRUE)

rd_spec(formula, data, subset = NULL, na.fail = TRUE)

obs_spec(formula, data, subset = NULL, na.fail = TRUE)

rct_specification(formula, data, subset = NULL, na.fail = TRUE)

rd_specification(formula, data, subset = NULL, na.fail = TRUE)

obs_specification(formula, data, subset = NULL, na.fail = TRUE)

obsstudy_spec(formula, data, subset = NULL, na.fail = TRUE)

obsstudy_specification(formula, data, subset = NULL, na.fail = TRUE)
```

## Arguments

- formula:

  a `formula` defining the `StudySpecification` components. See
  `Details` for specification.

- data:

  the data set from which to build the StudySpecification. Note that
  this data need not be the same as used to estimate the treatment
  effect; rather the `data` passed should contain information about the
  units of treatment assignment (as opposed to the units of analysis).

- subset:

  optional, subset the data before creating the `StudySpecification`
  object

- na.fail:

  If `TRUE` (default), any missing data found in the variables specified
  in `formula` (excluding treatment) will trigger an error. If `FALSE`,
  non-complete cases will be dropped before the creation of the
  `StudySpecification`

## Value

a `StudySpecification` object of the requested type for use in further
analysis.

## Details

The formula should include exactly one
[`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md)
to identify the units of assignment (one or more variables). (`uoa`,
`cluster`, or `unitid` are synonyms for `unit_of_assignment`; the choice
of which has no impact on the analysis. See below for a limited
exception in which the `unit_of_assignment` specification may be
omitted.) If defining an `rd_spec`, the formula must also include a
[`forcing()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md)
entry. The formula may optionally include a
[`block()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md)
as well. Each of these can take in multiple variables, e.g. to pass both
a household ID and individual ID as unit of assignment, use
`uoa(hhid, iid)` and not `uoa(hhid) + uoa(iid)`.

The treatment variable passed into the left-hand side of `formula` can
either be `logical`, `numeric`, or `character`. If it is anything else,
it attempts conversion to one of those types (for example, `factor` and
`ordered` are converted to `numeric` if the levels are `numeric`,
otherwise to `character`). If the treatment is not `logical` or
`numeric` with only values 0 and 1, in order to generate weights with
[`ate()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
or
[`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md),
the `dichotomy` argument must be used in those functions to identify the
treatment and control groups. See
[`ett()`](https://benbhansen-stats.github.io/propertee/reference/WeightCreators.md)
for more details on specifying a `dichotomy`.

There are a few aliases for each version.

If the formula excludes a
[`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md),
data merges are performed on row order. Such formulas can also be passed
as the specification argument to lmitt(), and that is their primary
intended use case. It is recommended that each formula argument passed
to \*\_specification() include a
[`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md),
[`uoa()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md)
or
[`cluster()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md)
term identifying the key variable(s) with which `StudySpecification`
data is to be merged with analysis data. Exceptions to this rule will be
met with a warning. To disable the warning, run
`options("propertee_warn_on_no_unit_of_assignment" = FALSE)`.

The units of assignment, blocks, and forcing variables must be `numeric`
or `character`. If they are otherwise, an attempt is made to cast them
into `character`.

## Examples

``` r
data(simdata)
spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2) + block(bid),
                  data = simdata)

data(schooldata)
spec <- obs_spec(treatment ~ unit_of_assignment(schoolid) + block(state),
                  data = schooldata)
```
