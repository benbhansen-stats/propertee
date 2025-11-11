# Accessors and Replacers for `StudySpecification` objects

Allows access to the elements which define a `StudySpecification`,
enabling their extraction or replacement.

## Usage

``` r
treatment(x, newdata = NULL, dichotomy = NULL, by = NULL, ...)

# S4 method for class 'StudySpecification'
treatment(x, newdata = NULL, dichotomy = NULL, by = NULL, ...)

treatment(x) <- value

# S4 method for class 'StudySpecification'
treatment(x) <- value

units_of_assignment(x, newdata = NULL, by = NULL)

# S4 method for class 'StudySpecification'
units_of_assignment(x, newdata = NULL, by = NULL)

units_of_assignment(x) <- value

# S4 method for class 'StudySpecification'
units_of_assignment(x) <- value

clusters(x, newdata = NULL, by = NULL)

# S4 method for class 'StudySpecification'
clusters(x, newdata = NULL, by = NULL)

clusters(x) <- value

# S4 method for class 'StudySpecification'
clusters(x) <- value

unitids(x)

# S4 method for class 'StudySpecification'
unitids(x)

unitids(x) <- value

# S4 method for class 'StudySpecification'
unitids(x) <- value

blocks(x, newdata = NULL, by = NULL, ...)

# S4 method for class 'StudySpecification'
blocks(x, newdata = NULL, by = NULL, ..., implicit = FALSE)

blocks(x) <- value

# S4 method for class 'StudySpecification'
blocks(x) <- value

has_blocks(x)

forcings(x, newdata = NULL, by = NULL)

# S4 method for class 'StudySpecification'
forcings(x, newdata = NULL, by = NULL)

forcings(x) <- value

# S4 method for class 'StudySpecification'
forcings(x) <- value
```

## Arguments

- x:

  a `StudySpecification` object

- newdata:

  optional; an additional `data.frame`. If passed, and the unit of
  assignment variable is found in `newdata`, then the requested variable
  type for each unit of `newdata` is returned. See `by` argument if the
  name of the unit of assignment differs.

- dichotomy:

  optional; a formula specifying how to dichotomize a non-binary
  treatment variable. See the Details section of the
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  or
  [`att()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  help pages for information on specifying this formula

- by:

  optional; named vector or list connecting names of unit of
  assignment/unitid/cluster variables in `x` to unit of
  assignment/unitid/cluster variables in `data`. Names represent
  variables in `x`; values represent variables in `newdata`. Only needed
  if variable names differ.

- ...:

  ignored.

- value:

  replacement. Either a `vector`/`matrix` of appropriate dimension, or a
  named `data.frame` if renaming variable as well. See `Details`.

- implicit:

  Should a block-less `StudySpecification` return a constant 1 when
  extracting `blocks`?

## Value

`data.frame` containing requested variable, or an updated
`StudySpecification`. `treatment()` works slightly differently, see
`Details`.

## Details

For `treatment()`, when argument `binary` is `FALSE`, the treatment
variable passed into the `StudySpecification` is returned as a
one-column `data.frame` regardless of whether it is binary or `x` has a
`dichotomy`

If a `dichotomy` is passed, a binary one-column `data.frame` will be
returned. If not and `binary` is `TRUE`, unless the `StudySpecification`
has a binary treatment, `treatment()` will error. If `binary` is
`"ifany"`, it will return the original treatment in this case.

The one-column `data.frame` returned by `treatment()` is named as
entered in the `StudySpecification` creation, but if a `dichotomy` is
passed, the column name is `"__z"` to try and avoid any name conflicts.

For the `value` when using replacers, the replacement must have the same
number of rows as the `StudySpecification` (the same number of units of
assignment). The number of columns can differ (e.g. if the
`StudySpecification` were defined with two variable uniquely identifying
blocks, you can replace that with a single variable uniquely identifying
blocks, as long as it respects other restrictions.)

If the replacement value is a `data.frame`, the name of the columns is
used as the new variable names. If the replacement is a `matrix` or
`vector`, the original names are retained. If reducing the number of
variables (e.g., moving from two variables uniquely identifying to a
single variable), the appropriate number of variable names are retained.
If increasing the number of variables, a `data.frame` with names must be
provided.

## Examples

``` r
data(simdata)
spec <- obs_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
blocks(spec) # empty
#> data frame with 0 columns and 10 rows
blocks(spec) <- data.frame(blks = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5))
blocks(spec)
#>    blks
#> 1     1
#> 2     1
#> 3     2
#> 4     2
#> 5     3
#> 6     3
#> 7     4
#> 8     4
#> 9     5
#> 10    5
blocks(spec) <- c(5, 5, 4, 4, 3, 3, 2, 2, 1, 1)
blocks(spec) # notice that variable is not renamed
#>    blks
#> 1     5
#> 2     5
#> 3     4
#> 4     4
#> 5     3
#> 6     3
#> 7     2
#> 8     2
#> 9     1
#> 10    1
```
