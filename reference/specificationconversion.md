# Convert `StudySpecification` between types

Convert a `StudySpecification` between a observational study, a
randomized control trial, and a regression discontinuity (created from
`obs_spec`, `rct_spec` and `rd_spec` respectively).

## Usage

``` r
as_rct_spec(StudySpecification, ..., loseforcing = FALSE)

as_obs_spec(StudySpecification, ..., loseforcing = FALSE)

as_rd_spec(StudySpecification, data, ..., forcing)
```

## Arguments

- StudySpecification:

  a `StudySpecification` to convert

- ...:

  Ignored.

- loseforcing:

  converting from RD to another `StudySpecification` type will error to
  avoid losing the forcing variable. Setting `loseforcing = TRUE` allows
  the conversion to automatically drop the forcing variable. Default
  `FALSE`.

- data:

  converting to an RD requires adding a `forcing` variable, which
  requires access to the original data.

- forcing:

  converting to an RD requires adding a `forcing` variable. This should
  be entered as a formula which would be passed to `update`, e.g.
  `forcing = . ~ . + forcing(forcevar)`.

## Value

`StudySpecification` of the updated type

## Examples

``` r
spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
spec
#> Randomized Control Trial
#> 
#>  Structure          Variables 
#>  ---------          --------- 
#>  Treatment          z         
#>  Unit of Assignment uoa1, uoa2
#> 
as_obs_spec(spec)
#> Observational Study
#> 
#>  Structure          Variables 
#>  ---------          --------- 
#>  Treatment          z         
#>  Unit of Assignment uoa1, uoa2
#> 
as_rd_spec(spec, simdata, forcing = ~ . + forcing(force))
#> Regression Discontinuity StudySpecification
#> 
#>  Structure          Variables 
#>  ---------          --------- 
#>  Treatment          z         
#>  Unit of Assignment uoa1, uoa2
#>  Forcing            force     
#> 
spec2 <- rd_spec(o ~ uoa(uoa1, uoa2) + forcing(force), data = simdata)
spec2
#> Regression Discontinuity StudySpecification
#> 
#>  Structure          Variables 
#>  ---------          --------- 
#>  Treatment          o         
#>  Unit of Assignment uoa1, uoa2
#>  Forcing            force     
#> 
# as_rct_spec(spec2) # this will produce an error
as_rct_spec(spec2, loseforcing = TRUE)
#> Randomized Control Trial
#> 
#>  Structure          Variables 
#>  ---------          --------- 
#>  Treatment          o         
#>  Unit of Assignment uoa1, uoa2
#> 
```
