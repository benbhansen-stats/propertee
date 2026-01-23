# Generate Direct Adjusted Weights for Treatment Effect Estimation

These should primarily be used inside models. See Details.

## Usage

``` r
ett(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)

att(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)

ate(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)

etc(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)

atc(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)

ato(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)

olw(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)

owt(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)

pwt(specification = NULL, dichotomy = NULL, by = NULL, data = NULL)
```

## Arguments

- specification:

  optional; a `StudySpecification` object created by one of
  [`rct_spec()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_objects.md),
  [`rd_spec()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_objects.md),
  or
  [`obs_spec()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecification_objects.md).

- dichotomy:

  optional; a formula defining the dichotomy of the treatment variable
  if it isn't already `0`/`1`. See details for more information. If
  `ett()` or `ate()` is called within a
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/reference/lmitt.md)
  call that specifies a `dichotomy` argument, that `dichotomy` will be
  used if the argument here has not been specified.

- by:

  optional; named vector or list connecting names of unit of assignment/
  variables in `specification` to unit of assignment/unitid/cluster
  variables in `data`. Names represent variables in the
  StudySpecification; values represent variables in the data. Only
  needed if variable names differ.

- data:

  optional; the data for the analysis to be performed on. May be
  excluded if these functions are included as the `weights` argument of
  a model.

## Value

a `WeightedStudySpecification` object, which is a vector of numeric
weights

## Details

These functions should primarily be used in the `weight` argument of
[`lmitt()`](https://benbhansen-stats.github.io/propertee/reference/lmitt.md)
or[`lm()`](https://rdrr.io/r/stats/lm.html). All arguments are optional
if used within those functions. If used on their own, `specification`
and `data` must be provided.

- `ate` - Average treatment effect. Aliases: `ate()`.

- `ett` - Effect of treatment on the treated. Aliases: `ett()`, `att()`.

- `etc` - Effect of treatment on controls. Aliases: `etc()`, `atc()`.

- `ato` - Overlap-weighted average effect. Aliases: `ato()`, `olw`,
  `owt`, `pwt`.

In a `StudySpecification` with `block`s, the weights are generated as a
function of the ratio of the number of treated units in a block versus
the total number of units in a block.

In any blocks where that ratio is 0 or 1 (that is, all units in the
block have the same treatment status), the weights will be 0. In effect
this removes from the target population any block in which there is no
basis for estimating either means under treatment or means under
control.

If block is missing for a given observation, a weight of 0 is applied.

A `dichotomy` is specified by a `formula` consisting of a conditional
statement on both the left-hand side (identifying treatment levels
associated with "treatment") and the right hand side (identifying
treatment levels associated with "control"). For example, if your
treatment variable was called `dose` and doses above 250 are considered
treatment, you might write:

`ate(..., dichotomy = dose > 250 ~ dose <= 250`

The period (`.`) can be used to assign all other units of assignment.
For example, we could have written the same treatment regime as either

`etc(..., dichotomy = dose > 250 ~ .`

or

`olw(..., dichotomy = . ~ dose <= 250`

The `dichotomy` formula supports Relational Operators (see
[Comparison](https://rdrr.io/r/base/Comparison.html)), Logical Operators
(see [Logic](https://rdrr.io/r/base/Logic.html)), and `%in%` (see
[`match()`](https://rdrr.io/r/base/match.html)).

The conditionals need not assign all values of treatment to control or
treatment, for example, `dose > 300 ~ dose < 200` does not assign
`200 <= dose <= 300` to either treatment or control. This would be
equivalent to manually generating a binary variable with `NA` whenever
`dose` is between 200 and 300. Standard errors will reflect the sizes of
the comparison groups specified by the `dichotomy`.

Tim Lycurgus contributed code for the computation of weights. The
‘overlap weight’ concept is due to Li, Morgan and Zaslavsky (2018),
although the current implementation differs from that discussed in their
paper in that it avoids estimated propensity scores.

## References

Li, Fan, Kari Lock Morgan, and Alan M. Zaslavsky. "Balancing covariates
via propensity score weighting." Journal of the American Statistical
Association 113, no. 521 (2018): 390-400.

## Examples

``` r
data(simdata)
spec <- rct_spec(z ~ unit_of_assignment(uoa1, uoa2), data = simdata)
summary(lmitt(y ~ 1, data = simdata, specification = spec, weights = ate()), vcov.type = "CR0")
#> 
#> Call:
#> lmitt(y ~ 1, data = simdata, specification = spec, weights = ate())
#> 
#>  Treatment Effects :
#>               Estimate Std. Error t value Pr(>|t|)
#> z.              0.1205     0.2301   0.524    0.613
#> y:(Intercept)  -0.0090     0.1734  -0.052    0.960
#> Std. Error calculated via type "CR0"
#> 
```
