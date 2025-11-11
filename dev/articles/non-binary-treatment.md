# Non-binary Treatment Specification

## Binary Treatment

When creating a StudySpecification, handling binary treatment variables
is straightforward. If the treatment variable is either `numeric` with
only values 0/1, or is `logical`, then
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
will estimate a treatment effect of the difference between the outcome
in the treated group (`1` or `TRUE`) versus the control group (`0` or
`FALSE`).

### Missing treatment status

In all cases (binary and non-binary), missing values are allowed and any
units of assignment with missing treatment values are excluded from
models fit via
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md).

## Non-binary Treatment

However, the
[`_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
functions can take in any (reasonable) form of treatment assignment.

If the treatment variable is a `numeric` with non-binary values, it is
treated as a continuous treatment effect and `lmitt(y ~ 1, ...` will
estimate a single coefficient on treatment.

If the treatment variable is a `character`, it is treated as a
multi-level treatment variable and `lmitt(y ~ 1, ...` will estimate
treatment effects against a reference category. The reference category
is the first level defined according to [R’s comparison of
characters](https://rdrr.io/r/base/Comparison.html).

`factor` and `ordered` objects are tricky to deal with, so while a
`StudySpecification` can be created with `factor` or `ordered` treatment
variables,
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
will refuse to estimate a model unless it is also provided a
[`dichotomy`](#dichotomzing-a-non-binary-treatment) (see below).

### Dichotomzing a Non-binary Treatment

Studies may offer treatment to units at different times or provide
treatment to units in varying intensities. Researchers may be interested
in estimating treatment effects at different times or given a certain
threshold of provided treatment, however. `propertee` accommodates these
wishes by storing the time or intensity of treatment for treated units
in the
[`StudySpecification`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
then offering a `dichotomy=` argument to the weights calculation
functions
[`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)/[`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
and the assginment creation function
[`assigned()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)
A `dichotomy` is presented as a formula, where the left-hand side is a
logical statement defining inclusion in the treatment group, and the
right-hand side is a logical statement defining inclusion in the control
group. For example, if `dose` represents the intensity of a given
treatment, we could set a threshold of 200, say, mg:

``` r
dose > 200 ~ dose <= 200
```

All units of assignment with `dose` above 200 are treated units, and all
units of assignment with `dose` of 200 or below are control units.

A `.` can be used to define either group as the inverse of the other.
For example, the above dichotomy could be defined as either of

``` r
dose > 200 ~ .
. ~ dose <= 200
```

Any units of assignment not assigned to either treatment or control are
assumed to have `NA` for a treatment status and will be ignored in the
estimation of treatment effects.

``` r
dose >= 300 ~ dose <= 100
```

In this `dichotomy`, units of assignment in the range (100,300) are
ignored.

### An Example

``` r
data(simdata)
table(simdata$dose)
#> 
#>  50 100 200 250 300 
#>  10  10  10  10  10
spec1 <- rct_spec(dose ~ uoa(uoa1, uoa2), data = simdata)
summary(spec1)
#> Randomized Control Trial
#> 
#>  Structure          Variables 
#>  ---------          --------- 
#>  Treatment          dose      
#>  Unit of Assignment uoa1, uoa2
#> 
#> Number of units per Treatment group: 
#>  Txt Grp Num Units
#>       50         2
#>      100         2
#>      200         2
#>      ...          
#> 2 smaller treatment groups excluded.
#> Use `specification_table` function to view full results.
```

``` r
head(ate(spec1, data = simdata, dichotomy = dose >= 300 ~ dose <= 100))
#> [1] 1.5 1.5 1.5 1.5 0.0 0.0
```

``` r
head(assigned(spec1, data = simdata, dichotomy = dose >= 300 ~ dose <= 100))
#> [1]  0  0  0  0 NA NA
```
