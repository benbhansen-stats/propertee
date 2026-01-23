# Introduction to propertee

## Main Features

The **propertee** package, Prognostic Regression Offsets with
Propagation of ERrors (for Treatment Effect Estimation), facilitates
direct adjustment or standardization for experiments and observational
studies, with the option of using a separately fitted prognostic
regression model for covariance adjustment. Propertee calls for a
self-standing specification of a study’s design and treatment
allocations, using these to create weights needed to align treatment and
control samples for estimation of treatment effects averaged across a
common standard population, with standard errors that appropriately
reflect the study design. For covariance adjustment it enables
offsetting the outcome against predictions from a dedicated covariance
model, with standard error calculations propagating error as appropriate
from the covariance model.

The main workflow consists of two main steps and one optional step:

1.  Encode the study’s design, including the unit of assignment,
    treatment status of each unit of assignment, and block membership
    information as appropriate, in a `StudySpecification` object. This
    is accomplished with the
    [`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
    [`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
    or
    [`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
    functions.
2.  Optionally, fit a covariate adjustment model. This is done with any
    of a number of modeling functions external to `propertee`, for
    example [`lm()`](https://rdrr.io/r/stats/lm.html) or
    [`glm()`](https://rdrr.io/r/stats/glm.html) from the `stats`
    package, or `lmrob` from the `robustbase` package.
3.  Use [`lm()`](https://rdrr.io/r/stats/lm.html) to contrast the
    separately reweighted treatment and control groups, incorporating
    optional covariance adjustments as an `offset` to
    [`lm()`](https://rdrr.io/r/stats/lm.html). The `StudySpecification`
    is unobtrusively retained, as is error propagation from the optional
    covariance adjustment, for use in subsequent standard error
    calculations. This is done using the
    [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
    function, which either takes and enriches a fitted `lm` or, in the
    more fully supported usage, internally calls
    [`lm()`](https://rdrr.io/r/stats/lm.html) to itself create a
    directly adjusted contrast of treatment conditions.

## Example Data

The example dataset comes from the state of Tennessee’s Student-Teacher
Achievement Ratio (STAR) experiment. Students were randomly assigned to
three possible classroom conditions: small (13 to 17 students per
teacher), regular class (22 to 25 students per teacher), and
regular-with-aide class (22 to 25 students with a full-time teacher’s
aide).

``` r
data("STARplus")
table(STARplus$cond_at_entry)
#> 
#>         regular    regular+aide           small nonexperimental 
#>            4328            4249            3024            1780
```

For simplicity for this first example, we will examine a single binary
treatment - “small” classrooms versus “regular” and “regular+aide”
classrooms.

``` r
STARplus$cond_small <- STARplus$cond_at_entry == "small"
table(STARplus$cond_small)
#> 
#> FALSE  TRUE 
#> 10357  3024
```

After this basic example, we will see how **propertee** makes it easy to
handle non-binary treatment variables by introducing `dichotomy`s.

The outcome of interest is a reading score at the end of kindergarten.

``` r
summary(STARplus$g1treadss)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#>   404.0   478.0   516.0   521.5   558.0   651.0    5805
```

We first take the random assignment scheme to have been either complete
or Bernoulli randomization of students, separately by school. (Later
we’ll change this to the more realistic assumption of complete or
Bernoulli randomization, separately by school *and* year of entry into
the study.) Accordingly, experimental blocks are given by the
`school_at_entry` variable.

``` r
length(unique(STARplus$school_at_entry))
#> [1] 81
head(table(STARplus$school_at_entry))
#> 
#> 112038 123056 128068 128076 128079 130085 
#>    129     85     54     96    135    135
```

We need a unique identifier for the unit by which treatment was
allocated. STAR treatments were assigned to students, as opposed to
classrooms or other aggregates of students; accordingly we need a unique
identifier per student. The `stdntid` variable fills this role.

``` r
length(unique(STARplus$stdntid))
#> [1] 13381
head(STARplus$stdntid)
#> [1] 10000 10001 10002 10003 10004 10005
```

## A Basic Example

### Defining the `StudySpecification`

The three `_spec` functions
([`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
[`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
and
[`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md))
operate similarly. The first argument is the most important, and encodes
all the specification information through the use of a formula. The
left-hand side of the formula identifies the treatment variable. The
right-hand side of the formula consists of the following potential
pieces of information:

1.  [`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md):
    This identifies the variable(s) which indicate the units of
    assignment. This is required for all specifications. The alias
    [`uoa()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
    can be used in its place.
2.  [`block()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md):
    The identifies the variable(s) which contain block information.
    Optional.
3.  [`forcing()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md):
    In regression discontinuity specifications
    ([`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)),
    this identifies the variable(s) which contain forcing information.

To define a `StudySpecification` in our example:

``` r
spec <- obs_spec(cond_small ~ unit_of_assignment(stdntid) + block(school_at_entry),
                  data = STARplus, na.fail = FALSE)
summary(spec)
#> Observational Study
#> 
#>  Structure          Variables      
#>  ---------          ---------      
#>  Treatment          cond_small     
#>  Unit of Assignment stdntid        
#>  Block              school_at_entry
#> 
#> Number of units per Treatment group: 
#>  Txt Grp Num Units
#>    FALSE      8577
#>     TRUE      3024
```

Should additional variables be needed to identify the unit of
assignment, block, or forcing, they can be included. Had the variable
identifying the years in which participants entered the study been
available, for example, we might have identified the experimental blocks
as combinations of school and year of study entry, rather than as the
school alone. To do this we would pass
`block(year_at_entry, school_at_entry)`, rather than
`block(school_at_entry)`, in the above.

### Estimating the treatment effect

The main function for estimating treatment effects is the
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
function. It takes in three main required arguments:

1.  A formula specifying the outcome and the desired treatment effect.
2.  The data set containing the outcome information.
3.  A [`StudySpecification`](#defining-the-specification).

Note that the data set does **not** need to be the same data set which
generated the `StudySpecification`; it does however need to include the
same variables to identify the units of assignment. (If the variable
names differ, the `by=` argument can be used to link them, though we
recommend renaming to reduce the likelihood of issues.)

For example, you may have one dataset containing school-level
information, and a separate dataset containing student-level
information. Assume school is the unit of assignment. While you could of
course merge those two data-sets, **propertee** can instead use the
school-level data to define the `StudySpecification`, and the
student-level data to estimate the treatment effect.

The formula entering
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
can take on one of two forms:

``` r
y ~ 1
```

will estimate the main treatment effect on outcome variable `y`, and

``` r
y ~ x
```

will estimate subgroup specific treatment effects for each level of `x`
for the outcome `y`, if `x` is categorical. For continuous `x`, a main
effect and a treatment-`x` interaction effect is estimated.

Therefore, to estimate the treatment effect in our example, we can run:

``` r
te <- lmitt(g1treadss ~ 1, data = STARplus, specification = spec)
#> The StudySpecification object contains block-level information, but it is not used in this model. Block information is used when weights are defined via `ate()` or `ett()` or if the `absorb=TRUE` argument is passed.
summary(te)
#> 
#> Call:
#> lmitt(g1treadss ~ 1, data = STARplus, specification = spec)
#> 
#>  Treatment Effects :
#>                       Estimate Std. Error t value Pr(>|t|)    
#> cond_small.TRUE        12.6024     1.5938   7.907 3.08e-15 ***
#> g1treadss:(Intercept) 517.4628     0.7898 655.150  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Std. Error calculated via type "HC0"
```

The data includes ethnicity; we can estimate subgroup effects by
ethnicity:

``` r
te_s <- lmitt(g1treadss ~ race, data = STARplus, specification = spec)
#> The StudySpecification object contains block-level information, but it is not used in this model. Block information is used when weights are defined via `ate()` or `ett()` or if the `absorb=TRUE` argument is passed.
summary(te_s)
#> Warning: The following subgroups do not have sufficient degrees of freedom for
#> standard error estimates and will be returned as NA: racemssng
#> 
#> Call:
#> lmitt(g1treadss ~ race, data = STARplus, specification = spec)
#> 
#> Treatment Effects: (1 not defined because of singularities)
#>                           Estimate Std. Error t value Pr(>|t|)    
#> cond_small.TRUE_racewhite   9.6229     1.9779   4.865 1.17e-06 ***
#> cond_small.TRUE_raceblack  16.4164     2.3296   7.047 2.02e-12 ***
#> cond_small.TRUE_raceasian -19.3846    29.7818  -0.651    0.515    
#> cond_small.TRUE_racehispa  43.1667    50.6159   0.853    0.394    
#> cond_small.TRUE_racenatam  12.5000    37.2522   0.336    0.737    
#> cond_small.TRUE_raceother   8.0000    30.4991   0.262    0.793    
#> cond_small.TRUE_racemssng       NA         NA      NA       NA    
#> g1treadss:racewhite       530.4149     0.9873 537.252  < 2e-16 ***
#> g1treadss:raceblack       492.0263     1.0665 461.335  < 2e-16 ***
#> g1treadss:raceasian       552.3846    16.8166  32.848  < 2e-16 ***
#> g1treadss:racehispa       526.5000    25.1480  20.936  < 2e-16 ***
#> g1treadss:racenatam       490.0000    26.1639  18.728  < 2e-16 ***
#> g1treadss:raceother       547.0000    19.0245  28.752  < 2e-16 ***
#> g1treadss:racemssng       480.0000         NA      NA       NA    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Std. Error calculated via type "HC0"
```

### Including specification weights

Study specification weights can be easily included in this estimation.
**propertee** supports average treatment effect (ATE) and effect of the
treatment on the treated (ETT) weights.

To include one of the weights, simply include the `weights = "ate"` or
`weights = "ett"` argument to
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md):

``` r
lmitt(g1treadss ~ 1, data = STARplus, specification = spec, weights = "ate")
#>       cond_small.TRUE g1treadss:(Intercept) 
#>              11.35214             517.91270
lmitt(g1treadss ~ 1, data = STARplus, specification = spec, weights = "ett")
#>       cond_small.TRUE g1treadss:(Intercept) 
#>              10.88851             519.17669
```

Internally, these call the
[`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
or
[`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
functions which can be used directly.

``` r
head(ate(spec, data = STARplus))
#> [1] 1.448276 1.237013 1.433071 1.413793 1.675676 1.361345
```

When included inside
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md),
you do not need to specify any additional arguments to
[`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
or
[`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md),
enabling easy functions of weights. For example if you had some other
weight variable, say `wgt`, you could include `weights = wgt*ate()` in
the
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
call.

### Covariance Adjustment models

By itself,
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
does not allow for other covariates; e.g. something like
`lmitt(y ~ 1 + control_var,...` will fail. To adjust for covariates, a
separate covariate model should be fit. Any model which supports a
[`predict()`](https://rdrr.io/r/stats/predict.html) function should
work.

``` r
camod <- lm(g1treadss ~ gender * dob + race, data = STARplus)
```

The
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
function can be used to process the covariance adjustment model and
produce the required values; and its output can be passed as an
`offset=`.

``` r
lmitt(g1treadss ~ 1, data = STARplus, specification = spec,
      weights = "ate", offset = cov_adj(camod))
#>       cond_small.TRUE g1treadss:(Intercept) 
#>              11.03055             517.91270
```

Similarly to the weight functions,
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
attempts to locate the correct arguments (in this case, mainly the
`data=` argument) to use in the model command; while
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
does fall back to using the data which is in the covariance model, its
safer to use the `newdata=` argument if calling
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
outside of the model.

``` r
head(cov_adj(camod, newdata = STARplus))
#> [1] 519.6613 528.8581 495.1975 531.1456 500.2958 490.3070
```

Also, similarly to weights,
[`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md)
can be used in normal modeling commands as well.

``` r
lm(g1treadss ~ cond_small, data = STARplus, weights = ate(spec),
   offset = cov_adj(camod))
#> 
#> Call:
#> lm(formula = g1treadss ~ cond_small, data = STARplus, weights = ate(spec), 
#>     offset = cov_adj(camod))
#> 
#> Coefficients:
#>    (Intercept)  cond_smallTRUE  
#>         -3.266          11.031
```

### Absorbing Blocks

If fixed effects for blocks are desired, which can be absorbed away to
avoid estimating, the `absorb=TRUE` argument can be passed.

``` r
lmitt(g1treadss ~ 1, data = STARplus, specification = spec, absorb = TRUE)
#>       cond_small.TRUE g1treadss:(Intercept) 
#>               11.0594              518.8194
```
