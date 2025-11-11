# Obtain Treatment from StudySpecification

When passing a `lm` object to
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md),
extract and use the treatment variable specified in the
`StudySpecification`.

## Usage

``` r
assigned(specification = NULL, data = NULL, dichotomy = NULL)

adopters(specification = NULL, data = NULL, dichotomy = NULL)

a.(specification = NULL, data = NULL, dichotomy = NULL)

z.(specification = NULL, data = NULL, dichotomy = NULL)
```

## Arguments

- specification:

  Optional `StudySpecification`. If the `StudySpecification` can't be
  identified in the model (usually because neither weights
  ([`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  or
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md))
  nor a covariate adjustment model
  ([`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md))
  are found), the `StudySpecification` can be passed diretly.

- data:

  Optional data set. By default `assigned()` will attempt to identify
  the appropriate data, if this fails (or you want to overwrite it), you
  can pass the data here.

- dichotomy:

  optional; a formula defining the dichotomy of the treatment variable
  if it isn't already `0`/`1`. See details for more information. If
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  or
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  is called within a
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  call that specifies a `dichotomy` argument, that `dichotomy` will be
  used if the argument here has not been specified.

## Value

The treatment variable to be placed in the regression formula.

## Details

When passing a `lm` object to
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md),
the treatment variable in the `formula` passed to
[`lm()`](https://rdrr.io/r/stats/lm.html) needs to be identifiable.
Rather than placing the treatment variable directly in the `formula`,
use one of these functions, to allow
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
to identify the treatment variable.

To keep the formula in the [`lm()`](https://rdrr.io/r/stats/lm.html)
call concise, instead of passing `specification` and `data` arguments to
these functions, one can pass a `WeightedStudySpecification` object to
the `weights` argument of the [`lm()`](https://rdrr.io/r/stats/lm.html)
call or a `SandwichLayer` object to the `offset` argument.

Alternatively, you can pass the `specification` and `data` arguments.

While `assigned()` can be used in any situation, it is most useful for
scenarios where the treatment variable is non-binary and the
`StudySpecification` contains a `Dichotomy`. For example, say `q` is a
3-level ordinal treatment variable, and the binary comparison of
interest is captured in `dichotomy = q == 3 ~ q < 3`. If you were to fit
a model including `q` as a predictor, e.g. `lm(y ~ q, ...)`, `lm` would
treat `q` as the full ordinal variable. On the other hand, by calling
`lm(y ~ assigned(), weights = ate(spec), ...)`, `assigned()` will
generate the appropriate binary variable to allow estimation of
treatment effects.

If called outside of a model call and without a `data` argument, this
will extract the treatment from the `specification`. If this is the
goal, the
[`treatment()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_extractreplace.md)
function is better suited for this purpose.

## Examples

``` r
data(simdata)
spec <- obs_spec(z ~ uoa(uoa1, uoa2), data = simdata)
mod <- lm(y ~ assigned(), data = simdata, weights = ate(spec))
lmittmod <- lmitt(mod)
summary(lmittmod, vcov.type = "CR0")
#> 
#> Call:
#> lmitt(mod)
#> 
#>  Coefficients :
#>               Estimate Std. Error t value Pr(>|t|)
#> (Intercept)    -0.0090     0.1734  -0.052    0.960
#> assigned()      0.1205     0.2301   0.524    0.613
#> y:(Intercept)  -0.0090     0.1734  -0.052    0.960
#> Std. Error calculated via type "CR0"
#> 
```
