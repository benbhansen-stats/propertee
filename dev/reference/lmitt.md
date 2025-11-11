# Linear Model for Intention To Treat

Generates a linear model object to estimate a treatment effect, with
proper estimation of variances accounting for the study specification.

## Usage

``` r
lmitt(obj, specification, data, ...)

# S3 method for class 'formula'
lmitt(
  obj,
  specification,
  data,
  absorb = FALSE,
  offset = NULL,
  weights = NULL,
  ...
)

# S3 method for class 'lm'
lmitt(obj, specification = NULL, ...)
```

## Arguments

- obj:

  A `formula` or a `lm` object. See Details.

- specification:

  The `StudySpecification` to be used. Alternatively, a formula creating
  a specification. (Of the type of that would be passed as the first
  argument to
  [`rd_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
  [`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
  or
  [`obs_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md),
  with the difference that
  [`cluster()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md),
  [`uoa()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  and
  [`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  terms can be omitted when each row of `data` represents a distinct
  unit of assignment.) If the formula includes a
  [`forcing()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecificationSpecials.md)
  element, an RD specification is created. Otherwise an observational
  specification is created. An RCT specification must be created
  manually using
  [`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md).

- data:

  A `data.frame` such as would be passed into
  [`lm()`](https://rdrr.io/r/stats/lm.html).

- ...:

  Additional arguments passed to
  [`lm()`](https://rdrr.io/r/stats/lm.html) and other functions. An
  example of the latter is `dichotomy=`, a formula passed to
  [`assigned()`](https://benbhansen-stats.github.io/propertee/dev/reference/AssignedAliases.md)
  and, as appropriate,
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md),
  [`att()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md),
  [`atc()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  or
  [`ato()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md).
  It is used to dichotomize a non-binary treatment variable in
  `specification`. See the Details section of the
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  help page for examples.

- absorb:

  If `TRUE`, fixed effects are included for blocks identified in the
  `StudySpecification`. Excluded in `FALSE`. Default is `FALSE`. The
  estimates of these fixed effects are suppressed from the returned
  object.

- offset:

  Offset of the kind which would be passed into
  [`lm()`](https://rdrr.io/r/stats/lm.html). Ideally, this should be the
  output of
  [`cov_adj()`](https://benbhansen-stats.github.io/propertee/dev/reference/cov_adj.md).

- weights:

  Which weights should be generated? Options are `"ate"` or `"ett"`.
  Alternatively, the output of a manually run
  [`ate()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  or
  [`ett()`](https://benbhansen-stats.github.io/propertee/dev/reference/WeightCreators.md)
  can be used.

## Value

`teeMod` object (see Details)

## Details

The first argument to `lmitt()` should be a formula specifying the
outcome on the left hand side. The right hand side of the formula can be
any of the following:

- `1`: Estimates a main treatment effect.

- a subgroup variable: Estimates a treatment effect within each level of
  your subgrouping variable.

- a continuous moderator: Estimates a main treatment effect as well as a
  treatment by moderator interaction. The moderator is not automatically
  centered.

Alternatively, `obj` can be a pre-created `lm` object. No modification
is made to the formula of the object. See the help for
[`as.lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/as_lmitt.md)
for details of this conversion.

The `lmitt()` function's `subset=` argument governs the subsetting of
`data` prior to model fitting, just as with
[`lm()`](https://rdrr.io/r/stats/lm.html). Functions such as
[`rct_spec()`](https://benbhansen-stats.github.io/propertee/dev/reference/StudySpecification_objects.md)
that create `StudySpecification`s also take an optional `subset=`
argument, but its role differs from that of the `subset=` argument of
[`lm()`](https://rdrr.io/r/stats/lm.html) or `lmitt()`. The `subset=`
argument when creating a `StudySpecification` restricts the data used to
generate the `StudySpecification`, but has no direct impact on the
future [`lm()`](https://rdrr.io/r/stats/lm.html) or `lmitt()` calls
using that `StudySpecification`. (It can have an indirect impact by
excluding particular units from receiving a treatment assignment or
weight. When treatment assignments or weights are reconstructed from the
`StudySpecification`, these units will receive NAs, and will be excluded
from the [`lm()`](https://rdrr.io/r/stats/lm.html) or `lmitt()` fit
under typical `na.action` settings.)

To avoid variable name collision, the treatment variable defined in the
`specification` will have a "`.`" appended to it. For example, if you
request a main treatment effect (with a formula of `~ 1`) with a
treatment variable named "txt", you can obtain its estimate from the
returned `teeMod` object via `$coefficients["txt."]`.

`lmitt()` will produce a message if the `StudySpecification` designates
treatment assignment by block but the blocking structure appears not to
be reflected in the `weights`, nor in a block fixed effect adjustment
(via `absorb=TRUE`). While not an error, this is at odds with intended
uses of `propertee`, so `lmitt()` flags it as a potential oversight on
the part of the analyst. To disable this message, run
`options("propertee_message_on_unused_blocks" = FALSE)`.

`lmitt()` returns objects of class ‘`teeMod`’, for Treatment Effect
Estimate Model, extending the lm class to add a summary of the response
distribution under control (the coefficients of a controls-only
regression of the response on an intercept and any moderator variable).
`teeMod` objects also record the underlying `StudySpecification` and
information about any externally fitted models `mod` that may have been
used for covariance adjustment by passing `offset=cov_adj(mod)`. In the
latter case, responses are offsetted by predictions from `mod` prior to
treatment effect estimation, but estimates of the response variable
distribution under control are calculated without reference to `mod`.

The response distribution under control is also characterized when
treatment effects are estimated with block fixed effects, i.e. for
`lmitt()` with a `formula` first argument with option `absorb=TRUE`.
Here as otherwise, the supplementary coefficients describe a regression
of the response on an intercept and moderator variables, to which only
control observations contribute; but in this case the weights are
modified for this supplementary regression. The treatment effect
estimates adjusted for block fixed effects can be seen to coincide with
estimates calculated without block effect but with weights multiplied by
an additional factor specific to the combination of block and treatment
condition. For block \\s\\ containing units with weights \\w_i\\ and
binary treatment assignments \\z_i\\, define \\\hat{\pi}\_s\\ by
\\\hat{\pi}\_s\sum_sw_i=\sum_sz_iw_i\\. If \\\hat{\pi}\_s\\ is 0 or 1,
the block doesn't contribute to effect estimation and the additional
weighting factor is 0; if \\0 \< \hat{\pi}\_s \< 1\\, the additional
weighting factor is \\1 - \hat{\pi}\_s\\ for treatment group members and
\\\hat{\pi}\_s\\ for controls. When estimating a main effect only or a
main effect with continuous moderator, supplementary coefficients under
option `absorb=TRUE` reflect regressions with additional weighting
factor equal to 0 or \\\hat{\pi}\_s\\, respectively, for treatment or
control group members of block \\s\\. With a categorical moderator and
`absorb=TRUE`, this additional weighting factor determining
supplementary coefficients is calculated separately for each level
\\\ell\\ of the moderator variable, with the sums defining
\\\hat{\pi}\_{s\ell}\\ restricted not only to block \\s\\ but also to
observations with moderator equal to \\\ell\\.

## Examples

``` r
data(simdata)
spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
mod1 <- lmitt(y ~ 1, data = simdata, specification = spec, weights = "ate")
mod2 <- lmitt(y ~ as.factor(o), data = simdata, specification = spec, weights = "ate")
### observational study with treatment z assigned row-wise within blocks:
mod3 <- lmitt(y ~ 1, data=simdata, specification=z ~ block(bid), weights="att")
### regression discontinuity study with units of assignment
### given by combinations of uoa1, uoa2:
mod4 <- lmitt(y ~ 1, data = simdata,
              specification = z ~ uoa(uoa1, uoa2) + forcing(force))
```
