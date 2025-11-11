# Add new variables to a model frame from a `teeMod` object

A variation of expand.model.frame which works for `teeMod` objects

## Usage

``` r
.expand.model.frame_teeMod(
  model,
  extras,
  envir = environment(formula(model)),
  na.expand = FALSE
)
```

## Arguments

- model:

  A `teeMod` object

- extras:

  one-sided formula or vector of character strings describing new
  variables to be added

- envir:

  an environment to evaluate things in

- na.expand:

  logical; see
  [`stats::expand.model.frame`](https://rdrr.io/r/stats/expand.model.frame.html)
  for details

## Value

A `data.frame`

## Details

When building a `teeMod` object inside
[`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md),
we do a lot of manipulation of the variables involved in the model such
that by the time the `teeMod` is produced, neither the outcome nor
predictors actually fit in the model exist in the `data` passed into the
call.

(E.g. to be specific, if a user calls
`myda <- lmitt(y ~ 1, data = mydata)`, then `model.frame(myda)` would
contain column names not found in `mydata`.)

This is a clone of
[`stats::expand.model.frame()`](https://rdrr.io/r/stats/expand.model.frame.html)
which has one addition

- after extracting the `model$call$data` from `model`, it adds columns
  from `model.frame(model)` to the object. This ensures that the
  additional variables created during
  [`lmitt()`](https://benbhansen-stats.github.io/propertee/dev/reference/lmitt.md)
  can be found.

Trivial modifications from
[`stats::expand.model.frame()`](https://rdrr.io/r/stats/expand.model.frame.html)
include ensuring `model` is a `teeMod` object, and using the `::` syntax
as appropriate.
