# Valid Weights

These track and access the weights we support as well as any aliases.

## Usage

``` r
.validWeights

.isValidWeightTarget(target)

.isValidWeightAlias(alias)

.listValidWeightTargets()

.listValidWeightAliases()
```

## Format

An object of class `list` of length 2.

## Arguments

- target:

  String

- alias:

  String

## Value

Logical for `.isValid*`, and a string for `.listValid*`.

## Details

"target" refers to weight calculations we support.

"alias" refers to all possible names. Every "target" has at least one
"alias" (itself) and may have more.

`.isValidWeightTarget()` and `.isValidWeightAlias()` identify whether a
given input (from a user) is a value weighting name.

`.listValidWeightTargets()` and `.listValidWeightAliases()` are for
returning nicely formatted strings for messages to users.

IMPORTANT: Adding new aliases MUST correspond to new functions defined
in weights_exported.R, with possible adjustments to calculations in
weights_internal.R.
