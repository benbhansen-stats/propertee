# (Internal) Updates `spec@call`'s formula with the currently defined variable names.

Helper function to update the `call` with the appropriate variable names
after they've been modified. Called within `StudySpecification`
replacers.

## Usage

``` r
.update_call_formula(specification)
```

## Arguments

- specification:

  A `StudySpecification`

## Value

An updated formula

## Details

It's return should be stuck into the specification via
`spec@call$formula <- .update_call_formula(spec)`
