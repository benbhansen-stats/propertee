# Special terms in `StudySpecification` creation formula

These are special functions used only in the definition of
`StudySpecification` objects. They identify the units of assignment,
blocks and forcing variables. They should never be used outside of the
`formula` argument to `obs_spec`, `rct_spec`, or `rd_spec`.

## Usage

``` r
unit_of_assignment(...)

unitid(...)

cluster(...)

uoa(...)

block(...)

forcing(...)
```

## Arguments

- ...:

  any number of variables of the same length.

## Value

the variables with appropriate labels. No use outside of their inclusion
in the `formula` argument to `obs_spec`, `rct_spec`, or `rd_spec`

## Details

These functions have no use outside of the formula in creating a
`StudySpecification`.

`unit_of_assignment`, `uoa`, `cluster` and `unitid` are synonyms; you
must include one and only one in each `StudySpecification`. The choice
of which to use will have no impact on any analysis, only on some output
and the name of the stored element in the `StudySpecification`.
Accessors/ replacers (`units_of_assignment`, `unitids`, `clusters`)
respect the choice made at the point of creation of the
`StudySpecification`, and only the appropriate function will work.

See `rct_spec`, `obs_spec`, or `rd_spec` for examples of their usage.
