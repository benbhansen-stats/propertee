# (Internal) Return ID's used to order observations in the direct adjustment sample

(Internal) Return ID's used to order observations in the direct
adjustment sample

## Usage

``` r
.sanitize_Q_ids(x, id_col = NULL, ...)
```

## Arguments

- x:

  a fitted `teeMod` model

- id_col:

  character vector; optional. Specifies column(s) whose ID's will be
  returned. The column must exist in the data that created the
  `StudySpecification` object. Default is NULL, in which case unit of
  assignment columns indicated in the specification will be used to
  generate ID's.

- ...:

  arguments passed to methods

## Value

A vector with length equal to the number of units of observation in the
direct adjustment sample
