# (Internal) Expand treatment variable from a `StudySpecification` to a dataframe with unit of assignment information

(Internal) Expand treatment variable from a `StudySpecification` to a
dataframe with unit of assignment information

## Usage

``` r
.expand_txt(txt, data, spec)
```

## Arguments

- txt:

  A dataframe with one column corresponding to the treatment. Can be
  dichotomized or as it's stored in `spec`

- data:

  A dataframe with unit of assignment information

- spec:

  A `StudySpecification`, used to align unit of assignment information
  with `txt`
