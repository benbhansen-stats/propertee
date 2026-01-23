# Make a dataframe that links units of assignment with clusters

Make a dataframe that links units of assignment with clusters

## Usage

``` r
.make_uoa_cluster_df(spec, cluster = NULL)
```

## Arguments

- spec:

  A `StudySpecification` object.

- cluster:

  A character vector of column names to use as clusters. Columns must
  exist in the dataframe used to create the `StudySpecification` object.
  Defaults to NULL, in which case the column names specified in the
  [`unitid()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md),
  [`unit_of_assignment()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md),
  or
  [`cluster()`](https://benbhansen-stats.github.io/propertee/reference/StudySpecificationSpecials.md)
  function in the `StudySpecification` formula will be used.

## Value

A dataframe where the number of rows coincides with the number of
distinct unit of assignment or cluster combinations (depending on
whether `cluster` is a more or less granular level than the assignment
level) and the columns correspond to the unit of assignment columns and
a "cluster" column
