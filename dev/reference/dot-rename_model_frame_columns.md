# (Internal) Rename columns to strip function calls

After calling
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) on the
formula input to
[`.new_StudySpecification()`](https://benbhansen-stats.github.io/propertee/dev/reference/dot-new_StudySpecification.md),
the names of the columns will include function names, e.g.
"block(blockvar)". This function strips all these.

## Usage

``` r
.rename_model_frame_columns(modframe)
```

## Arguments

- modframe:

  A `data.frame`.

## Value

The `data.frame` with function calls removed
