# (Internal) Rename cluster/unitid/uoa in a formula to unit_of_assignment for internal consistency

Internally, we always refer to uoa/cluster/unitid as
"unit_of_assignment"

## Usage

``` r
.update_form_to_unit_of_assignment(form)
```

## Arguments

- form:

  A formula passed to
  [`.new_StudySpecification()`](https://benbhansen-stats.github.io/propertee/dev/reference/dot-new_StudySpecification.md)

## Value

The formula with "cluster"/"unitid"/"uoa" replace with
"unit_of_assignment"
