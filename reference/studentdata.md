# Student data

An example of data sets stored at two levels.

## Usage

``` r
schooldata

studentdata
```

## Format

Two `data.frame`s, one with school-level data (`schooldata`) including
treatment assignment and a second with student-level data
(`studentdata`). `schoolata`:

- schoolid Unique school ID variable.

- treatment Was this school in the intervention group?

- state State which the school is in.

- pct_disadvantage Percent of student body flagged as "disadvantaged".

`studentdata`:

- id Unique student ID.

- schoolid Unique school ID variable.

- grade Student's grade, 3-5.

- gpa Student GPA in prior year.

- math Standarized math score (out of 100).

An object of class `data.frame` with 8713 rows and 5 columns.

## Details

In this hypothetical data, schools were randomly assignment to treatment
status, but the unit of analysis is students. Thus the two data sets,
one encoding school information (including treatment status) and one
encoding student information (which does not include treatment status).

## Examples

``` r
soec <- obs_spec(treatment ~ uoa(schoolid), data = schooldata)

# Treatment effect
mod1 <- lmitt(math ~ 1, specification = soec, data = studentdata)

# Treatment effect by grade
mod2 <- lmitt(math ~ as.factor(grade), specification = soec, data = studentdata)
```
