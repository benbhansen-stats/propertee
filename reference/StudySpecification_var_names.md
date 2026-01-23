# Extract Variable Names from `StudySpecification`

Methods to extract the variable names to the elements of the structure
of the `StudySpecification` (e.g. treatment, unit of analysis, etc)

## Usage

``` r
var_table(specification, compress = TRUE, report_all = FALSE)

var_names(specification, type, implicitBlocks = FALSE)
```

## Arguments

- specification:

  a `StudySpecification` object

- compress:

  should multiple variables be compressed into a comma-separated string?
  Default `TRUE`. If `FALSE`, multiple columns can be created instead.

- report_all:

  should we report all possible structures even if they don't exist in
  the `StudySpecification`? Default `FALSE`.

- type:

  one of "t", "u", "b", "f"; for "treatment", "unit_of_assignment",
  "block", and "forcing" respectively

- implicitBlocks:

  If the `StudySpecification` is created without blocks, setting this to
  `TRUE` will return "`.blocks_internal`" as the variable name
  corresponding to the blocks.

## Value

`var_table` returns the requested table. `var_names` returns a vector of
variable names.

## Details

When `compress` is `TRUE`, the result will always have two columns. When
`FALSE`, the result will have number of columns equal to the largest
number of variables in a particular role, plus one. E.g., a call such as
`rct_spec(z ~ unitid(a, b, c, d) ...` will have 4+1=5 columns in the
output matrix with `compress = FALSE`.

When `report_all` is `TRUE`, the matrix is guaranteed to have 3 rows
(when the `specification` is an RCT or Obs) or 4 rows (when the
`specification` is a RD), with empty variable entries as appropriate.
When `FALSE`, the matrix will have minimum 2 rows (treatment and unit of
assignment/unitid/cluster), with additional rows for blocks and forcing
if included in the `StudySpecification`.

## Examples

``` r
spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
var_table(spec)
#>      Structure            Variables   
#> [1,] "Treatment"          "z"         
#> [2,] "Unit of Assignment" "uoa1, uoa2"
#> [3,] "Block"              "bid"       
var_table(spec, compress = FALSE)
#>      Structure            Variable 1 Variable 2
#> [1,] "Treatment"          "z"        NA        
#> [2,] "Unit of Assignment" "uoa1"     "uoa2"    
#> [3,] "Block"              "bid"      NA        
var_names(spec, "t")
#> [1] "z"
var_names(spec, "u")
#> [1] "uoa1" "uoa2"
var_names(spec, "b")
#> [1] "bid"
```
