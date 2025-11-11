# Intervention data from a pair-matched study of schools in Michigan

Michigan high schools, with a plausible cluster RCT

## Usage

``` r
michigan_school_pairs
```

## Format

A `data.frame` with 14 rows and 13 columns.

- schoolid school id

- blk block

- z treatment variable

- MALE_G11_PERC percentage of G11 male students

- FEMALE_G11_PERC percentage of G11 female students

- AM_G11_PERC percentage of G11 American Indian/Alaska Native students

- ASIAN_G11_PERC percentage of G11 Asian students

- HISP_G11_PERC percentage of G11 Hispanic students

- BLACK_G11_PERC percentage of G11 Black students

- WHITE_G11_PERC percentage of G11 White students

- PACIFIC_G11_PERC percentage of G11 Hawaiian Native/Pacific Islander
  students

- TR_G11_PERC percentage of G11 Two or More Races students

- G11 Number of G11 students

## Details

Grade 11 demographics for all Michigan high schools in 2013, with mock
block and treatment assignments for 14 high schools within a large
county in the metro Detroit area. These schools were selected for this
demonstration based on their similarity to the 14 high schools from an
adjacent Michigan county that participated in the Pane et al (2013)
study. As a result, they serve as an example of what one might expect to
find as the state-specific school-level subsample in a multi-state
paired cluster randomized trial featuring random assignment at the
school level.

The mock experimental schools were selected by optimal matching of
experimental schools to adjacent county schools, with substitute schools
grouped into the same pairs or triples (‘fine strata’) as were their
experimental counterparts. The original pairs and triples had been
selected to reduce variation in baseline variables predictive of
outcomes, and the blocking structure the substitute sample inherits may
be expected to do this as well. The treatment/control distinction is
also inherited from the experimental sample, but there is of course no
treatment effect within the mock experiment.

The selection of mock experimental schools was based on both demographic
and student achievement variables, but the present data frame includes
only the demographic variables (as sourced from the Common Core of Data
\[CCD; U.S. Department of Education\]). School average outcomes in
student test scores are available separately, from Michigan's Center for
Education Performance Information. See the vignette ‘Real-data
demonstration with a finely stratified cluster RCT and a broader
administrative database’, available on the package website.

## References

Pane, John F., et al. "Effectiveness of cognitive tutor algebra I at
scale." *Educational Evaluation and Policy Analysis* 36.2 (2014):
127-144.

U.S. Department of Education. Public Elementary/Secondary School
Universe Survey Data, v.2a. Institute of Education Sciences, National
Center for Education Statistics.

## Examples

``` r
data(michigan_school_pairs)
mi_spec <- rct_spec(z ~ uoa(schoolid)+block(blk),
data=michigan_school_pairs)
mi_spec
#> Randomized Control Trial
#> 
#>  Structure          Variables
#>  ---------          ---------
#>  Treatment          z        
#>  Unit of Assignment schoolid 
#>  Block              blk      
#> 
table(is.na(michigan_school_pairs$blk))
#> 
#> FALSE  TRUE 
#>    14   680 
specification_table(mi_spec, "block", "treatment")
#>          blocks
#> treatment A B C D E F
#>         0 1 1 1 1 2 1
#>         1 1 1 1 1 1 2
```
