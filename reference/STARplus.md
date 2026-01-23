# STAR participants plus nonexperimental controls

Data from Tennessee’s Project STAR study. This data frame describes
student participants in the Project STAR (Student-Teacher Achievement
Ratio) field experiment conducted in Tennessee, USA beginning in the
mid-1980s, as well as an external control group consisting of the
contemporaneous cohort of students attending a matched sample of
Tennessee schools that did not participate in the STAR experiment.
Variables are as described in Project `STAR` data documentation (see
references), with five exceptions. Three `*_at_entry` variables were
constructed as follows: `grade_at_entry` indicates the grade of
student's first participation, while `school_at_entry` and
`cond_at_entry` reflect the school ID and classroom type corresponding
to the student's grade at entry to the study. Additionally, `read_yr1`
and `math_yr1` capture a student's scaled scores on the Scholastic
Assessment Test (SAT) administered to them during their
`grade_of_entry`, i.e. their earliest available post-treatment SAT
measurements.

## Usage

``` r
STARplus
```

## Format

A `data.frame` with 13,382 rows and 56 columns.

- stdnid Student ID

- gender Student gender

- race Student race

- birthmonth Student month of birth

- birthday Student day of birth

- birthyear Student year of birth

- read_yr1 SAT reading scaled score from grade at which student entered
  the study

- math_yr1 SAT math scaled score from grade at which student entered the
  study

- gktreadss Kindergarten reading scaled score (RCT participants only)

- gktmathss Kindergarten math scaled score (RCT participants only)

- gktlistss Kindergarten listening scaled score (RCT participants only)

- gkwordskillss Kindergarten word study skills scaled score (RCT
  participants only)

- g1schid Grade 1 School ID

- g1tchid Grade 1 Teacher ID

- g1classsize Class size of Grade 1

- g1treadss Grade 1 SAT reading scaled score

- g1tmathss Grade 1 SAT math scaled score

- g1tlistss Grade 1 total listening scale score in SAT

- g1wordskillss Grade 1 word study skills scale score in SAT

- g1readbsraw Grade 1 reading raw score in Basic Skills First (BSF)
  tests

- g1mathbsraw Grade 1 math raw score in BSF

- g1readbsobjpct Grade 1 reading percent objectives mastered in BSF
  tests

- g1mathbsobjpct Grade 1 math percent objectives mastered in BSF tests

- g2schid Grade 2 School ID

- g2tchid Grade 2 Teacher ID

- g2classsize Class size of Grade 2

- g2treadss Grade 2 total reading scale score in SAT

- g2tmathss Grade 2 total math scale score in SAT

- g2tlistss Grade 2 total listening scale score in SAT

- g2wordskillss Grade 2 word study skills scale score in SAT

- g2readbsraw Grade 2 reading raw score in BSF tests

- g2mathbsraw Grade 2 math raw score in BSF test

- g2readbsobjpct Grade 2 reading percent objectives mastered in BSF
  tests

- g3schid Grade 3 School ID

- g3tchid Grade 3 Teacher ID

- g3classsize Class size of Grade 3

- g3treadss Grade 3 total reading scale score in SAT

- g3tmathss Grade 3 total math scale score in SAT

- g3langss Grade 3 total language scale score in SAT

- g3tlistss Grade 3 total listening scale score in SAT

- g3socialsciss Grade 3 social science scale score in SAT

- g3spellss Grade 3 spelling scale score in SAT

- g3vocabss Grade 3 vocabulary scale score in SAT

- g3mathcomputss Grade 3 math computation scale score in SAT

- g3mathnumconcss Grade 3 concept of numbers scale score in SAT

- g3mathapplss Grade 3 math applications scale score in SAT

- g3wordskillss Grade 3 word study skills scale score in SAT

- g3readbsraw Grade 3 reading raw score in BSF tests

- g3mathbsraw Grade 3 math raw score in BSF tests

- g3readbsobjpct Grade 3 reading percent objectives mastered in BSF
  tests

- g3mathbsobjpct Grade 3 math percent objectives mastered in BSF tests

- dob Date of birth (with NAs imputed RCT participant median)

- dobNA Dat of birth not recorded

- grade_at_entry Grade at which each student first entered the study

- school_at_entry School ID corresponding to the student's grade at
  entry into the study

- cond_at_entry Classroom type corresponding to the student's grade at
  entry into the study

## Source

[doi:10.7910/DVN/SIWH9F](https://doi.org/10.7910/DVN/SIWH9F)

## Details

Note: This dataset bears a Creative Commons Zero license (v1.0).

## References

C.M. Achilles; Helen Pate Bain; Fred Bellott; Jayne Boyd-Zaharias;
Jeremy Finn; John Folger; John Johnston; Elizabeth Word, 2008,
"Tennessee's Student Teacher Achievement Ratio (STAR) project", Harvard
Dataverse, V1, https://doi.org/10.7910/DVN/SIWH9F
UNF:3:Ji2Q+9HCCZAbw3csOdMNdA
