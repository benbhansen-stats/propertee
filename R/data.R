#' @title Simulated data
#'
#' @description Simulated data to use with the \bold{propertee} package with unit
#'   of assignment level treatment assignment
#'
#' @format A data.frame with 100 rows and 7 columns.
#' \itemize{
#'   \item \code{uoa1} First level unit of assignment ID
#'   \item \code{uoa2} Second level unit of assignment ID
#'   \item \code{bid} Block ID
#'   \item \code{force} Forcing variable
#'   \item \code{z} Binary treatment indicator
#'   \item \code{o} 4-level ordered treatment variable
#'   \item \code{dose} Dose treatment variable
#'   \item \code{x} Some predictor
#'   \item \code{y} Some outcome
#' }
#' @keywords dataset
"simdata"

#' @title Synthethic Regression Discontinuity Data
#'
#' @description The data for this example were randomly simulated using the
#'   \bold{synthpop} package in R based on data originally collected by Lindo,
#'   Sanders, and Oreopoulos (2010).
#'
#' @details See the "Regression Discontinuity StudySpecifications" vignette on the
#'   \code{propertee} website for more details on the original data, a link to
#'   the code used to generate this synthethic data, and a detailed example.
#'
#' @format A \code{data.frame} with 40,403 rows and 11 columns.
#' \itemize{
#'   \item \code{R}
#'   \item \code{lhsgrade_pct}
#'   \item \code{nextGPA}
#'   \item \code{probation_year1}
#'   \item \code{totcredits_year1}
#'   \item \code{male}
#'   \item \code{loc_campus1}
#'   \item \code{loc_campus2}
#'   \item \code{bpl_north_america}
#'   \item \code{english}
#'   \item \code{age_at_entry}
#' }
#' @keywords dataset
"lsoSynth"

#' @title STAR participants plus nonexperimental controls
#'
#' @description Data from Tennessee’s Project STAR study.
#'   This data frame describes student participants in the Project
#'   STAR (Student-Teacher Achievement Ratio) field experiment
#'   conducted in Tennessee, USA beginning in the mid-1980s, as well
#'   as an external control group consisting of the contemporaneous
#'   cohort of students attending a matched sample of Tennessee
#'   schools that did not participate in the STAR experiment.
#'   Variables are as described in Project \code{STAR} data
#'   documentation (see references), with five exceptions.
#'   Three \code{*_at_entry} variables were constructed
#'   as follows: \code{grade_at_entry} indicates the grade of student's
#'   first participation, while \code{school_at_entry} and
#'   \code{cond_at_entry} reflect the school ID and classroom type corresponding
#'   to the student's grade at entry to the study. Additionally, \code{read_yr1}
#'   and \code{math_yr1} capture a student's scaled scores on the Scholastic
#'   Assessment Test (SAT) administered to them during their \code{grade_of_entry}, 
#'   i.e. their earliest available post-treatment SAT measurements.
#'
#' @format  A \code{data.frame} with 13,382 rows and 56 columns.
#' \itemize{
#'    \item {stdnid} Student ID
#'    \item {gender} Student gender
#'    \item {race} Student race
#'    \item {birthmonth} Student month of birth
#'    \item {birthday} Student day of birth
#'    \item {birthyear} Student year of birth
#'    \item {read_yr1} SAT reading scaled score from grade at which
#'          student entered the study
#'    \item {math_yr1} SAT math scaled score from grade at which 
#'          student entered the study
#'    \item {gktreadss} Kindergarten reading scaled score (RCT participants only)
#'    \item {gktmathss} Kindergarten math scaled score (RCT participants only)
#'    \item {gktlistss} Kindergarten listening scaled score (RCT participants only)
#'    \item {gkwordskillss} Kindergarten word study skills scaled score (RCT participants only)
#'    \item {g1schid} Grade 1 School ID
#'    \item {g1tchid} Grade 1 Teacher ID
#'    \item {g1classsize} Class size of Grade 1
#'    \item {g1treadss} Grade 1 SAT reading scaled score
#'    \item {g1tmathss} Grade 1 SAT math scaled score
#'    \item {g1tlistss} Grade 1 total listening scale score in SAT
#'    \item {g1wordskillss} Grade 1 word study skills scale score in SAT
#'    \item {g1readbsraw} Grade 1 reading raw score in Basic Skills First (BSF) tests
#'    \item {g1mathbsraw} Grade 1 math raw score in BSF
#'    \item {g1readbsobjpct} Grade 1 reading percent objectives mastered in BSF tests
#'    \item {g1mathbsobjpct} Grade 1 math percent objectives mastered in BSF tests
#'    \item {g2schid} Grade 2 School ID
#'    \item {g2tchid} Grade 2 Teacher ID
#'    \item {g2classsize} Class size of Grade 2
#'    \item {g2treadss} Grade 2 total reading scale score in SAT
#'    \item {g2tmathss} Grade 2 total math scale score in SAT
#'    \item {g2tlistss} Grade 2 total listening scale score in SAT
#'    \item {g2wordskillss} Grade 2 word study skills scale score in SAT
#'    \item {g2readbsraw} Grade 2 reading raw score in BSF tests
#'    \item {g2mathbsraw} Grade 2 math raw score in BSF test
#'    \item {g2readbsobjpct} Grade 2 reading percent objectives mastered in BSF tests
#'    \item {g3schid} Grade 3 School ID
#'    \item {g3tchid} Grade 3 Teacher ID
#'    \item {g3classsize} Class size of Grade 3
#'    \item {g3treadss} Grade 3 total reading scale score in SAT
#'    \item {g3tmathss} Grade 3 total math scale score in SAT
#'    \item {g3langss} Grade 3 total language scale score in SAT
#'    \item {g3tlistss} Grade 3 total listening scale score in SAT
#'    \item {g3socialsciss} Grade 3 social science scale score in SAT
#'    \item {g3spellss} Grade 3 spelling scale score in SAT
#'    \item {g3vocabss} Grade 3 vocabulary scale score in SAT
#'    \item {g3mathcomputss} Grade 3 math computation scale score in SAT
#'    \item {g3mathnumconcss} Grade 3 concept of numbers scale score in SAT
#'    \item {g3mathapplss} Grade 3 math applications scale score in SAT
#'    \item {g3wordskillss} Grade 3 word study skills scale score in SAT
#'    \item {g3readbsraw} Grade 3 reading raw score in BSF tests
#'    \item {g3mathbsraw} Grade 3 math raw score in BSF tests
#'    \item {g3readbsobjpct} Grade 3 reading percent objectives mastered in BSF tests
#'    \item {g3mathbsobjpct} Grade 3 math percent objectives mastered in BSF tests
#'    \item {dob} Date of birth (with NAs imputed RCT participant median)
#'    \item {dobNA} Dat of birth not recorded
#'    \item {grade_at_entry} Grade at which each student first entered the study
#'    \item {school_at_entry} School ID corresponding to the student's grade
#'          at entry into the study
#'    \item {cond_at_entry} Classroom type corresponding to the student's grade
#'          at entry into the study
#' }
#' @details Note: This dataset bears a Creative Commons Zero license (v1.0).
#'
#' @references C.M. Achilles; Helen Pate Bain; Fred Bellott; Jayne Boyd-Zaharias;
#' Jeremy Finn; John Folger; John Johnston; Elizabeth Word, 2008,
#' "Tennessee's Student Teacher Achievement Ratio (STAR) project",
#' Harvard Dataverse, V1, https://doi.org/10.7910/DVN/SIWH9F UNF:3:Ji2Q+9HCCZAbw3csOdMNdA
#'
#' @source \doi{doi:10.7910/DVN/SIWH9F}
#' @keywords dataset
"STARplus"

#' @title Student data
#'
#' @description An example of data sets stored at two levels.
#'
#' @details In this hypothetical data, schools were randomly assignment to
#'   treatment status, but the unit of analysis is students. Thus the two data
#'   sets, one encoding school information (including treatment status) and one
#'   encoding student information (which does not include treatment status).
#'
#' @rdname studentdata
#' @format Two \code{data.frame}s, one with school-level data (\code{schooldata})
#'   including treatment assignment and a second with student-level data
#'   (\code{studentdata}).
#' \code{schoolata}:
#' \itemize{
#'   \item {schoolid} Unique school ID variable.
#'   \item {treatment} Was this school in the intervention group?
#'   \item {state} State which the school is in.
#'   \item {pct}_disadvantage Percent of student body flagged as "disadvantaged".
#' }
#' \code{studentdata}:
#' \itemize{
#'   \item {id} Unique student ID.
#'   \item {schoolid} Unique school ID variable.
#'   \item {grade} Student's grade, 3-5.
#'   \item {gpa} Student GPA in prior year.
#'   \item {math} Standarized math score (out of 100).
#' }
#' @examples
#' soec <- obs_spec(treatment ~ uoa(schoolid), data = schooldata)
#'
#' # Treatment effect
#' mod1 <- lmitt(math ~ 1, specification = soec, data = studentdata)
#'
#' # Treatment effect by grade
#' mod2 <- lmitt(math ~ as.factor(grade), specification = soec, data = studentdata)
#' @keywords dataset
"schooldata"

#' @rdname studentdata
"studentdata"

#' @title Intervention data from a pair-matched study of schools in Michigan
#'
#' @description Michigan high schools, with a plausible cluster RCT
#'
#' @details
#'
#' Grade 11 demographics for all Michigan high schools in 2013,
#' with mock block and treatment assignments for 14 high schools within
#' a large county in the metro Detroit area.
#' These schools were selected for this demonstration based on their
#' similarity to the 14 high schools from an adjacent Michigan
#' county that participated in the Pane et al (2013) study.  As a result,
#' they serve as an example of what one might expect to find as the
#' state-specific school-level subsample in a multi-state
#' paired cluster randomized trial featuring random assignment at
#' the school level.
#'
#' The mock experimental schools were selected by optimal matching of
#' experimental schools to adjacent county schools, with substitute
#' schools grouped into the same pairs or triples (\sQuote{fine
#' strata}) as were their experimental counterparts.  The original
#' pairs and triples had been selected to reduce variation in baseline
#' variables predictive of outcomes, and the blocking structure the
#' substitute sample inherits may be expected to do this as well.  The
#' treatment/control distinction is also inherited from the
#' experimental sample, but there is of course no treatment effect
#' within the mock experiment.
#'
#' The selection of mock experimental schools was based on both demographic
#' and student achievement variables, but the present data frame
#' includes only the demographic variables (as sourced from the
#' Common Core of Data \[CCD; U.S. Department of Education\]).  School
#' average outcomes in student test scores are available separately, from
#' Michigan's Center for Education Performance Information.  See the
#' vignette \sQuote{Real-data demonstration with a finely stratified
#' cluster RCT and a broader administrative database}, available on
#' the package website.
#' @references Pane, John F., et al.
#'     "Effectiveness of cognitive tutor algebra I at scale."
#'     \emph{Educational Evaluation and Policy Analysis} 36.2 (2014):
#'     127-144.
#'
#' U.S. Department of Education. Public Elementary/Secondary School
#' Universe Survey Data, v.2a. Institute of Education Sciences,
#' National Center for Education Statistics.
#'
#' @format  A \code{data.frame} with 14 rows and 13 columns.
#' \itemize{
#'  \item{schoolid} school id
#'  \item{blk} block
#'  \item{z} treatment variable
#'  \item{MALE_G11_PERC} percentage of G11 male students
#'  \item{FEMALE_G11_PERC} percentage of G11 female students
#'  \item{AM_G11_PERC} percentage of G11 American Indian/Alaska Native students
#'  \item{ASIAN_G11_PERC} percentage of G11 Asian students
#'  \item{HISP_G11_PERC} percentage of G11 Hispanic students
#'  \item{BLACK_G11_PERC} percentage of G11 Black students
#'  \item{WHITE_G11_PERC} percentage of G11 White students
#'  \item{PACIFIC_G11_PERC} percentage of G11 Hawaiian Native/Pacific Islander students
#'  \item{TR_G11_PERC} percentage of G11 Two or More Races students
#'  \item{G11} Number of G11 students
#' }
#' @keywords dataset
#' @examples
#' data(michigan_school_pairs)
#' mi_spec <- rct_spec(z ~ uoa(schoolid)+block(blk),
#' data=michigan_school_pairs)
#' mi_spec
#' table(is.na(michigan_school_pairs$blk))
#' specification_table(mi_spec, "block", "treatment")
"michigan_school_pairs"

#' @title Cluster-randomized experiment data on voter turnout in cable system markets
#'
#' @description This dataset is a toy example derived from a cluster-randomized
#' field experiment that evaluates the effect of \dQuote{Rock the Vote} TV advertisements
#' on voter turnout rate. The original study included 23,869 first-time voters
#' across 85 cable television markets in 12 states. These markets were grouped
#' into matched sets based on their past voter turnout rates and then randomly
#' assigned to either a treatment or control condition. This toy dataset is
#' constructed by randomly sampling 10% of individuals from selected cable
#' television markets in the original dataset.
#'
#' @format A \code{data.frame} with 248 rows and 7 columns.
#' \itemize{
#'   \item {age} Age of participant
#'   \item {vote_04} Outcome variable indicating whether participant voted
#'   \item {tv_company} Cable system serving participant's residential area
#'   \item {treatment} Binary variable denoting treatment assignment
#'   \item {pairs} A numeric indicator for the strata or matched pair group to which a cable system belongs (1-3)
#'   \item {population_size} Total population size of residential area served by cable system
#'   \item {sample_size} Number of individuals sampled from the cable system cluster
#' }
#'
#' @details The original dataset was drawn from a randomized controlled trial in
#' which 85 cable system areas were first grouped into 40 matched sets based on
#' historical voter turnout. Within each matched set, one cable system area was
#' randomly assigned to the treatment condition, while the others served
#' as controls.
#'
#' This toy dataset includes a subset of the original replication data,
#' specifically individuals from matched sets 1–3, which encompass 7 of the 85
#' cable system areas. Within these selected clusters, a 10% random sample of
#' individuals was taken.
#'
#' The fuller Green-Vavreck dataset that this derives from bears a Creative
#' Commons BY-NC-ND license (v3.0) and is housed in Yale University's Institution
#' for Social and Policy Studies (ID: D005).
#'
#' @source <https://isps.yale.edu/research/data/d005>
#'
#' @references Green, Donald P. & Lynn Vavreck (2008) "Analysis of Cluster-Randomized
#' Experiments: A Comparison of Alternative Estimation Approaches." Political Analysis 16(2):138-152.
#'
#' @keywords dataset
"GV_data"
