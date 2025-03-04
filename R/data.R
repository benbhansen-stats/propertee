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

#' @title STAR data
#'
#' @description Data from Tennessee’s Project STAR study.
#'   This is a copy of \code{AER::STAR}, describing student
#'   participants in the Project STAR (Student-Teacher Achievement
#'   Ratio) study conducted in Tennessee, USA in the late 1980s.
#'   Variables are as described in the \code{AER} package's \code{STAR}
#'   help page, with the exception that the \code{row.names}
#'   of that data frame have been moved into a new first column,
#'   \code{studentid}.
#'
#' @format  A \code{data.frame} with 11,598 rows and 48 columns.
#' @keywords dataset
"STARdata"

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
#' @description Schools matched into pairs or triplets approximating the subset
#' of participating schools in Pane et al. (2014) from Michigan
#' @references Pane, John F., et al. "Effectiveness of cognitive tutor algebra I at scale."
#' \emph{Educational Evaluation and Policy Analysis} 36.2 (2014): 127-144.
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

