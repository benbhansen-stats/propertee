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
#' @details See the "Regression Discontinuity Designs" vignette on the
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
#' @description Tennessee’s Student-Teacher Achievement Ratio (STAR) data set.
#'   This is a copy of the data \code{AER::STAR}.
#'
#' @format  A \code{data.frame} with 11,598 rows and 47 columns.
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
#' des <- obs_design(treatment ~ uoa(schoolid), data = schooldata)
#'
#' # Treatment effect
#' mod1 <- lmitt(math ~ 1, design = des, data = studentdata)
#'
#' # Treatment effect by grade
#' mod2 <- lmitt(math ~ as.factor(grade), design = des, data = studentdata)
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
#' @format  A \code{data.frame} with 14 rows and 3 columns.
#' @keywords dataset
"michigan_school_pairs"