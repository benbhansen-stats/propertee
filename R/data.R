#' @title Simulated data with unit of assignment level treatment assignment
#' @format A data.frame with 100 rows and 7 columns.
#' \itemize{
#'   \item cid1 First level unit of assignment ID
#'   \item cid2 Second level unit of assignment ID
#'   \item bid Block ID
#'   \item force Forcing variable
#'   \item z Binary treatment indicator
#'   \item o 4-level ordered treatment variable
#'   \item dose Dose treatment variable
#'   \item x Some predictor
#'   \item y Some outcome
#' }
#' @keywords dataset
"simdata"

#' @title Fake data
#' @format A data.frame with ?? rows and ?? columns.
#' \itemize{
#'   \item ???
#' }
#' @keywords dataset
"lsoSynth"

#' @title STAR data
#' @format The STAR data set
#' \itemize{
#'   \item ???
#' }
#' @keywords dataset
"STARdata"

#' @title Student data
#' @rdname studentdata
#' @format Two data.frames, one with school-level data (\code{schooldata})
#'   including treatment assignment and a second with student-level data
#'   (\code{studentdata}).
#' \code{schoolata}:
#' \itemize{
#'   \item schoolid Unique school ID variable.
#'   \item treatment Was this school in the intervention group?
#'   \item state State which the school is in.
#'   \item pct_disadvantage Percent of student body flagged as "disadvantaged".
#' }
#' \code{studentdata}:
#' \itemize{
#'   \item id Unique student ID.
#'   \item schoolid Unique school ID variable.
#'   \item grade Student's grade, 3-5.
#'   \item gpa Student GPA in prior year.
#'   \item math Standarized math score (out of 100).
#' }
#' @examples
#' des <- obs_design(treatment ~ uoa(schoolid), data = schooldata)
#'
#' # Treatment effect
#' lmitt(math ~ 1, design = des, data = studentdata)
#'
#' # Treatment effect by grade
#' lmitt(math ~ as.factor(grade), design = des, data = studentdata)
#' @keywords dataset
"schooldata"

#' @rdname studentdata
"studentdata"
