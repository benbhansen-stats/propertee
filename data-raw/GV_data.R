# Load necessary libraries for data manipulation
library(dplyr)   

set.seed(123)

# Read in the dataset
GV <- read.csv(file="GreenVavreck_PolAnalysis_2008_PA_Replication.csv", head=TRUE, sep=",")

#' Green-Vavreck Dataset
#'
#' A subset of data from the '04 Rock-the-Vote data
#'
#' @format ## `GV`
#' A data frame with 23,869 rows and 49 columns:
#' \describe{
#'   \item{xage4}{Age of participant in decimal format}
#'   \item{tout1}{A binary variable indicating past turnout behavior for each participant}
#'   \item{syscode}{Sytem cable code}
#'   \item{treat}{A binary variable denoting treatment assignment for participants}
#'   \item{matched_pairs}{A numeric indicator for the strata or matched pair group to which a cable system belongs}
#'   \item{syspopall}{The total population size of the cable system to which the participant belongs}
#'   \item{ageint}{Integer value of participant's age}
#'   \item{wt}{The total sample size for the syscode cluster}
#' }
#' @source <https://isps.yale.edu/research/data/d005>
"GV"

# Filter the data for matched pairs 1-3 and select relevant columns
GV_temp <- GV[GV$matched_pairs %in% c(1, 2, 3), c("xage04", "tout1", "syscode", "treat", "matched_pairs", "syspopall", "ageint", "wt")]

# Randomly sample 10% based on "syscode"
GV_sample <- GV_temp %>%
  group_by(syscode) %>%
  sample_frac(0.1) %>%
  ungroup()

# Save the final dataset to a CSV file (if needed)
# write.csv(GV_sample, "GV_sampled_data.csv", row.names = FALSE)

# first few rows of the sampled data for review
head(GV_sample)

usethis::use_data(GV_sample, overwrite = TRUE)
