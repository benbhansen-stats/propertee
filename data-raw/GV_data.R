# Load necessary libraries for data manipulation
library(dplyr)   

set.seed(123)

# Read in the dataset
GV <- read.csv(file="GreenVavreck_PolAnalysis_2008_PA_Replication.csv", head=TRUE, sep=",")

# Filter the data for matched pairs 1-3 and select relevant columns
GV_temp <- GV[GV$matched_pairs %in% c(1, 2, 3), c("xage04", "tout1", "syscode", "treat", "matched_pairs", "syspopall", "ageint", "wt")]

# Randomly sample 10% based on "syscode"
GV_data <- GV_temp %>%
  group_by(syscode) %>%
  sample_frac(0.1) %>%
  ungroup()

# first few rows of the sampled data for review
head(GV_data)

# Save the data to the package's data directory
usethis::use_data(GV_data, overwrite = TRUE)

