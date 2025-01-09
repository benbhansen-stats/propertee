# this dataframe is used for vignettes/mi_vignette
data_temp <- data.frame(
  schoolid = c("6305000291", "6316006171", "6320001204", "6324009415",
               "6326003242", "6327000710", "6328002123", "6329004340",
               "6307005976", "6315004226", "6314002317", "6318000385",
               "6301004608", "6326005819"),
  blk = c("E", "B", "D", "C", "F", "F", "E", "B", "A", "A", "C", "D", "E", "F"),
  z = c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0)
)

source("../propertee/vignettes/mi_vignette/get_and_clean_external_data.R")

merged_data <- merge(analysis1data, data_temp, by = "schoolid")
michigan_school_pairs <- merged_data[, c("schoolid", "blk", "z", "male_g11_perc", "female_g11_perc", "am_g11_perc", "asian_g11_perc", "hisp_g11_perc", 
                            "black_g11_perc", "white_g11_perc", "pacific_g11_perc", "tr_g11_perc","g11")]

usethis::use_data(michigan_school_pairs, overwrite = TRUE)
