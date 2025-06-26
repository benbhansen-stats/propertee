# this dataframe is used for vignettes/mi_vignette
data_temp <- data.frame(
  schoolid = c("6305000291", "6316006171", "6320001204", "6324009415",
               "6326003242", "6327000710", "6328002123", "6329004340",
               "6307005976", "6315004226", "6314002317", "6318000385",
               "6301004608", "6326005819"),
  blk = c("E", "B", "D", "C", "F", "F", "E", "B", "A", "A", "C", "D", "E", "F"),
  z = c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0)
)

source("../vignettes/mi_vignette/get_and_clean_external_data.R")

merged_data <- merge(analysis1data, data_temp, by = "schoolid")
michigan_school_pairs <- merged_data[, c("schoolid", "blk", "z", "MALE_G11_PERC", "FEMALE_G11_PERC", "AM_G11_PERC", "ASIAN_G11_PERC", "HISP_G11_PERC", 
                            "BLACK_G11_PERC", "WHITE_G11_PERC", "PACIFIC_G11_PERC", "TR_G11_PERC","G11")]
colnames(michigan_school_pairs)  <-
    tolower(colnames(michigan_school_pairs))

usethis::use_data(michigan_school_pairs, overwrite = TRUE)
