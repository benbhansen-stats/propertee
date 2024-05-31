# this dataframe is used for vignettes/mi_vignette
paireddata <- data.frame(
  schoolid = c("6305000291", "6316006171", "6320001204", "6324009415",
               "6326003242", "6327000710", "6328002123", "6329004340",
               "6307005976", "6315004226", "6314002317", "6318000385",
               "6301004608", "6326005819"),
  blk = c("E", "B", "D", "C", "F", "F", "E", "B", "A", "A", "C", "D", "E", "F"),
  z = c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0)
)

usethis::use_data(paireddata, overwrite = TRUE)
