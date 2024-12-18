library(AER)
data(STAR)

STARdata <-
    data.frame(studentid=rownames(STAR), STAR, row.names=NULL)

usethis::use_data(STARdata, overwrite = TRUE)
