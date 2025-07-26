##+ setup, warning=FALSE, message=FALSE
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(dplyr)
### Optional: To run the following lines, please navigate to https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/10766
### and download STAR_Students.RData and Comparison_Students.RData in the data-raw directory
load("data-raw/STAR_Students.RData")
experimentsample  <- x
load("data-raw/Comparison_Students.RData")
extrasample  <- x
rm(x)

##' Both data sets appear to have been generated from files with
##' trailing line misrepresented as an observation with all NAs. 
all(is.na(experimentsample[nrow(experimentsample), ]))
all(is.na(extrasample[nrow(extrasample), ]))

##' Fixing...
experimentsample  <- experimentsample[-nrow(experimentsample), ]
extrasample <- extrasample[-nrow(extrasample), ]

##' Ensure commensurability of variables across the two
##' data frames
##+
levels(experimentsample$gender)  <- levels(extrasample$gender)  <-
    list("male"="MALE", "female"="FEMALE", "mssng"=NA)
experimentsample[is.na(experimentsample$gender), "gender"]  <- "mssng"
extrasample[is.na(extrasample$gender), "gender"]  <- "mssng"

experimentsample$birthmonth  <- as.integer(experimentsample$birthmonth)
levels(experimentsample$race)  <- list("white"="WHITE",
                                   "black"="BLACK",
                                   "asian"="ASIAN",
                                   "hispa"="HISPANIC",
                                   "natam"="NATIVE AMERICAN",
                                   "other"="OTHER",
                                   "mssng"=NA)
experimentsample[is.na(experimentsample$race), "race"] <- "mssng"
levels(extrasample$race)  <- list("white"="WHITE",
                                   "black"="BLACK",
                                  "asian"="ASIAN",
                                  "hispa"="HISPANIC",
                                  "natam"="NATIVE A",
                                  "other"="OTHER",
                                  "mssng"=c("MISSING", NA))
extrasample[is.na(extrasample$race), "race"] <- "mssng"

##' touching up levels of classroom type variables...
##+
ctvars <- paste0("g", c("k", 1:3), "classtype")
experimentsample[ctvars] |>
    sapply({\(x) all.equal(levels(x),
                           levels(experimentsample[[ctvars[1]]])
                           ) }
           ) |>
    all() |> stopifnot()
newlevs <- setNames(levels(experimentsample$gkclasstype),
                    levels(experimentsample$gkclasstype))
newlevs <- gsub(" CLASS", "", newlevs)
newlevs <- gsub(" ", "", newlevs)
newlevs <- tolower(newlevs)
newlevs <- sort(newlevs) #just by chance, gives desired levels ordering
newlevs <- c(newlevs,
             "nonexperimental"="nonexperimental"# label for external controls
             )
for (ctvar in ctvars)
{
    stopifnot(setequal(names(newlevs)[1:3], levels(experimentsample[[ctvar]])))
    experimentsample[[ctvar]] <-
        ordered(experimentsample[[ctvar]], levels=names(newlevs),
                labels=newlevs)
               }

##' We'll create the treatment assignment variable presently,
##' but for now let's record null values of it for the nonexperimentals
##+
extrasample$cond_at_entry <-
    rep("nonexperimental", nrow(extrasample)) |>
    ordered(levels=newlevs)

##' The gk* variables are available only for experimentals, but
##' we'll NA version of some them for the extrasample
##+
Filter( \(nm) substr(nm,1,2)=="gk", colnames(experimentsample))
Filter( \(nm) substr(nm,1,2)=="gk", colnames(extrasample))
extrasample$"gktreadss" <- extrasample$"gktmathss" <-
    extrasample$"gktlistss" <- extrasample$"gkwordskillss" <-
        rep(NA_real_, nrow(extrasample))
            
##' creating date of birth variables...
##+
experimentsample$dob <-
    with(experimentsample,
         paste(birthyear, birthmonth, birthday, sep="-")
         ) |> as.Date()
extrasample$dob <-
    with(extrasample,
         paste(birthyear, birthmonth, birthday, sep="-")
         ) |> as.Date()

dob_NA_fill_value  <- median(experimentsample$dob, na.rm=TRUE)
experimentsample$dobNA  <- is.na(experimentsample$dob)
experimentsample[experimentsample$dobNA, "dob"]  <-
    dob_NA_fill_value
extrasample$dobNA  <- is.na(extrasample$dob)
extrasample[extrasample$dobNA, "dob"] <-
    dob_NA_fill_value

##' Identifying grade and school in which each experimental participant
##' entered the study (removing the 1
##' student having NA class type in each of grades K-3)...
##+
experimentsample$grade_at_entry  <-
    cbind(!is.na(subset(experimentsample, select=c(gkclasstype:g3classtype))), 1) |>
    max.col(ties.method="first")
sum(experimentsample$grade_at_entry>=4)
experimentsample <- subset(experimentsample, grade_at_entry<=4)
experimentsample$school_at_entry <-
    as.matrix(subset(experimentsample, select=c(gkschid,g1schid,g2schid,g3schid))
              )[1L:nrow(experimentsample) + nrow(experimentsample)*(experimentsample$grade_at_entry-1L)]
extrasample$school_at_entry <- rep(NA_integer_, nrow(extrasample))

experimentsample$grade_at_entry <- c("k", "1", "2", "3")[experimentsample$grade_at_entry]
experimentsample$grade_at_entry <-
    ordered( experimentsample$grade_at_entry,
            levels=c("k", "1", "2", "3") )
extrasample$grade_at_entry <-
    rep(NA, nrow(extrasample)) |>
    ordered(levels=c("k", "1", "2", "3") )

##' Assembling treatment assignment variable...
##+
experimentsample$cond_at_entry  <- experimentsample$gkclasstype
for (gl in as.character(1:3))
{
    sset  <- experimentsample$grade_at_entry==gl
    sset  <- !is.na(sset) & sset
    experimentsample[sset, "cond_at_entry"] <-
        experimentsample[sset, paste0("g", gl, "classtype")]
}
with(experimentsample, table(cond_at_entry, grade_at_entry))
##' Combining experimental & supplemental data:
##+ echo=TRUE
common_cols  <- intersect(colnames(experimentsample), colnames(extrasample))
STARplus <- rbind(experimentsample[common_cols],
                  extrasample[common_cols])

##' Create a new variable read_at_entry_p1 and move column to front
##+
STARplus$read_yr1 <- with(STARplus, ifelse(grade_at_entry == "k", gktreadss,
                                    ifelse(grade_at_entry == "1", g1treadss,
                                    ifelse(grade_at_entry == "2", g2treadss,
                                    ifelse(grade_at_entry == "3", g3treadss, NA))))
                          )
STARplus <- STARplus %>%
  relocate(read_yr1, .after = 6)

STARplus$math_yr1 <- with(STARplus, ifelse(grade_at_entry == "k", gktmathss,
                                    ifelse(grade_at_entry == "1", g1tmathss,
                                    ifelse(grade_at_entry == "2", g2tmathss,
                                    ifelse(grade_at_entry == "3", g3tmathss, NA))))
                          )
STARplus <- STARplus %>%
  relocate(math_yr1, .after = 7)

usethis::use_data(STARplus, overwrite = TRUE)
