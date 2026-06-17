# Where to pull the data from
extdataURLs  <-
    list(CCD="https://nces.ed.gov/ccd/data/zip/sc132a_txt.zip",
         MME="https://www.michigan.gov/cepi/-/media/Project/Websites/cepi/MiSchoolData/historical/Historical_Assessments/2011-2014MME.zip")

# pads numerical values with leading zeros for fixed character length
zero_lpad <- function(vec, char_len) {
  vapply(
    vec,
    function(x) paste(c(rep(0, char_len - nchar(x)), x), collapse = ""),
    character(1L)
  )
}

CCD_ID_COLS <- c("STID", "LEANM", "SEASCH", "SCHNAM", "CONUM", "CONAME")
CCD_CAT_COLS <- c("TITLEISTAT", "MAGNET", "CHARTR", "TYPE")
# school grade from kindergarten to G12
CCD_ENROLLMENT_COLS <- c("KG", paste0("G", zero_lpad(seq(1, 12), 2)))
RACEETH_COLS <- c("AM", "ASIAN", "HISP", "BLACK", "WHITE", "PACIFIC", "TR")
RACEGENDER_PREFIXES <- c("AM", "AS", "HI", "BL", "WH", "HP", "TR")
# Total Free and Reduced-Price Lunch Eligible Student
FRL_COL <- "TOTFRL"

# recode variable into 0 or 1
recode_binary_variable <- function(col) {
  if (inherits(col, "character")) {
    nonmissing_levels <- unique(col[!(col %in% c("N", "M"))])
    if (length(setdiff(nonmissing_levels, c("1", "2"))) == 0) {
      col[col == "1"] <- "0"
      col[col == "2"] <- "1"
    }
  } else if (inherits(col, c("numeric", "integer"))) {
    nonmissing_levels <- unique(col[!(col %in% c(-2, -1))])
    if (length(setdiff(nonmissing_levels, c(1, 2))) == 0) {
      col[col == 1] <- 0
      col[col == 2] <- 1
    }
  } else {
    warning("Encountered an unexpected variable type")
  }
  
  return(col)
}

# sum of column values / number of columns 
# if number of columns or MEMBER equal 0, return 0 
# if any column has -1 (missing data), output -1
# if any column has -2 or MEMBER column contains -2 or -1, output -2 
make_demo_perc <- function(df, cols, total_col) {
  out <- rowSums(df[, cols, drop=FALSE]) / df[[total_col]]
  out[df[[total_col]] == 0 | df$MEMBER == 0] <- 0
  out[apply(df[, cols, drop=FALSE], 1, function(x) any(x == -1))] <- -1
  out[apply(df[, cols, drop=FALSE], 1,
            function(x) any(x == -2)) | df$MEMBER %in% c(-2, -1)] <- -2
  
  return(out)
}

clean_ccd <- function(ccd) {
  rownames(ccd) <- NULL
  # keep Michigan schools only 
  ccd <- ccd[ccd$LSTATE == "MI",]
  
  # calculate percentage of gender
  gender_percs <- mapply(
    make_demo_perc,
    lapply(c("ALM", "ALF"), function(pf) paste0(RACEGENDER_PREFIXES, pf)),
    rep("MEMBER", 2),
    MoreArgs = list(df = ccd),
    SIMPLIFY = TRUE
  )
  colnames(gender_percs) <- paste(c("MALE", "FEMALE"), "PERC", sep = "_")
  
  # calculate percentage of G11 gender
  gender_g11_percs <- mapply(
    make_demo_perc,
    lapply(paste0(11, c("M", "F")),
           function(pf) paste0(RACEGENDER_PREFIXES, pf)),
    rep("G11", 2),
    MoreArgs = list(df = ccd),
    SIMPLIFY = TRUE
  )
  colnames(gender_g11_percs) <- paste(c("MALE", "FEMALE"), "G11_PERC", sep = "_")
  
  # calculate percentage of ethinicity 
  raceeth_percs <- mapply(
    make_demo_perc,
    RACEETH_COLS,
    rep("MEMBER", length(RACEETH_COLS)),
    MoreArgs = list(df = ccd),
    SIMPLIFY = TRUE
  )
  colnames(raceeth_percs) <- paste(RACEETH_COLS, "PERC", sep = "_")
  
  # calculate percentage of G11 ethinicity
  raceeth_g11_percs <- mapply(
    make_demo_perc,
    lapply(RACEGENDER_PREFIXES,
           function(pf) paste0(pf, 11, c("M", "F"))),
    rep("G11", length(RACEGENDER_PREFIXES)),
    MoreArgs = list(df = ccd),
    SIMPLIFY = TRUE
  )
  colnames(raceeth_g11_percs) <- paste(RACEETH_COLS, "G11_PERC", sep = "_")
  
  frl_perc <- make_demo_perc(ccd, FRL_COL, "MEMBER")
  
  cleaned_cat_cols <- sapply(ccd[, CCD_CAT_COLS], recode_binary_variable)
  
  # create a new column called merge_id by concatenating STID and SEASCH
  # extract columns listed in CCD_ID_COLS as a data frame 
  # calculate total enrollment for each row by summing up the values in CCD_ENROLLMENT_COLS
  # adds the rest of the columns we processed earlier 
  cbind(merge_id = paste0(ccd$STID, ccd$SEASCH),
        ccd[, CCD_ID_COLS, drop=FALSE],
        sapply(ccd[, CCD_ENROLLMENT_COLS, drop=FALSE], as.numeric),
        TOTAL_ENROLLMENT = as.numeric(apply(
          ccd[, CCD_ENROLLMENT_COLS], 1,
          function(r) sum(r[!(r %in% c(-2, -1))]))),
        cleaned_cat_cols,
        gender_percs, raceeth_percs, TOTFRL_PERC = frl_perc,
        gender_g11_percs, raceeth_g11_percs)
}

clean_scores <- function(all_scores, ccd) {
  all_scores <- all_scores[
    !is.na(all_scores$Building.Code) &
      !(all_scores$Building.Name %in% c("STATEWIDE", "ISDWIDE", "DISTRICTWIDE")),]

  stopifnot(sum(substr(colnames(all_scores),1,24)=="Average.Scale.Score.2011")==2)
  colnames(all_scores)[substr(colnames(all_scores),1,24
                              )=="Average.Scale.Score.2011"]  <-
      paste0("Average.Scale.Score.201", 2:1)
  
  out_cols <- Reduce(
    cbind,
    lapply(
      seq(2011, 2014),
      function(df, yr) {
        colnms <- paste0(c("Total.Tested.", "Average.Scale.Score."), yr)
        out_cols <- df[, colnms, drop=FALSE]
        out_cols[, 3] <- as.numeric(out_cols[, 1] == "<10")
        out_cols[, 4] <- as.numeric(out_cols[, 1] == "NULL")
        out_cols[, 5] <- as.numeric(out_cols[, 2] == "NULL" &
                                      !(out_cols[, 1] %in% c("<10", "NULL")))
        out_cols[out_cols[, 1] %in% c("<10", "NULL"), 1] <- NA_character_
        out_cols[out_cols[, 2] == "NULL", 2] <- NA_character_
        colnames(out_cols) <- c(colnms,
                                paste0(c("Masked.Average.Scale.Score.",
                                         "Missing.Students.Tested.",
                                         "Missing.Average.Scale.Score."), yr))
        sapply(out_cols, as.numeric)
      },
      df = all_scores
    )
  )
  all_scores[, grepl("(.Tested|Average.Scale)", colnames(all_scores))] <- NULL
  all_scores <- cbind(all_scores, out_cols)
  
  nchar_ccd_distid <- max(nchar(ccd$STID[!grepl("^[A-Z]", ccd$STID)]))
  nchar_ccd_schid <- max(nchar(ccd$SEASCH[!grepl("^[A-Z]", ccd$SEASCH)]))
  all_scores$District.Code <- zero_lpad(all_scores$District.Code, nchar_ccd_distid)
  all_scores$Building.Code <- zero_lpad(all_scores$Building.Code, nchar_ccd_schid)
  all_scores$merge_id <- paste0(all_scores$District.Code, all_scores$Building.Code)
  
  all_scores
}

##' Retrieve external data files from NCES and MI's CEPI
stopifnot(exists("extdataURLs"),
          require("readxl"), require("httr")
          )

.tf <- tempfile(fileext=".zip")
download.file(extdataURLs$CCD, .tf)
ccd <- read.delim(unz(.tf, "sc132a.txt"))
unlink(.tf)
cleaned_ccd <- clean_ccd(ccd)

.tf  <- tempfile(fileext=".zip")
headers = c(
  `user-agent` = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.5005.61 Safari/537.36'
)
httr::GET(url=extdataURLs$MME,
                 httr::add_headers(.headers=headers)
          )$content |>
writeBin(.tf)

.td  <- tempdir()
mme_xls  <-
    unzip(.tf,
          files="Spring2011-2014MMEFourYearDemographicDataFile-Sortable.xls",
          exdir=.td)

mme  <-  read_excel(mme_xls, .name_repair="universal")
rm(.tf, mme_xls)

all_schools <- clean_scores(mme, cleaned_ccd) |>
               merge(cleaned_ccd, by = "merge_id", all = TRUE)

analysis1data <- all_schools[
  !is.na(all_schools$DemographicGroup) &
    all_schools$DemographicGroup == "All Students" &
    all_schools$TYPE == "1" &
    all_schools$Subject == "M",]
rownames(analysis1data) <- NULL
analysis1data$schoolid <- analysis1data$merge_id
analysis1data$merge_id <- NULL
analysis2data <- all_schools[
  !is.na(all_schools$DemographicGroup) &
    all_schools$DemographicGroup %in% c("White", "Black or African American",
                                        "Hispanic of any race", "Asian",
                                        "American Indian or Alaska Native",
                                        "Two or more races") &
    all_schools$TYPE == "1" &
    all_schools$Subject == "M",]
analysis2data$schoolid <- analysis2data$merge_id
analysis2data$merge_id <- NULL
analysis2data$DemographicGroup[!(analysis2data$DemographicGroup %in%
                                   c("White", "Black or African American"))] <-
  "Other Race/Ethnicity"
rownames(analysis2data) <- NULL
analysis2data$uniqueid <- seq_len(nrow(analysis2data))
