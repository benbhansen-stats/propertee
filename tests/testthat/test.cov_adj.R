## Stop tidyverse from spamming the test output display
options(tidyverse.quiet = TRUE) 
options(conflicts.policy = list(warn = FALSE))
library(tidyverse)

## Loads the STAR data adn does preprocessing on it
## Returns only outcomes for kindergarten
load_data <- function() {
  library(AER, quietly = TRUE)
  data(STAR)
  STAR$treatment <- STAR$stark == "small"
  STAR$treatment[is.na(STAR$treatment)] <- FALSE
  STAR$studentid <- as.character(1:nrow(STAR))
  STAR_pre <- STAR[, c("studentid", "treatment",
                     "gender", "ethnicity", "birth", "lunchk",  # individual demographics
                     "schoolk", "degreek", "ladderk", "experiencek", "tethnicityk", # school and teacher demographics
                     "systemk", "schoolidk" # school ID information
                     )]
  
  STAR_school <- group_by(STAR_pre, schoolidk) %>% 
    summarize(school_n = n(), 
              school_n1 = sum(treatment), 
              school_n0 = school_n - school_n1)
  
  STAR_pre <- left_join(STAR_pre, STAR_school, by = "schoolidk") %>%
   mutate(E_Z         = school_n1 / school_n,
          weight_ate  = treatment / E_Z + (1 - treatment) / (1 - E_Z),
          weight_ett  = 1 + (1 - treatment) * E_Z / (1 - E_Z),
          weight_etc  = treatment * (1 - E_Z) / E_Z + (1 - treatment))
  
  STAR_post <- rbind(
    data.frame(studentid = STAR$studentid, year = "k", read = STAR$readk, math = STAR$mathk, strings.as.factors = FALSE),
    data.frame(studentid = STAR$studentid, year = "1", read = STAR$read1, math = STAR$math1, strings.as.factors = FALSE),
    data.frame(studentid = STAR$studentid, year = "2", read = STAR$read2, math = STAR$math2, strings.as.factors = FALSE),
    data.frame(studentid = STAR$studentid, year = "3", read = STAR$read3, math = STAR$math3, strings.as.factors = FALSE))

  STAR_pre_post <- inner_join(STAR_pre, STAR_post, by = "studentid")
  
  rhs <- ~ gender + ethnicity + birth + lunchk + 
    ladderk + experiencek + tethnicityk + year

  covariance_y0_read <- lm(update(rhs, read ~ .), data = STAR_pre_post, subset = !treatment)
  covariance_y0_math <- lm(update(rhs, math ~ .), data = STAR_pre_post, subset = !treatment)

  ## of !! is to turn numeric into logical
  covariance_y1_read <- lm(update(rhs, read ~ .), data = STAR_pre_post, subset = !!treatment)
  covariance_y1_math <- lm(update(rhs, math ~ .), data = STAR_pre_post, subset = !!treatment)
  
  pred <- function(mod) {
    predict(mod, newdata = STAR_pre_post, type = "response")
  }

  STAR_pre_post <- 
    mutate(STAR_pre_post,
         read_y0_hat = pred(covariance_y0_read),
         read_y1_hat = pred(covariance_y1_read),
         math_y0_hat = pred(covariance_y0_math),
         math_y1_hat = pred(covariance_y1_math))
  
  return(list(STAR_pre = STAR_pre,
              STAR_post = STAR_post,
              STAR_pre_post = STAR_pre_post))
}


dw <- function(w, z, y, y0_hat) {
  tmp <- data.frame(w, z, y, y0_hat) %>% na.omit
  with(tmp, 
    sum(w * z * (y - y0_hat)) / sum(w * z) - 
      sum(w * (1 - z) * (y - y0_hat)) / sum(w * (1 - z))
  )
}


test_that("Effect point estimates", { 
  
  ## load the TN STAR data set, with our modifications
  STAR_data       <- load_data()
  STAR_pre        <- STAR_data$STAR_pre
  STAR_post       <- STAR_data$STAR_post
  STAR_pre_post   <- STAR_data$STAR_pre_post
  STAR_pre_post_k <- filter(STAR_pre_post, year == "k") 
  
  ## Get the ETT and ATE estimates using both `lm` and our direct implementation
  ett_read_lm <- lm(read ~ treatment, 
                    data = STAR_pre_post_k,
                    weights = weight_ett,
                    offset = read_y0_hat)
  
  ett_read_dw <- with(STAR_pre_post_k, dw(weight_ett, treatment, read, read_y0_hat))
  
  ## coef will have names, so ignore attributes
  expect_equal(coef(ett_read_lm)[2], ett_read_dw, ignore_attr = TRUE)
  
  ate_read_lm <- lm(read ~ treatment, 
                  data = STAR_pre_post_k, 
                  weights = weight_ate,
                  offset = read_y0_hat)

  ate_read_dw <- with(STAR_pre_post_k, dw(weight_ate, treatment, read, read_y0_hat))
  expect_equal(coef(ate_read_lm)[2], ate_read_dw, ignore_attr = TRUE)
  
  ## now for conditional ATEs by ethnicity subgroups
  ## to keep things simpler, we'll do the ATE estimates
  
  ate_read_ethnicity_lm <- lm(read ~ treatment * ethnicity, 
                            data = STAR_pre_post_k,
                            weights = weight_ate,
                            offset = read_y0_hat)
  
  ate_1 <- c(coef(ate_read_ethnicity_lm)[2] + c(cauc = 0, coef(ate_read_ethnicity_lm)[8:12]), NaN)
  
  ate_2 <- STAR_pre_post_k %>% 
    group_by(ethnicity) %>%
    summarize(ate = dw(weight_ate, treatment, read, read_y0_hat)) %>%
    select(ate) %>% first
  
  expect_equal(round(ate_1, 5), round(ate_2, 5), ignore_attr = TRUE)
  
  ## so everything so far is checking our checks!
  ## now actually test flexida routines
  
  ## recreated from the set up code above
  fmla <- read ~ gender + ethnicity + birth + lunchk + 
    ladderk + experiencek + tethnicityk + year

  ## notice, this model is built on all years, not just k 
  covariance_y0_read <- lm(fmla, data = STAR_pre_post, subset = !treatment)
  
  ## NB: treatment must be numeric right now, hence `1 * treatment`
  ## NB: including a strata() arg causes an error
  STAR_design <- RCT_Design(I(1 * treatment) ~ cluster(studentid), data = STAR_pre) 
  STAR_ate    <- ate(STAR_design)
  STAR_ett    <- ett(STAR_design)
  
  ett_read_flex <- lm(read ~ treatment, 
                      data = STAR_pre_post_k,
                      weights = cov_adj(STAR_ett, covariance_y0_read))
  expect_equal(coef(ate_read_lm)[2], ate_read_dw, ignore_attr = TRUE)
  
  ## notice the formula does not mention the blocking factor and no need for the offset
  ate_read_flex <- lm(read ~ treatment, 
                      data = STAR_pre_post_k,
                      weights = cov_adj(STAR_ate, covariance_y0_read))
  
  ate_read_eth_flex <- lm(read ~ treatment*ethnicity,
                          data = STAR_pre_post_k,
                          weights = cov_adj(STAR_ate, covariance_y0_read))
  
  ate_3 <- c(coef(ate_read_eth_flex)[2] + c(cauc = 0, coef(ate_read_ethnicity_lm)[8:12]), NaN)
  expect_equal(ate_3, ate_1, ignore_attr = TRUE)
})