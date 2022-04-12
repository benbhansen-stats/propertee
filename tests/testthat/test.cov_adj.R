old_opt <- options()
on.exit(options(old_opt))
options(flexida_warn_on_conditional_treatment = FALSE)

test_that("cov_adj basics", {
  data(STARdata)
  STARdata$id <- seq_len(nrow(STARdata))
  des <- rct_design(stark == "small" ~ unitid(id), data = STARdata)
  cmod <- lm(readk ~ gender + ethnicity, data = STARdata)


  ca <- cov_adj(cmod, design = des)
  expect_true(is(ca, "numeric"))
  expect_true(is(ca, "vector"))
  expect_equal(length(ca), nrow(cmod$model))

  expect_error(cov_adj(cmod), "Unable to locate")
})

test_that("cov_adj as offset", {
  data(simdata)
  des <- obs_design(z ~ uoa(cid1, cid2), data = simdata)
  cmod <- lm(y ~ x, data = simdata)
  m <- lm(y ~ z, data = simdata,
          offset = cov_adj(cmod),
          weights = ate(des))
  expect_true(is(m$offset, "numeric"))
  expect_true(is(m$offset, "vector"))


  data(STARdata)
  STARdata$id <- seq_len(nrow(STARdata))
  des <- rct_design(stark == "small" ~ unitid(id), data = STARdata)
  cmod <- lm(readk ~ gender + ethnicity, data = STARdata)


  # This currently errors with
  # Error in model.frame.default(formula = readk ~ birth + lunchk, data =
  #STARdata, : variable lengths differ (found for '(weights)')
  #mod <- lm(readk ~ birth + lunchk, data = STARdata,
  #          offset = cov_adj(cmod, newdata = STARdata), weights = ate(des))

})





## ## Stop tidyverse from spamming the test output display
## options(tidyverse.quiet = TRUE)
## options(conflicts.policy = list(warn = FALSE))
## library(tidyverse)

## ## Loads the STAR data adn does preprocessing on it
## ## Returns only outcomes for kindergarten
## load_data <- function() {
##   data(STARdata)
##   STARdata$treatment <- STARdata$stark == "small"
##   STARdata$treatment[is.na(STARdata$treatment)] <- FALSE
##   STARdata$studentid <- as.character(1:nrow(STARdata))
##   STAR_pre <- STARdata[, c("studentid", "treatment",
##                      "gender", "ethnicity", "birth", "lunchk",  # individual demographics
##                      "schoolk", "degreek", "ladderk", "experiencek", "tethnicityk", # school and teacher demographics
##                      "systemk", "schoolidk" # school ID information
##                      )]

##   STAR_school <- group_by(STAR_pre, schoolidk) %>%
##     summarize(school_n = n(),
##               school_n1 = sum(treatment),
##               school_n0 = school_n - school_n1)

##   STAR_pre <- left_join(STAR_pre, STAR_school, by = "schoolidk") %>%
##    mutate(E_Z         = school_n1 / school_n,
##           weight_ate  = treatment / E_Z + (1 - treatment) / (1 - E_Z),
##           weight_ett  = treatment + (1 - treatment) * E_Z / (1 - E_Z),
##           weight_etc  = treatment * (1 - E_Z) / E_Z + (1 - treatment))

##   STAR_post <- rbind(
##     data.frame(studentid = STARdata$studentid, year = "k", read = STARdata$readk, math = STARdata$mathk, strings.as.factors = FALSE),
##     data.frame(studentid = STARdata$studentid, year = "1", read = STARdata$read1, math = STARdata$math1, strings.as.factors = FALSE),
##     data.frame(studentid = STARdata$studentid, year = "2", read = STARdata$read2, math = STARdata$math2, strings.as.factors = FALSE),
##     data.frame(studentid = STARdata$studentid, year = "3", read = STARdata$read3, math = STARdata$math3, strings.as.factors = FALSE))

##   STAR_pre_post <- inner_join(STAR_pre, STAR_post, by = "studentid")

##   rhs <- ~ gender + ethnicity + birth + lunchk +
##     ladderk + experiencek + tethnicityk + year

##   covariance_y0_read <- lm(update(rhs, read ~ .), data = STAR_pre_post, subset = !treatment)
##   covariance_y0_math <- lm(update(rhs, math ~ .), data = STAR_pre_post, subset = !treatment)

##   ## of !! is to turn numeric into logical
##   covariance_y1_read <- lm(update(rhs, read ~ .), data = STAR_pre_post, subset = !!treatment)
##   covariance_y1_math <- lm(update(rhs, math ~ .), data = STAR_pre_post, subset = !!treatment)

##   pred <- function(mod) {
##     predict(mod, newdata = STAR_pre_post, type = "response")
##   }

##   STAR_pre_post <-
##     mutate(STAR_pre_post,
##          read_y0_hat = pred(covariance_y0_read),
##          read_y1_hat = pred(covariance_y1_read),
##          math_y0_hat = pred(covariance_y0_math),
##          math_y1_hat = pred(covariance_y1_math))

##   return(list(STAR_pre = STAR_pre,
##               STAR_post = STAR_post,
##               STAR_pre_post = STAR_pre_post))
## }


## dw <- function(w, z, y, y0_hat) {
##   tmp <- data.frame(w, z, y, y0_hat) %>% na.omit
##   with(tmp,
##     sum(w * z * (y - y0_hat)) / sum(w * z) -
##       sum(w * (1 - z) * (y - y0_hat)) / sum(w * (1 - z))
##   )
## }


## test_that("Effect point estimates", {

##   ## load the TN STAR data set, with our modifications
##   STAR_data       <- load_data()
##   STAR_pre        <- STAR_data$STAR_pre
##   STAR_post       <- STAR_data$STAR_post
##   STAR_pre_post   <- STAR_data$STAR_pre_post
##   STAR_pre_post_k <- filter(STAR_pre_post, year == "k")

##   ## Get the ETT and ATE estimates using both `lm` and our direct implementation
##   ett_read_lm <- lm(read ~ treatment,
##                     data = STAR_pre_post_k,
##                     weights = weight_ett,
##                     offset = read_y0_hat)

##   ett_read_dw <- with(STAR_pre_post_k, dw(weight_ett, treatment, read, read_y0_hat))

##   ## coef will have names, so ignore attributes
##   expect_equal(coef(ett_read_lm)[2], ett_read_dw, ignore_attr = TRUE)

##   ate_read_lm <- lm(read ~ treatment,
##                   data = STAR_pre_post_k,
##                   weights = weight_ate,
##                   offset = read_y0_hat)

##   ate_read_dw <- with(STAR_pre_post_k, dw(weight_ate, treatment, read, read_y0_hat))
##   expect_equal(coef(ate_read_lm)[2], ate_read_dw, ignore_attr = TRUE)

##   ## now for conditional ATEs by ethnicity subgroups
##   ## to keep things simpler, we'll do the ATE estimates

##   ate_read_ethnicity_lm <- lm(read ~ treatment * ethnicity,
##                             data = STAR_pre_post_k,
##                             weights = weight_ate,
##                             offset = read_y0_hat)

##   ate_1 <- c(coef(ate_read_ethnicity_lm)[2] + c(cauc = 0, coef(ate_read_ethnicity_lm)[8:12]), NaN)

##   ate_2 <- STAR_pre_post_k %>%
##     group_by(ethnicity) %>%
##     summarize(ate = dw(weight_ate, treatment, read, read_y0_hat)) %>%
##     select(ate) %>% first

##   expect_equal(round(ate_1, 5), round(ate_2, 5), ignore_attr = TRUE)

##   ## so everything so far is checking our checks!
##   ## now actually test flexida routines

##   ## recreated from the set up code above
##   fmla <- read ~ gender + ethnicity + birth + lunchk +
##     ladderk + experiencek + tethnicityk + year

##   ## notice, this model is built on all years, not just k
##   covariance_y0_read <- lm(fmla, data = STAR_pre_post, subset = !treatment)

##   ## NB: treatment must be numeric right now, hence `1 * treatment`
##   ## NB: including a strata() arg causes an error
##   STAR_design <- rct_design(I(1 * treatment) ~ cluster(studentid), data = STAR_pre)
##   STAR_ate    <- ate(STAR_design)
##   STAR_ett    <- ett(STAR_design)

##   ett_read_flex <- lm(read ~ treatment,
##                       offset = cov_adj(covariance_y0_read),
##                       data = STAR_pre_post_k,
##                       weights = STAR_ett)
##   expect_equal(coef(ate_read_lm)[2], ate_read_dw, ignore_attr = TRUE)

##   ## notice the formula does not mention the blocking factor
##   ate_read_flex <- lm(read ~ treatment,
##                       offset = cov_adj(covariance_y0_read),
##                       data = STAR_pre_post_k,
##                       weights = STAR_ate)

##   ate_read_eth_flex <- lm(read ~ treatment*ethnicity,
##                           offset = cov_adj(covariance_y0_read),
##                           data = STAR_pre_post_k,
##                           weights = STAR_ate)

##   ate_3 <- c(coef(ate_read_eth_flex)[2] + c(cauc = 0, coef(ate_read_ethnicity_lm)[8:12]), NaN)
##   expect_equal(ate_3, ate_1, ignore_attr = TRUE)
## })

test_that("cov_adj variance estimates for orthogonal predictors", {
  library(sandwich)
  set.seed(230274373)
  k <- 25
  n <- 4 * k
  x1 <- c(rep(1, 2 * k), rep(-1, 2 * k))
  x2 <- c(rep(1, k), rep(-1, 2 * k), rep(1, k))
  z <- c(rep(1, k), rep(-1, k), rep(1, k), rep(-1, k))

  # all orthogonal by design
  expect_equal((t(x1) %*% x2)[1], 0)
  expect_equal((t(x1) %*% z)[1], 0)
  expect_equal((t(x2) %*% z)[1], 0)

  xy <- function(b0 =  1, b1 = 2, b2 = 3, b3 = -1, sigma2 = 4) {
    y <- b0 + b1 * x1 + b2 * x2 + b3 * z + rnorm(n, sd = sqrt(sigma2))
    data.frame(id = 1:n, x1 = x1, x2 = x2, z = z, y = y)
  }

  df <- xy()
  mboth <- glm(y ~ x1 + x2 + z, data = df)
  m1 <- glm(y ~ x1, data = df)

  # we'll start by manually doing the offest Y - \hat Y
  m2 <- glm(y ~ x2 + z, data = df, offset = predict(m1))

  # Because of orthogonality, we should get $\hat \beta_0= \hat \alpha_0$,
  # $\hat \beta_1 = \hat \alpha_1$, and $\hat beta_2 = \hat gamma_1$.

  expect_equal(coef(mboth), c(coef(m1), coef(m2)[2:3]))

  # now repeat with cov_adj
  design <- rct_design(z == 1 ~ unitid(id), data = df)
  m2ca <- glm(y ~ x2 + z, data = df, offset = cov_adj(m1, design = design))

  expect_equal(coef(mboth), c(coef(m1), coef(m2ca)[2:3]))

  ## naive case
  m2naive <- glm(y ~ x2 + z, data = df)

  expect_false(all(vcov(m2naive) == vcov(m2ca)))

  hc_both  <- vcovHC(mboth, type = "HC0")
  hc_naive <- vcovHC(m2naive, type = "HC0")
  hc_m2ca  <- vcovHC(m2ca, type = "HC0")

  # trim down the bigger model to just the variables in the second stage
  in_two <- c("(Intercept)", "x2", "z")
  hc_both_trimmed <- hc_both[in_two, in_two]

  expect_equal(hc_m2ca, hc_both_trimmed)
  expect_false(all(hc_naive == hc_m2ca))


})


test_that("cov_adj variance estimates for correlated predictors", {
  library(sandwich)
  set.seed(230274373)
  k <- 25
  n <- 4 * k
  x1 <- c(rep(1, 2 * k), rep(-1, 2 * k))
  z <- sample(c(rep(1, 2 * k), rep(0, 2 * k)), prob = x1 + 2)

  xy <- function(b0 =  1, b1 = 2, b2 = 3, sigma2 = 4) {
    y <- b0 + b1 * x1 + b2 * z + rnorm(n, sd = sqrt(sigma2))
    data.frame(id = 1:n, x1 = x1,  z = z, y = y)
  }

  df <- xy()
  mboth <- glm(y ~ x1 + z, data = df)
  m1 <- glm(y ~ x1, data = df)

  # we'll start by manually doing the offest Y - \hat Y
  m2 <- glm(y ~ z, data = df, offset = predict(m1))

  # now repeat with cov_adj
  design <- rct_design(z == 1 ~ unitid(id), data = df)
  m2ca <- glm(y ~  z, data = df, offset = cov_adj(m1, design = design))
  ## naive case
  m2naive <- glm(y ~ z, data = df)

  expect_false(all(vcov(m2naive) == vcov(m2ca)))

  hc_both  <- vcovHC(mboth, type = "HC0")
  hc_naive <- vcovHC(m2naive, type = "HC0")
  hc_m2ca  <- vcovHC(m2ca, type = "HC0")

  # trim down the bigger model to just the variables in the second stage
  in_two <- c("(Intercept)", "z")
  hc_both_trimmed <- hc_both[in_two, in_two]

#  expect_equal(hc_m2ca, hc_both_trimmed)
#  expect_false(all(hc_naive == hc_m2ca))


})

options(old_opt)
