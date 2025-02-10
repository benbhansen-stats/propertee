test_that("dichotomy args where specification has no blocks, no time-varying assignment", {
  data(simdata)
  spec1 <- rct_spec(dose ~ cluster(uoa1, uoa2), simdata)

  lm1 <- lm(y ~ assigned(spec1, dichotomy = dose >250~dose<=250), simdata)
  expect_equal(lm1$model[[2]], as.numeric(simdata$dose>250))

  lm2 <- lm(y ~ assigned(spec1, dichotomy = dose >250~dose<=250), simdata,
            weights = ate(spec1, dichotomy = dose >250~dose<=250))
  lm3 <- lm(y ~ assigned(spec1), simdata,
            weights = ate(spec1, dichotomy = dose >250~dose<=250))
  wts <- ate(spec1, data = simdata, dichotomy = dose >250~dose<=250)
  lm4 <- lm(y ~ assigned(spec1, dichotomy = dose >250~dose<=250), simdata,
            weights = wts)

  expect_true(all.equal(lm2$coefficients, lm3$coefficients, check.attributes = FALSE))
  expect_true(all.equal(lm2$coefficients, lm4$coefficients, check.attributes = FALSE))
  expect_equal(lm2$model[[2]], lm1$model[[2]])
  expect_equal(lm3$model[[2]], lm1$model[[2]])
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(lm2$model$`(weights)`, slot)),
                           deparse1(methods::slot(wts, slot)))
               } else {
                 identical(methods::slot(lm2$model$`(weights)`, slot),
                           methods::slot(wts, slot))
               },
               logical(1L)))
  )
  expected_wts <- numeric(nrow(simdata))
  expected_wts[simdata$dose >250] <- 1 / mean(simdata$dose>250)
  expected_wts[simdata$dose <=250] <- 1 / mean(simdata$dose<=250)
  expect_equal(lm2$model$`(weights)`@.Data, expected_wts)

  lmitt1 <- lmitt(y ~ 1, specification = spec1, simdata, dichotomy = dose>250~dose<=250)
  expect_equal(lmitt1$model$dose., as.numeric(simdata$dose>250))
  expect_true(all.equal(lmitt1$coefficients[!grepl("^y:", names(lmitt1$coefficients))],
                        lm1$coefficients, check.attributes = FALSE))

  lmitt2 <- lmitt(y ~ 1, specification = spec1, simdata,
                  weights = ate(spec1, dichotomy = dose >250~dose<=250))
  lmitt3 <- lmitt(y ~ 1, specification = spec1, simdata, dichotomy = dose >250~dose<=250,
                  weights = ate(spec1, dichotomy = dose >250~dose<=250))
  lmitt4 <- lmitt(y ~ 1, specification = spec1, simdata, dichotomy = dose >250~dose<=250,
                  weights = ate(spec1))
  lmitt5 <- lmitt(y ~ 1, specification = spec1, simdata, dichotomy = dose >250~dose<=250,
                  weights = wts)

  expect_equal(lmitt2$coefficients, lmitt3$coefficients)
  expect_equal(lmitt2$coefficients, lmitt4$coefficients)
  expect_equal(lmitt2$coefficients, lmitt5$coefficients)
  expect_true(all.equal(lmitt2$coefficients[!grepl("^y:", names(lmitt2$coefficients))],
                        lm2$coefficients, check.attributes = FALSE))
  expect_equal(lmitt2$model[[2]], lmitt1$model[[2]])
  expect_equal(lmitt3$model[[2]], lm1$model[[2]])
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(lmitt2$model$`(weights)`, slot)),
                           deparse1(methods::slot(wts, slot)))
               } else {
                 identical(methods::slot(lmitt2$model$`(weights)`, slot),
                           methods::slot(wts, slot))
               },
               logical(1L)))
  )
  expect_equal(lmitt2$model$`(weights)`@.Data, expected_wts)

})

test_that("dichotomy args where specification has blocks, no time-varying assignment", {
  data(simdata)
  spec1 <- rct_spec(dose ~ cluster(uoa1, uoa2) + block(bid), simdata)

  lm1 <- lm(y ~ assigned(spec1, dichotomy = dose >250~dose<=250), simdata)
  expect_equal(lm1$model[[2]], as.numeric(simdata$dose>250))

  lm2 <- lm(y ~ assigned(spec1, dichotomy = dose >250~dose<=250), simdata,
            weights = ate(spec1, dichotomy = dose >250~dose<=250))
  lm3 <- lm(y ~ assigned(spec1), simdata,
            weights = ate(spec1, dichotomy = dose >250~dose<=250))
  wts <- ate(spec1, data = simdata, dichotomy = dose >250~dose<=250)
  lm4 <- lm(y ~ assigned(spec1, dichotomy = dose >250~dose<=250), simdata,
            weights = wts)

  expect_true(all.equal(lm2$coefficients, lm3$coefficients, check.attributes = FALSE))
  expect_true(all.equal(lm2$coefficients, lm4$coefficients, check.attributes = FALSE))
  expect_equal(lm2$model[[2]], lm1$model[[2]])
  expect_equal(lm3$model[[2]], lm1$model[[2]])
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(lm2$model$`(weights)`, slot)),
                           deparse1(methods::slot(wts, slot)))
               } else {
                 identical(methods::slot(lm2$model$`(weights)`, slot),
                           methods::slot(wts, slot))
               },
               logical(1L)))
  )
  expected_wts <- numeric(nrow(simdata))
  # in block 1, no uoa's have dose > 250, while in blocks 2 and 1/3 do
  expected_wts[simdata$bid == 1] <- 0 #
  expected_wts[simdata$dose >250 & simdata$bid %in% c(2, 3)] <- 3
  expected_wts[simdata$dose <=250 & simdata$bid %in% c(2, 3)] <- 3/2
  expect_equal(lm2$model$`(weights)`@.Data, expected_wts)

  suppressMessages(
    lmitt1 <- lmitt(y ~ 1, specification = spec1, simdata, dichotomy = dose>250~dose<=250)
  )
  expect_equal(lmitt1$model$dose., as.numeric(simdata$dose>250))

  lmitt2 <- lmitt(y ~ 1, specification = spec1, simdata,
                  weights = ate(spec1, dichotomy = dose >250~dose<=250))
  lmitt3 <- lmitt(y ~ 1, specification = spec1, simdata, dichotomy = dose >250~dose<=250,
                  weights = ate(spec1, dichotomy = dose >250~dose<=250))
  lmitt4 <- lmitt(y ~ 1, specification = spec1, simdata, dichotomy = dose >250~dose<=250,
                  weights = ate(spec1))
  lmitt5 <- lmitt(y ~ 1, specification = spec1, simdata, dichotomy = dose >250~dose<=250,
                  weights = wts)

  expect_equal(lmitt2$coefficients, lmitt3$coefficients)
  expect_equal(lmitt2$coefficients, lmitt4$coefficients)
  expect_equal(lmitt2$coefficients, lmitt5$coefficients)
  expect_equal(lmitt2$model[[2]], lmitt1$model[[2]])
  expect_equal(lmitt3$model[[2]], lm1$model[[2]])
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(lmitt2$model$`(weights)`, slot)),
                           deparse1(methods::slot(wts, slot)))
               } else {
                 identical(methods::slot(lmitt2$model$`(weights)`, slot),
                           methods::slot(wts, slot))
               },
               logical(1L)))
  )
  expect_equal(lmitt2$model$`(weights)`@.Data, expected_wts)

})

test_that("dichotomy args where specification has blocks, time-varying assignment", {
  set.seed(202021)
  analysis_dat  <- data.frame(
    id = rep(letters[1:6], each = 2),
    year = ordered(rep(c("AY20", "AY21"), 6), levels=c("AY20", "AY21", ".")),
    y = rnorm(12)
  )
  specification_dat  <- data.frame(id=letters[1:6], blk=character(6),
                            year_trt = ordered(rep(c("AY21", "AY20", "."),
                                                   each = 2),
                                               levels=c("AY20", "AY21", ".")
                            )
  )
  specification_dat[   specification_dat$id %in% c("a", "e") , "blk"]  <- "A"
  specification_dat[ !(specification_dat$id %in% c("a", "e")), "blk"]  <- "B"
  spec  <- obs_spec(year_trt~uoa(id)+block(blk), data=specification_dat)

  wts <-
    lapply(c("AY20", "AY21"),
           {\(yr) ett(spec, data=subset(analysis_dat,year==yr),
                      dichotomy= as.formula(paste0("year_trt<=\"", yr, "\"~."))) }
    ) |> unsplit(analysis_dat$year)

  lm1 <- lm(y ~ assigned(spec, dichotomy = year_trt<=year~.), data = analysis_dat)
  analysis_dat_year_trts <- setNames(specification_dat$year_trt, specification_dat$id)[analysis_dat$id]
  expect_equal(lm1$model[[2]], as.numeric(analysis_dat_year_trts <= analysis_dat$year))

  lm2 <- lm(y ~ assigned(spec, dichotomy = year_trt<=year~.), data = analysis_dat,
            weights = wts)
  expected_wts <- numeric(nrow(analysis_dat))
  expected_wts[analysis_dat$id %in% c("a", "e") & analysis_dat$year == "AY20"] <- 0
  expected_wts[analysis_dat$id %in% c("a", "e") & analysis_dat$year == "AY21"] <- 1
  expected_wts[analysis_dat$id %in% c("c", "d") & analysis_dat$year == "AY20"] <- 1
  expected_wts[analysis_dat$id %in% c("b", "f") & analysis_dat$year == "AY20"] <- 1
  expected_wts[analysis_dat$id %in% c("b", "c", "d") & analysis_dat$year == "AY21"] <- 1
  expected_wts[analysis_dat$id %in% c("f") & analysis_dat$year == "AY21"] <- 3
  expect_equal(wts@.Data, expected_wts)

  suppressMessages(
    lmitt1 <- lmitt(y ~ 1, specification = spec, analysis_dat, dichotomy = year_trt<=year ~.)
  )
  expect_equal(lmitt1$model$year_trt., as.numeric(analysis_dat$year>=analysis_dat_year_trts))

  lmitt2 <- lmitt(y ~ 1, specification = spec, analysis_dat, dichotomy = year_trt<=year ~.,
                  weights = wts)
  expect_equal(lmitt2$model$year_trt., as.numeric(analysis_dat$year>=analysis_dat_year_trts))
})
