test_that("Combining weighted specifications", {
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)

  w1 <- ate(spec, data = simdata[1:30,])
  w2 <- ate(spec, data = simdata[31:40,])
  w3 <- ate(spec, data = simdata[41:50,])

  c_w <- c(w1, w2, w3)
  expect_true(inherits(c_w, "numeric"))
  expect_length(c_w, 50)

  w1e <- ett(spec, data = simdata[1:30,])
  w2e <- ett(spec, data = simdata[31:40,])
  w3e <- ett(spec, data = simdata[41:50,])

  c_we <- c(w1e, w2e, w3e)
  expect_true(inherits(c_we, "numeric"))
  expect_length(c_we, 50)

  expect_message(c(w1, 1:5), "with a non-WeightedStudySpecification")
  expect_error(c(w1, w1e), "same target")

  spec2 <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  alt_w1 <- ate(spec2, data = simdata)

  expect_error(c(w1, alt_w1), "differing StudySpecification")

  # if the first argument is compatible with WeightedStudySpecification but isn't one (e.g.
  # numeric vector), c() will return a numeric vector
  #expect_true(inherits(c(1:5, w1), "WeightedStudySpecification"))
})

test_that("Combining weighted specifications with different dichotomys ", {
  spec <- rct_spec(dose ~ uoa(uoa1, uoa2), data = simdata)

  w1 <- ate(spec, data = simdata[1:10, ], dichotomy = dose >= 300 ~ .)
  w2 <- ate(spec, data = simdata[11:30, ], dichotomy = dose >= 200 ~ .)
  w3 <- ate(spec, data = simdata[31:50, ], dichotomy = dose >= 100 ~ .)

  c_w <- c(w1, w2, w3)
  expect_true(inherits(c_w, "numeric"))
  expect_length(c_w, 50)

  expect_warning(c(w1, w2, w3, warn_dichotomy_not_equal = TRUE),
                 "Concatenating")

  data_w  <- rbind(cbind(simdata[1:10, ], w=w1),
                   cbind(simdata[11:30, ], w=w2),
                   cbind(simdata[31:50, ], w=w3))
  expect_s4_class(data_w$w, "numeric")
  expect_true(inherits(data_w$w, "WeightedStudySpecification"))
})

test_that("Combine WeightedStudySpecifications & align weights with analysis data",{
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

    w20  <- ett(spec, data=subset(analysis_dat,year=="AY20"),
                dichotomy= year_trt<="AY20" ~ .)
    w20[which(subset(analysis_dat,year=="AY20")$id %in% c("a", "e"))] |>
    as.numeric() |> expect_equal(rep(0,2))
    ## ...confirming that blocks without treatment variation
    ## receive 0 weight. This is also tested by the test that
    ## "#130 zero weights with non-varying treatment in a block",
    ## in test.WeightedStudySpecification, but here we go on to check that
    ## in the next year, when there is variation in treatment status,
    ## the same pair has nonzero weights:
    w21  <- ett(spec, data=subset(analysis_dat,year=="AY21"),
               dichotomy= year_trt<="AY21" ~ .)
    w21[which(subset(analysis_dat,year=="AY21")$id %in% c("a", "e"))] |>
    as.numeric() |> expect_equal(rep(1,2))

    w0 <- c(w20, w21)
    expect_length(w0, length(w20)+length(w21))
    expect_true(inherits(w0, "numeric"))
    expect_false(inherits(w0, "WeightedStudySpecification"))

### Bringing the weights back into the data frame is easier
### if you're happy to reorder the data.
    mf_dat <- cbind(analysis_dat[order(analysis_dat$year),], w0)
    expect_identical(mf_dat$w0, w0)

### Users might try the following to put w0 into the same order as the
### rows of analysis_dat.
    analysis_dat$w0 <- numeric(12)
    analysis_dat[analysis_dat$year=="AY20","w0"] <- w20
    analysis_dat[analysis_dat$year=="AY21","w0"] <- w21
    expect_true(inherits(analysis_dat$w0, "numeric"))
    expect_false(inherits(analysis_dat$w0, "WeightedStudySpecification"))

    analysis_dat$w0 <- w0
    expect_true(inherits(analysis_dat$w0, "numeric"))
    analysis_dat[which(analysis_dat$year=="AY20"), "w0"]  <- w20
    analysis_dat[which(analysis_dat$year=="AY21"), "w0"]  <- w21
    expect_true(inherits(analysis_dat$w0, "numeric"))

### Alternatively, use lapply and unsplit:
    analysis_dat$w1 <-
        lapply(c("AY20", "AY21"),
    {\(yr) ett(spec, data=subset(analysis_dat,year==yr),
               dichotomy= as.formula(paste0("year_trt<=\"", yr, "\"~."))) }
    ) |> unsplit(analysis_dat$year)
    expect_s4_class(analysis_dat$w1, "WeightedStudySpecification")
    expect_equal(analysis_dat$w1@.Data, analysis_dat$w0@.Data)

### (One might hope to use tapply instead of lapply, but I couldn't find
###    a `dichotomy=` formula that's compatible with that mechanism.)
### Doesn't work:
###    tapply(analysis_dat, ~year,
###    {\(dat) ett(spec, data=dat, dichotomy = year_trt <= dat[1L, "year"] ~ .)}
###    )

})
