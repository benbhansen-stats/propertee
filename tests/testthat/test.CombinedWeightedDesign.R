test_that("Combining weighted designs", {
  data(simdata)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)

  w1 <- ate(des, data = simdata[1:30,])
  w2 <- ate(des, data = simdata[31:40,])
  w3 <- ate(des, data = simdata[41:50,])

  c_w <- c(w1, w2, w3)
  expect_true(inherits(c_w, "numeric"))
  expect_length(c_w, 50)

  w1e <- ett(des, data = simdata[1:30,])
  w2e <- ett(des, data = simdata[31:40,])
  w3e <- ett(des, data = simdata[41:50,])

  c_we <- c(w1e, w2e, w3e)
  expect_true(inherits(c_we, "numeric"))
  expect_length(c_we, 50)

  expect_message(c(w1, 1:5), "with a non-WeightedDesign")
  expect_error(c(w1, w1e), "same target")

  des2 <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  alt_w1 <- ate(des2, data = simdata)

  expect_error(c(w1, alt_w1), "differing Design")

  # if the first argument is compatible with WeightedDesign but isn't one (e.g.
  # numeric vector), c() will return a numeric vector
  #expect_true(inherits(c(1:5, w1), "WeightedDesign"))
})

test_that("Combining weighted designs with different dichotomys ", {
  des <- rct_design(dose ~ uoa(uoa1, uoa2), data = simdata)

  w1 <- ate(des, data = simdata[1:10, ], dichotomy = dose >= 300 ~ .)
  w2 <- ate(des, data = simdata[11:30, ], dichotomy = dose >= 200 ~ .)
  w3 <- ate(des, data = simdata[31:50, ], dichotomy = dose >= 100 ~ .)

  c_w <- c(w1, w2, w3)
  expect_true(inherits(c_w, "numeric"))
  expect_length(c_w, 50)
  
  expect_warning(c(w1, w2, w3, warn_dichotomy_not_equal = TRUE),
                 "Concatenating")

  data_w  <- rbind(cbind(simdata[1:10, ], w=w1),
                   cbind(simdata[11:30, ], w=w2),
                   cbind(simdata[31:50, ], w=w3))
  expect_s4_class(data_w$w, "numeric")
  ##would prefer the following were true:
  expect_true(inherits(data_w$w, "WeightedDesign"))
  expect_false(inherits(data_w$w, "CombinedWeightedDesign"))
  ##ToDo: WD's combined using `rbind()` should give CWD's
})

test_that("Combine WeightedDesigns & align weights with analysis data",{
    set.seed(202021)
    analysis_dat  <- data.frame(
        id = rep(letters[1:6], each = 2),
        year = ordered(rep(c("AY20", "AY21"), 6), levels=c("AY20", "AY21", ".")),
        y = rnorm(12)
    )
    design_dat  <- data.frame(id=letters[1:6], blk=character(6),
                              year_trt = ordered(rep(c("AY21", "AY20", "."),
                                                     each = 2),
                                                 levels=c("AY20", "AY21", ".")
                                                 )
                              )
    design_dat[   design_dat$id %in% c("a", "e") , "blk"]  <- "A"
    design_dat[ !(design_dat$id %in% c("a", "e")), "blk"]  <- "B"
    des  <- obs_design(year_trt~uoa(id)+block(blk), data=design_dat)

    w20  <- ett(des, data=subset(analysis_dat,year=="AY20"),
                dichotomy= year_trt<="AY20" ~ .)
    w20[which(subset(analysis_dat,year=="AY20")$id %in% c("a", "e"))] |>
    as.numeric() |> expect_equal(rep(0,2))
    ## ...confirming that blocks without treatment variation
    ## receive 0 weight. This is also tested by the test that
    ## "#130 zero weights with non-varying treatment in a block",
    ## in test.WeightedDesign, but here we go on to check that
    ## in the next year, when there is variation in treatment status,
    ## the same pair has nonzero weights:
    w21  <- ett(des, data=subset(analysis_dat,year=="AY21"),
               dichotomy= year_trt<="AY21" ~ .)
    w21[which(subset(analysis_dat,year=="AY21")$id %in% c("a", "e"))] |>
    as.numeric() |> expect_equal(rep(1,2))

    w0 <- c(w20, w21)
    expect_length(w0, length(w20)+length(w21))
    expect_true(inherits(w0, "numeric"))
    expect_false(inherits(w0, "WeightedDesign"))
    expect_false(inherits(w0, "CombinedWeightedDesign"))


### Bringing the weights back into the data frame is easier
### if you're happy to reorder the data.
    mf_dat <- cbind(analysis_dat[order(analysis_dat$year),], w0)
    expect_identical(mf_dat$w0, w0)

### Let obvious how to put w0 into the same order as the
### rows of analysis_dat.  Users might try the following.
    analysis_dat$w0 <- numeric(12)
    analysis_dat[analysis_dat$year=="AY20","w0"] <- w20
    analysis_dat[analysis_dat$year=="AY21","w0"] <- w21
    expect_true(inherits(analysis_dat$w0, "numeric"))
    expect_false(inherits(analysis_dat$w0, "WeightedDesign"))
    expect_false(inherits(analysis_dat$w0, "CombinedWeightedDesign"))
    ## rather, this workflow requires us to make sure a_d$w0 is a
    ## WeightedDesign &/or CombinedWeightedDesign from the get-go.
    ## One way to do this:
    analysis_dat$w0 <- w0
    expect_true(inherits(analysis_dat$w0, "numeric"))
    analysis_dat[which(analysis_dat$year=="AY20"), "w0"]  <- w20
    analysis_dat[which(analysis_dat$year=="AY21"), "w0"]  <- w21
    expect_true(inherits(analysis_dat$w0, "numeric"))

### Alternatively, use lapply and unsplit:
    analysis_dat$w1 <-
        lapply(c("AY20", "AY21"),
    {\(yr) ett(des, data=subset(analysis_dat,year==yr),
               dichotomy= as.formula(paste0("year_trt<=\"", yr, "\"~."))) }
    ) |> unsplit(analysis_dat$year)
    expect_s4_class(analysis_dat$w1, "WeightedDesign")
    expect_equal(analysis_dat$w1@.Data, analysis_dat$w0@.Data)
    
    ## ...except that the below should really be `TRUE`:
    expect_false(inherits(analysis_dat$w1, "CombinedWeightedDesign"))
    ## ToDo: Ensure that `unsplit()`, `split<-` etc create
    ## CWDs carrying appropriate info

### (One might hope to use tapply instead of lapply, but I couldn't find
###    a `dichotomy=` formula that's compatible with that mechanism.)
### Doesn't work:
###    tapply(analysis_dat, ~year,
###    {\(dat) ett(des, data=dat, dichotomy = year_trt <= dat[1L, "year"] ~ .)}
###    )

})
