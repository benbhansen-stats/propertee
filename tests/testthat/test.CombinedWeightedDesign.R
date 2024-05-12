test_that("Combining weighted designs", {
  data(simdata)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)

  w1 <- ate(des, data = simdata[1:30,])
  w2 <- ate(des, data = simdata[31:40,])
  w3 <- ate(des, data = simdata[41:50,])

  c_w <- c(w1, w2, w3)
  expect_true(inherits(c_w, "WeightedDesign"))
  expect_length(c_w, 50)
  expect_identical(c_w, ate(des, data = simdata))

  w1e <- ett(des, data = simdata[1:30,])
  w2e <- ett(des, data = simdata[31:40,])
  w3e <- ett(des, data = simdata[41:50,])

  c_we <- c(w1e, w2e, w3e)
  expect_true(inherits(c_we, "WeightedDesign"))
  expect_length(c_we, 50)
  expect_identical(c_we, ett(des, data = simdata))

  expect_error(c(w1, 1:5), "with other")
  expect_error(c(w1, w1e), "same target")

  des2 <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  alt_w1 <- ate(des2, data = simdata)

  expect_error(c(w1, alt_w1), "which differ on elements")

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
  expect_true(inherits(c_w, "WeightedDesign"))
  expect_length(c_w, 50)

  expect_error(c(w1, w2, w3, force_dichotomy_equal = TRUE),
               "must be identical")

  data_w  <- rbind(cbind(simdata[1:10, ], w=w1),
                   cbind(simdata[11:30, ], w=w2),
                   cbind(simdata[31:50, ], w=w3))
  expect_s4_class(data_w$w, "WeightedDesign")
  ##would prefer the following were true:
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
                              year_trt = ordered(rep(c("AY20", "AY21", "."),
                                                     each = 2),
                                                 levels=c("AY20", "AY21", ".")
                                                 )
                              )
    design_dat[   design_dat$id %in% c("a", "e") , "blk"]  <- "A"
    design_dat[ !(design_dat$id %in% c("a", "e")), "blk"]  <- "B"
    des  <- obs_design(year_trt~uoa(id)+block(blk), data=design_dat)

    w20  <- ett(des, data=subset(analysis_dat,year=="AY20"),
               dichotomy= year_trt<="AY20" ~ .)
    w21  <- ett(des, data=subset(analysis_dat,year=="AY21"),
               dichotomy= year_trt<="AY21" ~ .)
    w0 <- c(w20, w21)
    expect_length(w0, length(w20)+length(w21))
    expect_true(inherits(w0, "WeightedDesign"))
    expect_true(inherits(w0, "CombinedWeightedDesign"))


### Bringing the weights back into the data frame is easier
### if you're happy to reorder the data.
    mf_dat <- cbind(analysis_dat[order(analysis_dat$year),], w0)
    expect_true(inherits(mf_dat$w0, "CombinedWeightedDesign"))
    expect_equal(mf_dat$w0, w0)

### Let obvious how to put w0 into the same order as the
### rows of analysis_dat.  Users might try the following.
    analysis_dat$w0 <- numeric(12)
    analysis_dat[analysis_dat$year=="AY20","w0"] <- w20
    analysis_dat[analysis_dat$year=="AY21","w0"] <- w21
    ## would have been preferable that the below 2 assertions be true:
    expect_false(inherits(analysis_dat$w, "WeightedDesign"))
    expect_false(inherits(analysis_dat$w, "CombinedWeightedDesign"))
    ## rather, this workflow requires us to make sure a_d$w0 is a
    ## WeightedDesign &/or CombinedWeightedDesign from the get-go.
    ## One way to do this:
    analysis_dat$w0 <- w0
    expect_s4_class(analysis_dat$w0, "CombinedWeightedDesign")
    analysis_dat[which(analysis_dat$year=="AY20"), "w0"]  <- w20
    analysis_dat[which(analysis_dat$year=="AY21"), "w0"]  <- w21
    expect_s4_class(analysis_dat$w0, "CombinedWeightedDesign")
    expect_identical(w0@Design, analysis_dat$w0@Design)
    expect_equal(as.numeric(w20), as.numeric(analysis_dat$w0)[analysis_dat$year=="AY20"])
    expect_equal(as.numeric(w21), as.numeric(analysis_dat$w0)[analysis_dat$year=="AY21"])
    ## Now let's confirm that the CombinedWeightedDesign internals
    ## are as they should be:
    expect_equal(mf_dat$w0@dichotomies, analysis_dat$w0@dichotomies)
    ## (`expect_setequal()` would work too, but `expect_equal()`
    ## sets up the next test).
    ## This ought to have been true:
    expect_false(all(analysis_dat[analysis_dat$w0@keys[[1]], "year"]=="AY20"))
    ## ... just as it is when we've reordered the data to match
    ##the CWD, rather than the reverse.
    expect_true(all(mf_dat[mf_dat$w0@keys[[1]], "year"]=="AY20"))
    ## ToDo: Ensure CWD@keys get reordered upon reorder of the .Data

### Alternatively, use lapply and unsplit:
    analysis_dat$w1 <-
        lapply(c("AY20", "AY21"),
    {\(yr) ett(des, data=subset(analysis_dat,year==yr),
               dichotomy= year_trt<=yr ~ .) }
    ) |> unsplit(analysis_dat$year)
    expect_equal(as.numeric(analysis_dat$w1),
                 as.numeric(analysis_dat$w0))
    expect_identical(analysis_dat$w1@Design@structure,
                     analysis_dat$w0@Design@structure)
    expect_s4_class(analysis_dat$w1, "WeightedDesign")
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
