test_that("Treatment inputs", {
  data(simdata)

  # numeric 0/1 treatment
  d <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.numeric(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))
  expect_silent(mod <- lmitt(y ~ 1, data = simdata, design = d))

  # numeric non-0/1
  d <- rct_design(o ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.numeric(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))
  expect_silent(mod <- lmitt(y ~ 1, data = simdata, design = d))

  # factor 0/1
  simdata$foo <- as.factor(simdata$z)
  d <- rct_design(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))
  expect_error(lmitt(y ~ 1, data = simdata, design = d),
               "Factor treatment")

  # factor numeric non 0/1
  simdata$foo <- as.factor(simdata$o)
  d <- rct_design(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))
  expect_error(lmitt(y ~ 1, data = simdata, design = d),
               "Factor treatment")

  # factor string levels
  simdata$foo <- factor(letters[simdata$o])
  d <- rct_design(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))
  expect_error(lmitt(y ~ 1, data = simdata, design = d),
               "Factor treatment")

  # ordinal 0/1
  simdata$foo <- as.ordered(simdata$z)
  d <- rct_design(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(d)[, 1]))
  expect_true(is.ordered(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))
  expect_error(lmitt(y ~ 1, data = simdata, design = d),
               "Ordered treatment")

  # ordinal numeric non 0/1
  simdata$foo <- as.ordered(simdata$o)
  d <- rct_design(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(d)[, 1]))
  expect_true(is.ordered(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))
  expect_error(lmitt(y ~ 1, data = simdata, design = d),
               "Ordered treatment")

  # ordinal string levels
  simdata$foo <- ordered(letters[simdata$o])
  d <- rct_design(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(d)[, 1]))
  expect_true(is.ordered(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))
  expect_error(lmitt(y ~ 1, data = simdata, design = d),
               "Ordered treatment")

  # character
  simdata$foo <- letters[simdata$o]
  d <- rct_design(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.character(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))

  # logical
  simdata$foo <- as.logical(simdata$z)
  d <- rct_design(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.logical(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))

  # conditional
  d <- rct_design(o > 2 ~ uoa(uoa1, uoa2), data = simdata)
  expect_true(is.logical(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))

  old_opt <- options()
  on.exit(options(old_opt))
  options(propertee_warn_on_conditional_treatment = FALSE)
  d <- rct_design(o > 2 ~ uoa(uoa1, uoa2), data = simdata)
  expect_true(is.logical(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))
  options(old_opt)

  # Non-conditional but bad variable name
  names(simdata)[5] <- "z==1"
  expect_error(d <- rct_design(`z==1`  ~ uoa(uoa1, uoa2),
                                              data = simdata),
               "isn't logical")


  # tests for other things...?

  z <- c(10485849600, 10477641600, 10561104000, 10562745600)
  d <- data.frame(z = as.Date(as.POSIXct(z, origin = "1582-10-14", tz = "GMT")),
                  unit = 1:4)

  expect_warning(rct_design(z ~ unitid(unit), data = d),
                 "STRONGLY")

})
