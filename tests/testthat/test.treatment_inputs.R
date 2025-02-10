test_that("Treatment inputs", {
  data(simdata)

  # numeric 0/1 treatment
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.numeric(treatment(spec)[, 1]))
  expect_true(has_binary_treatment(spec))
  expect_silent(mod <- lmitt(y ~ 1, data = simdata, specification = spec))

  # numeric non-0/1
  spec <- rct_spec(o ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.numeric(treatment(spec)[, 1]))
  expect_false(has_binary_treatment(spec))
  expect_error(mod <- lmitt(y ~ 1, data = simdata, specification = spec),
               "continuous treatment")

  # factor 0/1
  simdata$foo <- as.factor(simdata$z)
  spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(spec)[, 1]))
  expect_false(has_binary_treatment(spec))
  expect_error(lmitt(y ~ 1, data = simdata, specification = spec),
               "Factor treatment")

  # factor numeric non 0/1
  simdata$foo <- as.factor(simdata$o)
  spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(spec)[, 1]))
  expect_false(has_binary_treatment(spec))
  expect_error(lmitt(y ~ 1, data = simdata, specification = spec),
               "Factor treatment")

  # factor string levels
  simdata$foo <- factor(letters[simdata$o])
  spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(spec)[, 1]))
  expect_false(has_binary_treatment(spec))
  expect_error(lmitt(y ~ 1, data = simdata, specification = spec),
               "Factor treatment")

  # ordinal 0/1
  simdata$foo <- as.ordered(simdata$z)
  spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(spec)[, 1]))
  expect_true(is.ordered(treatment(spec)[, 1]))
  expect_false(has_binary_treatment(spec))
  expect_error(lmitt(y ~ 1, data = simdata, specification = spec),
               "Ordered treatment")

  # ordinal numeric non 0/1
  simdata$foo <- as.ordered(simdata$o)
  spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(spec)[, 1]))
  expect_true(is.ordered(treatment(spec)[, 1]))
  expect_false(has_binary_treatment(spec))
  expect_error(lmitt(y ~ 1, data = simdata, specification = spec),
               "Ordered treatment")

  # ordinal string levels
  simdata$foo <- ordered(letters[simdata$o])
  spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.factor(treatment(spec)[, 1]))
  expect_true(is.ordered(treatment(spec)[, 1]))
  expect_false(has_binary_treatment(spec))
  expect_error(lmitt(y ~ 1, data = simdata, specification = spec),
               "Ordered treatment")

  # character
  simdata$foo <- letters[simdata$o]
  spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.character(treatment(spec)[, 1]))
  expect_false(has_binary_treatment(spec))

  # logical
  simdata$foo <- as.logical(simdata$z)
  spec <- rct_spec(foo ~ cluster(uoa1, uoa2), data = simdata)
  expect_true(is.logical(treatment(spec)[, 1]))
  expect_true(has_binary_treatment(spec))

  # conditional
  spec <- rct_spec(o > 2 ~ uoa(uoa1, uoa2), data = simdata)
  expect_true(is.logical(treatment(spec)[, 1]))
  expect_true(has_binary_treatment(spec))

  old_opt <- options()
  on.exit(options(old_opt))
  options(propertee_warn_on_conditional_treatment = FALSE)
  spec <- rct_spec(o > 2 ~ uoa(uoa1, uoa2), data = simdata)
  expect_true(is.logical(treatment(spec)[, 1]))
  expect_true(has_binary_treatment(spec))
  options(old_opt)

  # Non-conditional but bad variable name
  names(simdata)[5] <- "z==1"
  expect_error(spec <- rct_spec(`z==1`  ~ uoa(uoa1, uoa2),
                                data = simdata),
               "isn't logical")


  # tests for other things...?

  z <- c(10485849600, 10477641600, 10561104000, 10562745600)
  dat <- data.frame(z = as.Date(as.POSIXct(z, origin = "1582-10-14", tz = "GMT")),
                    unit = 1:4)

  expect_warning(rct_spec(z ~ unitid(unit), data = dat),
                 "STRONGLY")

})
