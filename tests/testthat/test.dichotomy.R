test_that("find and validate dichotomies with no upstream dichotomies and list of NULLs argument", {
  expect_true(is.null(.validate_dichotomy(list(lmitt = NULL, weights = NULL))))
})

test_that("find and validate dichotomies with no upstream dichotomies and non-null argument", {
  expect_identical(deparse1(.validate_dichotomy(z ~ 1)), deparse1(z ~ 1))
})

test_that("find and validate dichotomies doesn't get valid option", {
  expect_error(.validate_dichotomy(character()), "formulas provided")
})

test_that("find and validate dichotomies with upstream dichotomy in lmitt.formula()", {
  data(simdata)
  spec <- rct_spec(dose ~ cluster(uoa1, uoa2), simdata)
  expect_warning(
    expect_warning(
      expect_warning(
        mod1 <- lmitt(y ~ 1, specification = spec, data = simdata, dichotomy = dose > 50 ~ dose == 50,
                      weights = ate(spec, dichotomy = dose <= 250 ~ dose > 250)),
        "passed to `.weights_calc()"
      ),
      "passed to `a.()"
    ),
    "passed to `a.()"
  )
  expect_equal(mod1$model$dose., as.numeric(simdata$dose > 50))
  expect_equal(mod1$model$`(weights)`, ate(spec, dose <= 250 ~ dose > 250, data = simdata))

  mod2 <- lmitt(y ~ 1, specification = spec, data = simdata, dichotomy = dose <= 250 ~ dose > 250,
                weights = ate(spec, dichotomy = dose <= 250 ~ dose > 250))
  expect_equal(mod2$model$dose., as.numeric(simdata$dose <= 250))
  expect_equal(mod2$model$`(weights)`, ate(spec, dose <= 250 ~ dose > 250, data = simdata))

  mod3 <- lmitt(y ~ 1, specification = spec, data = simdata, dichotomy = dose <= 250 ~ dose > 250,
                weights = ate(spec))
  expect_equal(mod3$model$dose., mod2$model$dose.)
  expect_equal(mod3$model$`(weights)`@.Data, mod2$model$`(weights)`@.Data)
})

test_that("find and validate dichotomies with no upstream lmitt.formula()", {
  data(simdata)
  spec <- rct_spec(dose ~ cluster(uoa1, uoa2), simdata)
  expect_warning(
    mod1 <- lm(y ~ assigned(dichotomy = dose > 50 ~ dose == 50), data = simdata,
               weights = ate(spec, dichotomy = dose <= 250 ~ dose > 250)),
    "passed to `assigned()"
  )
  expect_equal(mod1$model[[2]], as.numeric(simdata$dose > 50))
  expect_equal(mod1$model$`(weights)`@.Data,
               ate(spec, dose <= 250 ~ dose > 250, data = simdata)@.Data)
  expect_equal(mod1$model$`(weights)`@target,
               ate(spec, dose <= 250 ~ dose > 250, data = simdata)@target)
  expect_equal(deparse1(mod1$model$`(weights)`@dichotomy),
               deparse1(ate(spec, dose <= 250 ~ dose > 250, data = simdata)@dichotomy))


  mod2 <- lm(y ~ assigned(dichotomy = dose <= 250 ~ dose > 250), data = simdata,
             weights = ate(spec, dichotomy = dose <= 250 ~ dose > 250))
  expect_equal(mod2$model[[2]], as.numeric(simdata$dose <= 250))
  expect_equal(mod2$model$`(weights)`@.Data,
               ate(spec, dose <= 250 ~ dose > 250, data = simdata)@.Data)
  expect_equal(mod2$model$`(weights)`@target,
               ate(spec, dose <= 250 ~ dose > 250, data = simdata)@target)
  expect_equal(deparse1(mod2$model$`(weights)`@dichotomy),
               deparse1(ate(spec, dose <= 250 ~ dose > 250, data = simdata)@dichotomy))

  expect_error(
    mod3 <- lm(y ~ assigned(dichotomy = dose > 50 ~ dose == 50), data = simdata,
               weights = ate(spec)),
    "Must provide a dichotomy"
  )
})

test_that("no warning when dichotomy is passed as a stored object", {
  set.seed(33)
  uoadata <- data.frame(cid = seq_len(8),
                        bid = c(rep(seq_len(2), 2), rep(1, 2), rep(2, 2)),
                        dose = c(rep(0, 4), rep(1, 2), rep(2, 2)))
  individdata <- data.frame(cid = rep(seq_len(8), each = 10), y = rnorm(80))
  specification <- rct_spec(dose ~ uoa(cid) + block(bid), data = uoadata)
  
  # store the dichotomy in an object and create weights ahead of time with stored dichotomy
  dich <- dose == 1 ~ dose == 0
  wts <- ett(specification, data = individdata, dichotomy = dich)
  
  # first, check when dichotomy is passed directly and through the weights
  expect_silent(om1 <- lmitt(y~1, specification, individdata, weights = wts, dichotomy = dich))
  
  # then check when only passed through the weights
  expect_silent(om2 <- lmitt(y~1, specification, individdata,
                             weights = ett(specification, data = individdata, dichotomy = dich)))
  
  # also, check when the dichotomy or an object storing a dichotomy is passed to a function
  fit_lmitt <- function(data, spec, d) {
    wts <- ett(spec, data = data, dichotomy = d)
    m1 <- lmitt(y~1, spec, data, weights = wts, dichotomy = d)
    m1
  }
  expect_silent(om3 <- fit_lmitt(individdata, specification, dose == 1 ~ dose == 0))
  expect_silent(om4 <- fit_lmitt(individdata, specification, dich))
  
  expect_equal(om1$coef, om2$coef)
  expect_equal(om1$coef, om3$coef)
  expect_equal(om1$coef, om4$coef)
  
  # and finally, non-lmitt calls
  expect_silent(om5 <- lm(y ~ assigned(specification, individdata, dich), individdata, weights = wts))
  expect_equal(unname(om1$coef[1:2]), unname(om5$coef))
})

test_that("proceed with NULL dichotomy when passed as a named object", {
  fit_lmitt <- function(data, spec, dich) {
    m1 <- lmitt(y~1, spec, data, dichotomy = dich)
    m1
  }
  
  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)), y = rnorm(10))
  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)
  
  expect_silent(om1 <- lmitt(y~1, spec, testdata))
  expect_silent(om2 <- lmitt(y~1, spec, testdata, dichotomy = NULL))
  expect_silent(om3 <- fit_lmitt(testdata, spec, z == 1 ~ z == 0))
  expect_silent(om4 <- fit_lmitt(testdata, spec, NULL))
  expect_equal(om1$coefficients, om2$coefficients)
  expect_equal(om1$coefficients, om3$coefficients)
  expect_equal(om1$coefficients, om4$coefficients)
})

test_that("fail when dichotomy is a named but non-existent object", {
  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)), y = rnorm(10))
  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)
  expect_error(ett(spec, data = testdata, dichotomy = quote(d)),
               "Could not find d in call stack")
})
