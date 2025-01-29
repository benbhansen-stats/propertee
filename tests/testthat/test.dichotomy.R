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
