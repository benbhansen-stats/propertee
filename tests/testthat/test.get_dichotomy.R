test_that(".get_dichotomy with no lmitt.formula() and no argument", {
  expect_true(is.null(.get_dichotomy()))
})

test_that(".get_dichotomy with no lmitt.formula()", {
  expect_true(all.equal(.get_dichotomy(z ~ 1), z ~ 1))
})

test_that(".get_dichotomy with lmitt.formula()", {
  data(simdata)
  des <- rct_design(dose ~ cluster(uoa1, uoa2), simdata)
  expect_warning(
    mod1 <- lmitt(y ~ 1, design = des, data = simdata, dichotomy = dose > 50 ~ dose == 50,
                  weights = ate(des, dichotomy = dose <= 250 ~ dose > 250)),
    "passed to `.weights_calc()"
  )
  expect_equal(mod1$model$dose., as.numeric(simdata$dose > 50))
  expect_equal(mod1$model$`(weights)`, ate(des, dose <= 250 ~ dose > 250, data = simdata))
  
  mod2 <- lmitt(y ~ 1, design = des, data = simdata, dichotomy = dose <= 250 ~ dose > 250,
                weights = ate(des, dichotomy = dose <= 250 ~ dose > 250))
  expect_equal(mod2$model$dose., as.numeric(simdata$dose <= 250))
  expect_equal(mod2$model$`(weights)`, ate(des, dose <= 250 ~ dose > 250, data = simdata))
})

test_that(".get_dichotomy with weights", {
  data(simdata)
  des <- rct_design(dose ~ cluster(uoa1, uoa2), simdata)
  expect_warning(
    mod1 <- lm(y ~ assigned(dichotomy = dose > 50 ~ dose == 50), data = simdata,
               weights = ate(des, dichotomy = dose <= 250 ~ dose > 250)),
    "passed to `assigned()"
  )
  expect_equal(mod1$model[[2]], as.numeric(simdata$dose > 50))
  expect_equal(mod1$model$`(weights)`, ate(des, dose <= 250 ~ dose > 250, data = simdata))
  
  mod2 <- lm(y ~ assigned(dichotomy = dose <= 250 ~ dose > 250), data = simdata,
             weights = ate(des, dichotomy = dose <= 250 ~ dose > 250))
  expect_equal(mod2$model[[2]], as.numeric(simdata$dose <= 250))
  expect_equal(mod2$model$`(weights)`, ate(des, dose <= 250 ~ dose > 250, data = simdata))
})
