test_that("Design creation", {
  data(simdata)
  des <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata, weights = ate(des))

  #Handling multiple models in the stack
  mod1 <- lm(predict(lm(y~x, data = simdata, weights = ate(des))) ~ 1)
  mod2 <- lm(predict(lm(y~x, data = simdata)) ~ 1, weights = ate(des), data = simdata)
  expect_true(mod1$coefficients[1] != mod2$coefficients[1])

  # linear glm and lm are equivalent
  mod3 <- glm(y ~ x, data = simdata, weights = ate(des))
  expect_equal(mod$coefficients,
               mod3$coefficients)

  # Works using `with`
  mod4 <- with(simdata,
               lm(y ~ x, weights = ate(des)))
  expect_equal(mod$coefficients,
               mod4$coefficients)

  # Works using `attach`
  attach(simdata)
  mod5 <- lm(y ~ x, weights = ate(des))
  expect_equal(mod$coefficients,
               mod5$coefficients)
  detach(simdata)

  # Allow simple manipulation of weights
  mod6 <- lm(y ~ x, data = simdata, weights = sqrt(ate(des)))


  # Test for fallback if no model.frame call

  expect_warning(with(simdata, ate(des)), "No call")
})
