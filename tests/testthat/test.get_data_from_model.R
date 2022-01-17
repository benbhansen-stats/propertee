test_that("Design creation", {
  data(simdata)
  des <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata, weights = ate(des))

  # Error if multiple models in stack
  expect_error(lm(lm(y~x, data = simdata, weights = ate(des)) ~ 1),
               "Multiple")

  # glm is right now an "unsupported model", will need to update this
  # once we add glm support
  expect_warning(mod2 <- glm(y ~ x, data = simdata, weights = ate(des)),
                 "fallback")
  expect_equal(mod$coefficients,
               mod2$coefficients)


  expect_warning(mod3 <- with(simdata,
                              lm(y ~ x, weights = ate(des))),
                 "normal means")
  expect_equal(mod$coefficients,
               mod3$coefficients)

  # again, glm as an "unsupported model"
  expect_warning(mod4 <- with(simdata,
                              glm(y ~ x, weights = ate(des))),
                 "No supported models")
  expect_equal(mod$coefficients,
               mod3$coefficients)




})
