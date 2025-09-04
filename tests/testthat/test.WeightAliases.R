test_that("validWeights object", {
  expect_true(is.logical(.isValidWeightAlias("ate")))
  expect_true(is.logical(.isValidWeightTarget("ate")))
  expect_true(is.logical(.isValidWeightAlias("xxx")))
  expect_true(is.logical(.isValidWeightTarget("xxx")))

  expect_true(is.character(.listValidWeightAliases()))
  expect_true(is.character(.listValidWeightTargets()))

  expect_true(length(.validWeights) == 2)
  expect_true(all(sapply(.validWeights, length) >= 1))
})

test_that("weight aliases", {
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, spec = spec, weights = ate())
  mod2 <- lmitt(y ~ 1, data = simdata, spec = spec, weights = ett())
  mod3 <- lmitt(y ~ 1, data = simdata, spec = spec, weights = att())
  expect_false(isTRUE(all.equal(summary(mod1)$coeff, summary(mod2)$coeff)))
  expect_true(isTRUE(all.equal(summary(mod2)$coeff, summary(mod3)$coeff)))

  mod4 <- lmitt(y ~ 1, data = simdata, spec = spec, weights = "ate")
  mod5 <- lmitt(y ~ 1, data = simdata, spec = spec, weights = "ett")
  mod6 <- lmitt(y ~ 1, data = simdata, spec = spec, weights = "att")
  expect_false(isTRUE(all.equal(summary(mod4)$coeff, summary(mod5)$coeff)))
  expect_true(isTRUE(all.equal(summary(mod5)$coeff, summary(mod6)$coeff)))

})
