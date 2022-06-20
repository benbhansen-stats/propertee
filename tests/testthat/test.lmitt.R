test_that("lmitt", {

  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  da <- lmitt(y ~ dose + x, weights = ate(des), data = simdata)

  expect_s4_class(da, "DirectAdjusted")
  expect_true(is(da, "lm"))
  expect_equal(da@target, "ate")

  da_ett <- lmitt(y ~ z + x, weights = ett(des), data = simdata)

  expect_equal(da_ett@target, "ett")
})

test_that("lmitt and lm return the same in simple cases", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid),
                   data = simdata)

  ml <- lm(y ~ z, data = simdata, weights = ate(des))
  da <- lmitt(y ~ z, data = simdata, weights = ate(des))

  expect_true(all(da$coef == ml$coef))
  expect_identical(da$model$`(weights)`, ml$model$`(weights)`)

  ml <- lm(y ~ z, data = simdata, weights = ett(des))
  da <- lmitt(y ~ z, data = simdata, weights = ett(des))

  expect_true(all(da$coef == ml$coef))
  expect_identical(da$model$`(weights)`, ml$model$`(weights)`)

})

test_that("covariate adjustment", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  camod <- lm(y ~ x, data = simdata)

  da <- lmitt(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(camod))
  expect_true(!is.null(da$model$"(offset)"))
  expect_true(is(da$model$"(offset)", "CovAdjPrediction"))

})


test_that("Design argument", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid1, cid2), data = simdata)

  mod1 <- lmitt(y ~ z, data = simdata, design = des, target = "ate")
  mod2 <- lmitt(y ~ z, data = simdata, weights = ate(des))
  mod3 <- lmitt(y ~ z, data = simdata, design = z ~ cluster(cid1, cid2),
                target = "ate")
  expect_true(mod3@Design@type == "Obs")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

  des2 <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), data = simdata)

  mod1 <- lmitt(y ~ z, data = simdata, design = des2, target = "ate")
  mod2 <- lmitt(y ~ z, data = simdata, weights = ate(des2))
  mod3 <- lmitt(y ~ z, data = simdata,
                design = z ~ cluster(cid1, cid2) + forcing(force),
                target = "ate")
  expect_true(mod3@Design@type == "RD")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

  des3 <- obs_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod1 <- lmitt(y ~ z, data = simdata, design = des3, target = "ate",
                subset = simdata$dose < 300)
  mod2 <- lmitt(y ~ z, data = simdata, weights = ate(des3),
                subset = simdata$dose < 300)
  mod3 <- lmitt(y ~ z, data = simdata,
                design = z ~ cluster(cid1, cid2) + block(bid),
                target = "ate", subset = simdata$dose < 300)
  expect_true(mod3@Design@type == "Obs")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

})
