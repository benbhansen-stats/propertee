 test_that("lmitt", {

  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  da <- lmitt(y ~ 1, weights = ate(), data = simdata, design = des)

  expect_s4_class(da, "DirectAdjusted")
  expect_true(inherits(da, "lm"))

  da_ett <- lmitt(y ~ 1, weights = ett(), data = simdata, design = des)

  expect_s4_class(da_ett, "DirectAdjusted")
})

test_that("lmitt and lm return the same in simple cases", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid),
                   data = simdata)

  ml <- lm(y ~ z, data = simdata, weights = ate(des))
  da <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)
  ml2da <- lmitt(ml)

  expect_true(all(da$coef == ml$coef))
  expect_true(all(da$coef == ml2da$coef))
  expect_identical(da$model$`(weights)`, ml$model$`(weights)`)
  expect_identical(da$model$`(weights)`, ml2da$model$`(weights)`)

  ml <- lm(y ~ z, data = simdata, weights = ett(des))
  da <- lmitt(y ~ 1, data = simdata, weights = ett(), design = des)

  expect_true(all(da$coef == ml$coef))
  expect_identical(da$model$`(weights)`, ml$model$`(weights)`)

})

test_that("covariate adjustment", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  camod <- lm(y ~ x, data = simdata)

  da <- lmitt(y ~ 1, data = simdata, weights = ate(),
              offset = cov_adj(camod), design = des)
  expect_true(!is.null(da$model$"(offset)"))
  expect_true(inherits(da$model$"(offset)", "SandwichLayer"))

})


test_that("Design argument", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid1, cid2), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, design = des)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)
  mod3 <- lmitt(y ~ 1, data = simdata, design = z ~ cluster(cid1, cid2))
  expect_true(mod3@Design@type == "Obs")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

  des2 <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, design = des2)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des2)
  mod3 <- lmitt(y ~ 1, data = simdata,
                design = z ~ cluster(cid1, cid2) + forcing(force))
  expect_true(mod3@Design@type == "RD")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

  des3 <- obs_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, design = des3,
                subset = simdata$dose < 300)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(),
                subset = simdata$dose < 300, design = des3)
  mod3 <- lmitt(y ~ 1, data = simdata,
                design = z ~ cluster(cid1, cid2) + block(bid),
                subset = simdata$dose < 300)
  expect_true(mod3@Design@type == "Obs")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

})

test_that("Dichotomy option in Design creation", {

  data(simdata)
  des <- obs_design(dose ~ cluster(cid1, cid2), data = simdata,
                    dichotomy = dose > 200 ~ .)
  mod1 <- lmitt(y ~ assigned(), data = simdata, design = des)
  mod2 <- lmitt(y ~ assigned(), data = simdata,
                design = dose ~ cluster(cid1, cid2),
                dichotomy = dose > 200 ~ .)
  expect_identical(mod1@Design, mod2@Design)
  expect_true(all.equal(mod1$coefficients, mod2$coefficients,
                        check.attributes =  FALSE))


})

test_that("Allow non-binary treatment", {

  data(simdata)
  des <- obs_design(dose ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des)
  expect_true(any(!model.frame(mod1)$`assigned()` %in% 0:1))

})

test_that("Aliases for assigned", {

  data(simdata)
  des <- obs_design(dose ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lmitt(y ~ assigned(), data = simdata, design = des)
  mod2 <- lmitt(y ~ adopters(), data = simdata, design = des)
  mod3 <- lmitt(y ~ a.(), data = simdata, design = des)
  mod4 <- lmitt(y ~ z.(), data = simdata, design = des)
  expect_identical(mod1$coeff, mod2$coeff)
  expect_identical(mod1$coeff, mod3$coeff)
  expect_identical(mod1$coeff, mod4$coeff)
})

test_that("lmitt finds Design wherever it's stored", {
  data(simdata)
  des_form <- z ~ uoa(cid1, cid2) + block(bid)
  des <- obs_design(des_form, data = simdata)
  camod <- lm(y ~ x, data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, weights = ate(),
                offset = cov_adj(camod, design = des),
                design = des)
  mod2 <- lmitt(y ~ 1, design = des_form, data = simdata, weights = ate(),
                offset = cov_adj(camod))
  mod3 <- lmitt(y ~ 1, design = des, data = simdata, weights = ate(),
                offset = cov_adj(camod))

  expect_equal(coef(mod1), coef(mod2))
  expect_equal(coef(mod1), coef(mod3))
  expect_equal(vcovDA(mod1), vcovDA(mod2))
  expect_equal(vcovDA(mod1), vcovDA(mod3))
})

test_that("Allowed inputs to lmitt #73", {
  data(simdata)
  des <- obs_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)

  ### Allowed versions

  expect_no_error(l1 <- lmitt(y ~ 1, design = des, data = simdata))
  expect_no_error(l2 <- lmitt(y ~ dose, design = des, data = simdata))
  expect_no_error(l3 <- lmitt(y ~ 1 + dose, design = des, data = simdata))
  expect_no_error(l4 <- lmitt(y ~ dose + 1, design = des, data = simdata))

  ### equality of allowed versions

  expect_identical(l2$coeff, l3$coeff)
  expect_identical(l2$coeff, l4$coeff)

  ### Disallowed versions

  expect_error(lmitt(y ~ 0, design = des, data = simdata))
  expect_error(lmitt(y ~ dose + 0, design = des, data = simdata))
  expect_error(lmitt(y ~ 0 + dose, design = des, data = simdata))
  expect_error(lmitt(y ~ dose + bid, design = des, data = simdata))
  expect_error(lmitt(y ~ 1 + dose + bid, design = des, data = simdata))
  expect_error(lmitt(y ~ 1 + z, design = des, data = simdata))
  expect_error(lmitt(y ~ z, design = des, data = simdata))
  expect_error(lmitt(y ~ z + 1, design = des, data = simdata))



})
