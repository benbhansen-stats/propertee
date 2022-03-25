test_that("DirectAdjusted object", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- DirectAdjusted(lm(y ~ z, data = simdata, weights = ate(des)),
                         Design = des, target = "ett")

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(is(dalm, "lm"))

  expect_identical(dalm$model$"(weights)"@Design, des)
  expect_identical(dalm$model$"(weights)"@Design, dalm@Design)

  des <- obs_design(o ~ cluster(cid2, cid1) + block(bid), data = simdata)

  # issue18 - extra expect_warning
  expect_warning(dalm <- DirectAdjusted(lm(y ~ o, data = simdata,
                                           weights = ate(des)),
                                        Design = des, target = "ett"))
  # issue18 end

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(is(dalm, "lm"))

  expect_identical(dalm$model$"(weights)"@Design, des)
  expect_identical(dalm$model$"(weights)"@Design, dalm@Design)


  # issue18 - extra expect_warning
  expect_warning(expect_error(DirectAdjusted(lm(y ~ z, data = simdata, weights = ate(des)),
                                             Design = des, target = "abc"),
                              "must be one of"))
  # issue18 end
})

test_that("DirectAdjusted print/show", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- DirectAdjusted(lm(y ~ z, data = simdata, weights = ate(des)),
                         Design = des, target = "ett")

  aslm <- as(dalm, "lm")

  expect_silent(invisible(capture.output(expect_identical(print(dalm), dalm))))
  expect_silent(invisible(capture.output(expect_identical(show(dalm), dalm))))

  expect_output(print(dalm), "Coeff")
  expect_output(show(dalm), "Coeff")
})

test_that("Conversion from lm to DirectAdjusted", {

  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  mod <- lm(y~z, data = simdata, weights = ate(des))

  mod_da <- as.DirectAdjusted(mod)

  expect_s4_class(mod_da, "DirectAdjusted")

  mod_itt <- ittestimate(des, simdata, "y")

  expect_true(all(mod_da$coef == mod_itt$coef))
  expect_identical(mod_da@Design, mod_itt@Design)

  expect_error(as.DirectAdjusted(1),
               "lm object")
  expect_error(as.DirectAdjusted(lm(y ~ z, data = simdata, weights = seq_len(nrow(simdata)))),
               "WeightedDesign weights")

  expect_error(as.DirectAdjusted(lm(y ~ x, data = simdata, weights = ate(des))),
               "treatment")

  expect_error(as.DirectAdjusted(lm(y ~ 1, data = simdata, weights = ate(des))),
               "treatment")

  expect_error(as.DirectAdjusted(lm(y ~ rep(0, nrow(simdata)), data = simdata, weights = ate(des))),
               "constant")
})

test_that("vcov, confint, etc", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- DirectAdjusted(lm(y ~ z, data = simdata, weights = ate(des)),
                         Design = des, target = "ett")

  expect_true(is.matrix(vcov(dalm)))
  expect_equal(dim(vcov(dalm)), c(2, 2))

  expect_true(is.matrix(confint(dalm)))
  expect_equal(dim(confint(dalm)), c(2, 2))


})
