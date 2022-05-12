test_that("DirectAdjusted object", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- new("DirectAdjusted",
              lm(y ~ z, data = simdata, weights = ate(des)),
              Design = des, target = "ett")

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(is(dalm, "lm"))

  expect_identical(dalm$model$"(weights)"@Design, des)
  expect_identical(dalm$model$"(weights)"@Design, dalm@Design)

  expect_error(new("DirectAdjusted",
                   lm(y ~ z, data = simdata, weights = ate(des)),
                   Design = des, target = "abc"),
                              "must be one of")

})

test_that("DA ensure treatment is found", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- as.DirectAdjusted(lm(y ~ z, data = simdata, weights = ate(des)))

  expect_type(treatment(dalm), "character")
  expect_length(treatment(dalm), 1)
  expect_identical(treatment(dalm), var_names(des, "t"))
  expect_true(!is.na(coef(dalm)[treatment(dalm)]))

  dalm <- as.DirectAdjusted(lm(y ~ adopters(), data = simdata,
                               weights = ate(des)))

  expect_type(treatment(dalm), "character")
  expect_length(treatment(dalm), 1)
  expect_identical(treatment(dalm), "adopters()")
  expect_true(!is.na(coef(dalm)[treatment(dalm)]))

  des2 <- obs_design(o ~ cluster(cid2, cid1) + block(bid), data = simdata,
                     dichotomy = o > 2 ~ . )

  dalm2 <- as.DirectAdjusted(lm(y ~ adopters(), data = simdata,
                                weights = ate(des2)))

  expect_type(treatment(dalm2), "character")
  expect_length(treatment(dalm2), 1)
  expect_identical(treatment(dalm2), "adopters()")
  expect_true(!is.na(coef(dalm2)[treatment(dalm2)]))

  expect_error(as.DirectAdjusted(lm(y ~ o, data = simdata,
                                    weights = ate(des2))),
               "treatment not found")


  daitt <- ittestimate(des, simdata, "y")

  expect_type(treatment(daitt), "character")
  expect_length(treatment(daitt), 1)
  expect_identical(treatment(daitt), "z__")
  # internal name used in ittestimate
  expect_true(!is.na(coef(daitt)[treatment(daitt)]))
})

test_that("DirectAdjusted print/show", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- new("DirectAdjusted",
              lm(y ~ z, data = simdata, weights = ate(des)),
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

  mod <- lm(y ~ z, data = simdata, weights = ate(des))

  mod_da <- as.DirectAdjusted(mod)

  expect_s4_class(mod_da, "DirectAdjusted")

  mod_itt <- ittestimate(des, simdata, "y")

  expect_true(all(mod_da$coef == mod_itt$coef))
  expect_identical(mod_da@Design, mod_itt@Design)

  expect_error(as.DirectAdjusted(1),
               "lm object")
  expect_error(as.DirectAdjusted(lm(y ~ z, data = simdata,
                                    weights = seq_len(nrow(simdata)))),
               "Cannot locate `Design`")

  expect_error(as.DirectAdjusted(lm(y ~ rep(0, nrow(simdata)),
                                    data = simdata, weights = ate(des))),
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
