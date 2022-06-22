test_that("DirectAdjusted model created with invalid target", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)
  expect_error(new("DirectAdjusted",
                   lm(y ~ z, data = simdata, weights = ate(des)),
                   Design = des, target = "abc"),
               "must be one of")

})

test_that(paste("DirectAdjusted object created correctly with weights and no",
                "SandwichLayer in the lm call"), {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- new("DirectAdjusted",
              lm(y ~ z, data = simdata, weights = ate(des)),
              Design = des, target = "ett")

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(is(dalm, "lm"))

  expect_identical(dalm$model$"(weights)"@Design, des)
  expect_identical(dalm$model$"(weights)"@Design, dalm@Design)
})

test_that(paste("DirectAdjusted object created correctly with weights and ",
                "SandwichLayer in the lm call"), {
  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)
  cmod <- lm(y ~ x, data = simdata)
  dalm <- new("DirectAdjusted",
              lm(y ~ z, data = simdata, weights = ate(des),
                 offset = cov_adj(cmod)),
              Design = des, target = "ett")

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(is(dalm, "lm"))

  expect_equal(dalm$model$`(offset)`@.Data, as.numeric(cmod$fitted.values))
  expect_identical(dalm$model$`(offset)`@fitted_covariance_model, cmod)
  expect_equal(dalm$model$`(offset)`@prediction_gradient,
               stats::model.matrix(cmod))
  expect_identical(dalm$model$`(offset)`@keys, simdata[,var_names(des, "u")])
})

test_that("DA ensure treatment is found", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)
  cmod <- lm(y ~ x, data = simdata)

  dalm <- lmitt(y ~ adopters(), data = simdata, weights = ate(des),
                offset = cov_adj(cmod))

  expect_type(treatment_name(dalm), "character")
  expect_length(treatment_name(dalm), 1)
  expect_identical(treatment_name(dalm), "adopters()")
  expect_true(!is.na(coef(dalm)[treatment_name(dalm)]))

  dalm2 <- as.DirectAdjusted(lm(y ~ z, data = simdata, weights = ate(des),
                                offset = cov_adj(cmod)))
  expect_true(all.equal(dalm$coefficients,
                        dalm2$coefficients,
                        check.attributes = FALSE))

  # two identical adopters works silently
  dalm3 <- as.DirectAdjusted(lm(y ~ adopters() + adopters(), data = simdata,
                                weights = ate(des), offset = cov_adj(cmod)))
  expect_true(all.equal(dalm$coefficients,
                        dalm3$coefficients,
                        check.names = FALSE))

  # two different adopters fails
  expect_error(as.DirectAdjusted(lm(y ~ adopters(des) + adopters(),
                                    data = simdata,
                                    weights = ate(des),
                                    offset = cov_adj(cmod))),
               "Differing `adopters")
  # No treatment
  expect_error(as.DirectAdjusted(lm(y ~ x, data = simdata, weights = ate(des))),
               "Treatment z")

  sd2 <- simdata
  sd2$z <- NULL

  expect_error(as.DirectAdjusted(lm(y ~ z, data = sd2,
                                    weights = ate(des))),
               "'z' not found")
  dalm <- as.DirectAdjusted(lm(y ~ adopters(), data = sd2,
                               weights = ate(des)))

  expect_type(treatment_name(dalm), "character")
  expect_length(treatment_name(dalm), 1)
  expect_identical(treatment_name(dalm), "adopters()")
  expect_true(!is.na(coef(dalm)[treatment_name(dalm)]))

  des2 <- obs_design(o ~ cluster(cid2, cid1) + block(bid), data = simdata,
                     dichotomy = o > 2 ~ . )

  dalm2 <- as.DirectAdjusted(lm(y ~ adopters(), data = simdata,
                                weights = ate(des2)))

  expect_type(treatment_name(dalm2), "character")
  expect_length(treatment_name(dalm2), 1)
  expect_identical(treatment_name(dalm2), "adopters()")
  expect_true(!is.na(coef(dalm2)[treatment_name(dalm2)]))

  expect_error(as.DirectAdjusted(lm(y ~ o, data = simdata,
                                    weights = ate(des2))),
               "non-binary treatment")


  dalm_direct <- lmitt(y ~ adopters(), data = simdata, weights = ate(des2))

  expect_type(treatment_name(dalm_direct), "character")
  expect_length(treatment_name(dalm_direct), 1)
  expect_identical(treatment_name(dalm_direct), "adopters()")
  # internal name used in ittestimate
  expect_true(!is.na(coef(dalm_direct)[treatment_name(dalm_direct)]))
})

test_that("DirectAdjusted print/show", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)
  cmod <- lm(y ~ z, data = simdata)

  dalm <- new("DirectAdjusted",
              lm(y ~ z, data = simdata, weights = ate(des),
                 offset = cov_adj(cmod)),
              Design = des, target = "ett")

  aslm <- as(dalm, "lm")

  expect_silent(invisible(capture.output(expect_identical(print(dalm), dalm))))
  expect_silent(invisible(capture.output(expect_identical(show(dalm), dalm))))

  expect_output(print(dalm), "Coeff")
  expect_output(show(dalm), "Coeff")
})

test_that("lm to DirectAdjusted succeeds with weights and no SandwichLayer", {

  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  mod <- lm(y ~ z, data = simdata, weights = ate(des))

  mod_da <- as.DirectAdjusted(mod)

  expect_s4_class(mod_da, "DirectAdjusted")
  expect_true(is(mod_da, "lm"))

  expect_identical(mod_da$model$"(weights)"@Design, des)

  mod_lmitt <- lmitt(y ~ z, data = simdata, weights = ate(des))

  expect_true(all(mod_da$coef == mod_lmitt$coef))
  expect_identical(mod_da@Design, mod_lmitt@Design)
})

test_that("lm to DirectAdjusted with weights and SandwichLayer", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  cmod <- lm(y ~ x, data = simdata)

  mod <- lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod))

  mod_da <- as.DirectAdjusted(mod)

  expect_s4_class(mod_da, "DirectAdjusted")
  expect_true(is(mod_da, "lm"))

  expect_equal(mod_da$model$`(offset)`@.Data, as.numeric(cmod$fitted.values))
  expect_identical(mod_da$model$`(offset)`@fitted_covariance_model, cmod)
  expect_equal(mod_da$model$`(offset)`@prediction_gradient,
               stats::model.matrix(cmod))
  expect_identical(mod_da$model$`(offset)`@keys, simdata[,var_names(des, "u")])

  mod_lmitt <- lmitt(y ~ z, data = simdata, weights = ate(des),
                     offset = cov_adj(cmod))

  expect_true(all(mod_da$coef == mod_lmitt$coef))
  expect_identical(mod_da@Design, mod_lmitt@Design)
})

test_that("Conversion from lm to DirectAdjusted fails without an lm object", {
  expect_error(as.DirectAdjusted(1),
               "lm object")
})

test_that("lm to DirectAdjusted fails without a Design object", {
  data(simdata)

  expect_error(as.DirectAdjusted(lm(y ~ z, data = simdata,
                                    weights = seq_len(nrow(simdata)))),
               "Cannot locate `Design`")
})

test_that("Conversion from lm to DirectAdjusted fails without a target", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  camod <- lm(y ~ x, data = simdata)
  mod2 <- lm(y ~ z, data = simdata, offset = cov_adj(camod, design = des))

  expect_error(as.DirectAdjusted(mod2), "Cannot locate `target`")

  mod3 <- lm(y ~ z + offset(cov_adj(camod, design = des)), data = simdata)

  expect_error(as.DirectAdjusted(mod3), "Cannot locate `target`")
})

test_that("Conversion from lm to DirectAdjusted succeeds with a target", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  camod <- lm(y ~ x, data = simdata)
  mod2 <- lm(y ~ z, data = simdata, offset = cov_adj(camod, design = des))

  damod <- as.DirectAdjusted(mod2, target = "ate")
  expect_true(is(damod, "DirectAdjusted"))
  expect_identical(damod@target, "ate")

  mod3 <- lm(y ~ z + offset(cov_adj(camod, design = des)), data = simdata)

  damod <- as.DirectAdjusted(mod3, target = "ate")
  expect_true(is(damod, "DirectAdjusted"))
  expect_identical(damod@target, "ate")
})

test_that("vcov, confint, etc", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- as.DirectAdjusted(lm(y ~ z, data = simdata, weights = ate(des)))

  expect_true(is.matrix(vcov(dalm)))
  expect_equal(dim(vcov(dalm)), c(2, 2))

  expect_true(is.matrix(confint(dalm)))
  expect_equal(dim(confint(dalm)), c(2, 2))


})

test_that("subsetting model with weights and/or cov_adj", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid1, cid2), data = simdata)

  mod1 <- lm(y ~ adopters(),
             data = simdata,
             weights = ate(des),
             offset = cov_adj(lm(y ~ x, data = simdata)))

  expect_true(is(mod1$model$`(weights)`, "WeightedDesign"))
  expect_true(is(mod1$model$`(offset)`, "SandwichLayer"))

  # add subsetting
  mod2 <- lm(y ~ adopters(),
             data = simdata,
             weights = ate(des),
             offset = cov_adj(lm(y ~ x, data = simdata)),
             subset = simdata$dose < 300)

  expect_true(is(mod2$model$`(weights)`, "WeightedDesign"))
  expect_true(is(mod2$model$`(offset)`, "SandwichLayer"))

})
