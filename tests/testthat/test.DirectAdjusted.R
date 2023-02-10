test_that(paste("DirectAdjusted object created correctly with weights and no",
                "SandwichLayer in the lm call"), {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- new("DirectAdjusted",
              lm(y ~ assigned(), data = simdata, weights = ate(des)),
              Design = des)

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(inherits(dalm, "lm"))

  expect_identical(dalm$model$"(weights)"@Design, des)
  expect_identical(dalm$model$"(weights)"@Design, dalm@Design)
})

test_that(paste("DirectAdjusted object created correctly with weights and ",
                "SandwichLayer in the lm call"), {
  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)
  cmod <- lm(y ~ x, data = simdata)
  dalm <- new("DirectAdjusted",
              lm(y ~ assigned(), data = simdata, weights = ate(des),
                 offset = cov_adj(cmod)),
              Design = des)

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(inherits(dalm, "lm"))

  expect_equal(dalm$model$`(offset)`@.Data, as.numeric(cmod$fitted.values))
  expect_identical(dalm$model$`(offset)`@fitted_covariance_model, cmod)
  expect_equal(dalm$model$`(offset)`@prediction_gradient,
               stats::model.matrix(cmod))
  expect_identical(dalm$model$`(offset)`@keys, simdata[,var_names(des, "u")])
})

test_that("DirectAdjusted print/show", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)
  cmod <- lm(y ~ z, data = simdata)

  dalm <- new("DirectAdjusted",
              lm(y ~ assigned(), data = simdata, weights = ate(des),
                 offset = cov_adj(cmod)),
              Design = des)

  aslm <- as(dalm, "lm")

  expect_silent(invisible(capture.output(expect_identical(print(dalm), dalm))))
  expect_silent(invisible(capture.output(expect_identical(show(dalm), dalm))))

  # Expect "assigned()"
  expect_output(print(dalm), "assigned()")
  expect_output(show(dalm), "assigned()")
})

test_that("lm to DirectAdjusted succeeds with weights and no SandwichLayer", {

  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  mod <- lm(y ~ assigned(), data = simdata, weights = ate(des))

  mod_da <- as.lmitt(mod)

  expect_s4_class(mod_da, "DirectAdjusted")
  expect_true(inherits(mod_da, "lm"))

  expect_identical(mod_da$model$"(weights)"@Design, des)

  mod_lmitt <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)

  expect_true(all(round(mod_lmitt$coef, 4) %in% round(mod_da$coef, 4)))
  expect_identical(mod_da@Design, mod_lmitt@Design)
})

test_that("lm to DirectAdjusted with weights and SandwichLayer", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  cmod <- lm(y ~ x, data = simdata)

  mod <- lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))

  mod_da <- as.lmitt(mod)

  expect_s4_class(mod_da, "DirectAdjusted")
  expect_true(inherits(mod_da, "lm"))

  expect_equal(mod_da$model$`(offset)`@.Data, as.numeric(cmod$fitted.values))
  expect_identical(mod_da$model$`(offset)`@fitted_covariance_model, cmod)
  expect_equal(mod_da$model$`(offset)`@prediction_gradient,
               stats::model.matrix(cmod))
  expect_identical(mod_da$model$`(offset)`@keys, simdata[,var_names(des, "u")])

  mod_lmitt <- lmitt(y ~ 1, data = simdata, weights = ate(),
                     offset = cov_adj(cmod), design = des)

  expect_true(all(mod_da$coef == mod_lmitt$coef))
  expect_identical(mod_da@Design, mod_lmitt@Design)
})

test_that("Conversion from lm to DirectAdjusted fails without an lm object", {
  expect_error(as.lmitt(1),
               "lm object")
})

test_that("lm to DirectAdjusted fails without a Design object", {
  data(simdata)

  expect_error(as.lmitt(lm(y ~ assigned(), data = simdata,
                                    weights = seq_len(nrow(simdata)))),
               "Unable to locate Design")
})

test_that("vcov, confint, etc", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))

  expect_true(is.matrix(vcov(dalm)))
  expect_equal(dim(vcov(dalm)), c(2, 2))

  expect_true(is.matrix(confint(dalm)))
  expect_equal(dim(confint(dalm)), c(2, 2))


})

test_that("subsetting model with weights and/or cov_adj", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid1, cid2), data = simdata)

  mod1 <- lm(y ~ assigned(),
             data = simdata,
             weights = ate(des),
             offset = cov_adj(lm(y ~ x, data = simdata)))

  expect_true(inherits(mod1$model$`(weights)`, "WeightedDesign"))
  expect_true(inherits(mod1$model$`(offset)`, "SandwichLayer"))

  # add subsetting
  mod2 <- lm(y ~ assigned(),
             data = simdata,
             weights = ate(des),
             offset = cov_adj(lm(y ~ x, data = simdata)),
             subset = simdata$dose < 300)

  expect_true(inherits(mod2$model$`(weights)`, "WeightedDesign"))
  expect_true(inherits(mod2$model$`(offset)`, "SandwichLayer"))

})

test_that("Differing designs found", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid1, cid2), data = simdata)
  des2 <- obs_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  expect_error(lm(y ~ assigned(),
                  data = simdata,
                  weights = ate(des),
                  offset = cov_adj(lm(y ~ x, data = simdata),
                                   design = des2)),
               "Multiple differing")

  expect_error(lmitt(y ~ 1,
                     data = simdata,
                     design = des,
                     weights = ate(des),
                     offset = cov_adj(lm(y ~ x, data = simdata),
                                      design = des2)))

  mod <- lm(y ~ assigned(), data = simdata,
            offset = cov_adj(lm(y ~ x, data = simdata),
                             design = des2))

  expect_error(as.lmitt(mod, design = des), "Multiple differing")

  expect_error(lmitt(y ~ 1, data = simdata,
                     offset = cov_adj(lm(y ~ x, data = simdata),
                                      design = des2),
                     design = des),
               "Multiple differing")

  # Checking that things work when passing multiple of the same designs
  mod1 <- lmitt(y ~ x, data = simdata, weights = ate(), design = des)
  mod2 <- lmitt(y ~ x, data = simdata, weights = ate(des), design = des)
  expect_equal(mod1$coefficients, mod2$coefficients)

  mod3 <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)
  mod4 <- lmitt(y ~ 1, data = simdata, weights = ate(des), design = des)
  expect_equal(mod3$coefficients, mod4$coefficients)


})

test_that("DirectAdjusted object has its own evaluation environment", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des)
  mod2 <- lmitt(lm(y ~ 1, simdata, weights = ate(des)))

  expect_false(identical(environment(), environment(formula(mod1))))
  expect_false(identical(environment(), environment(formula(mod2))))
  expect_equal(environment(formula(mod1))$data, simdata)
  expect_equal(environment(formula(mod1))$data, environment(formula(mod2))$data)
  expect_equal(environment(formula(mod1))$design, des)
  expect_equal(environment(formula(mod1))$design, environment(formula(mod2))$design)
})

test_that("vcov.DirectAdjusted handles vcovDA `type` arguments and non-SL offsets", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  damod1 <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des),
                     offset = cov_adj(cmod)))
  damod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)

  vmat1 <- vcov(damod1)
  vmat2 <- vcov(damod1, type = "CR0")

  expect_error(vcov(damod1, type = "not_a_type"), "should be")
  expect_identical(vmat1, vmat2)
  expect_identical(vmat1, vcovDA(damod1))

  vmat3 <- vcov(damod2)
  expect_identical(vmat3, do.call(getS3method("vcov", "lm"),
                                  list(object = damod2)))
})

test_that("confint.DirectAdjusted handles vcovDA `type` arguments and non-SL offsets", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  damod1 <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des,
                     offset = cov_adj(cmod))
  damod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)

  expect_error(confint(damod1, type = "not_a_type"), "should be")

  vcovDA_ci.95 <- damod1$coefficients + sqrt(diag(vcovDA(damod1))) %o%
    qt(c(0.025, 0.975), damod1$df.residual)
  dimnames(vcovDA_ci.95) <- list(names(damod1$coefficients), c("2.5 %", "97.5 %"))
  ci1 <- confint(damod1, type = "CR0")
  ci2 <- confint(damod1)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcovDA_ci.95)

  vcovDA_ci.9 <- damod1$coefficients + sqrt(diag(vcovDA(damod1))) %o%
    qt(c(0.05, 0.95), damod1$df.residual)
  dimnames(vcovDA_ci.9) <- list(names(damod1$coefficients), c("5 %", "95 %"))
  ci1 <- confint(damod1, level = 0.9)
  expect_equal(ci1, vcovDA_ci.9)

  vcovDA_z.95 <- matrix(damod1$coefficients["assigned()"] +
                          sqrt(vcovDA(damod1)["assigned()", "assigned()"]) *
                          qt(c(0.025, 0.975), damod1$df.residual), nrow = 1)
  dimnames(vcovDA_z.95) <- list(c("assigned()"), c("2.5 %", "97.5 %"))
  ci1 <- confint(damod1, "assigned()")
  ci2 <- confint(damod1, 2)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcovDA_z.95)

  vcovDA_z.9 <- matrix(damod1$coefficients["assigned()"] +
                         sqrt(vcovDA(damod1)["assigned()", "assigned()"]) *
                         qt(c(0.05, 0.95), damod1$df.residual), nrow = 1)
  dimnames(vcovDA_z.9) <- list(c("assigned()"), c("5 %", "95 %"))
  ci1 <- confint(damod1, "assigned()", level = 0.9)
  ci2 <- confint(damod1, 2, level = 0.9)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcovDA_z.9)

  vcovlm.9 <- damod2$coefficients + sqrt(diag(do.call(getS3method("vcov", "lm"),
                                                      list(object = damod2)))) %o%
    qt(c(0.05, 0.95), damod2$df.residual)
  dimnames(vcovlm.9) <- list(names(damod2$coefficients), c("5 %", "95 %"))
  ci1 <- confint(damod2, level = 0.9)
  expect_equal(ci1, vcovlm.9)

  vcovlm_z.95 <- matrix(damod2$coefficients["assigned()"] +
                          sqrt(do.call(getS3method("vcov", "lm"),
                                       list(object = damod2))["assigned()", "assigned()"]) *
                          qt(c(0.025, 0.975), damod2$df.residual), nrow = 1)
  dimnames(vcovlm_z.95) <- list(c("assigned()"), c("2.5 %", "97.5 %"))
  ci1 <- confint(damod2, "assigned()")
  ci2 <- confint(damod2, 2)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcovlm_z.95)
})
