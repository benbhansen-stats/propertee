test_that(paste("DirectAdjusted object created correctly with weights and no",
                "SandwichLayer in the lm call"), {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- new("DirectAdjusted",
              lm(y ~ z, data = simdata, weights = ate(des)),
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
              lm(y ~ z, data = simdata, weights = ate(des),
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

  dalm2 <- as.lmitt(lm(y ~ z, data = simdata, weights = ate(des),
                                offset = cov_adj(cmod)))
  expect_true(all.equal(dalm$coefficients,
                        dalm2$coefficients,
                        check.attributes = FALSE))

  # two identical adopters works silently
  dalm3 <- as.lmitt(lm(y ~ adopters() + adopters(), data = simdata,
                                weights = ate(des), offset = cov_adj(cmod)))
  expect_true(all.equal(dalm$coefficients,
                        dalm3$coefficients,
                        check.names = FALSE))

  # two different adopters fails
  expect_error(as.lmitt(lm(y ~ adopters(des) + adopters(),
                                    data = simdata,
                                    weights = ate(des),
                                    offset = cov_adj(cmod))),
               "Differing `adopters")
  # No treatment
  expect_error(as.lmitt(lm(y ~ x, data = simdata, weights = ate(des))),
               "Treatment z")

  sd2 <- simdata
  sd2$z <- NULL

  expect_error(as.lmitt(lm(y ~ z, data = sd2,
                                    weights = ate(des))),
               "'z' not found")
  dalm <- as.lmitt(lm(y ~ adopters(), data = sd2,
                               weights = ate(des)))

  expect_type(treatment_name(dalm), "character")
  expect_length(treatment_name(dalm), 1)
  expect_identical(treatment_name(dalm), "adopters()")
  expect_true(!is.na(coef(dalm)[treatment_name(dalm)]))

  des2 <- obs_design(o ~ cluster(cid2, cid1) + block(bid), data = simdata,
                     dichotomy = o > 2 ~ . )

  dalm2 <- as.lmitt(lm(y ~ adopters(), data = simdata,
                                weights = ate(des2)))

  expect_type(treatment_name(dalm2), "character")
  expect_length(treatment_name(dalm2), 1)
  expect_identical(treatment_name(dalm2), "adopters()")
  expect_true(!is.na(coef(dalm2)[treatment_name(dalm2)]))

  expect_error(as.lmitt(lm(y ~ o, data = simdata,
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
              Design = des)

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

  mod_da <- as.lmitt(mod)

  expect_s4_class(mod_da, "DirectAdjusted")
  expect_true(inherits(mod_da, "lm"))

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

  mod_da <- as.lmitt(mod)

  expect_s4_class(mod_da, "DirectAdjusted")
  expect_true(inherits(mod_da, "lm"))

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
  expect_error(as.lmitt(1),
               "lm object")
})

test_that("lm to DirectAdjusted fails without a Design object", {
  data(simdata)

  expect_error(as.lmitt(lm(y ~ z, data = simdata,
                                    weights = seq_len(nrow(simdata)))),
               "Cannot locate a `Design`")
})

test_that("vcov, confint, etc", {
  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- as.lmitt(lm(y ~ z, data = simdata, weights = ate(des)))

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

  expect_true(inherits(mod1$model$`(weights)`, "WeightedDesign"))
  expect_true(inherits(mod1$model$`(offset)`, "SandwichLayer"))

  # add subsetting
  mod2 <- lm(y ~ adopters(),
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

  mod <- lm(y ~ z,
            data = simdata,
            weights = ate(des),
            offset = cov_adj(lm(y ~ x, data = simdata),
                             design = des2))
  expect_error(as.lmitt(mod), "Multiple differing")

  expect_error(lmitt(y ~ z,
                     data = simdata,
                     weights = ate(des),
                     offset = cov_adj(lm(y ~ x, data = simdata),
                                      design = des2)))

  mod <- lm(y ~ z, data = simdata,
            offset = cov_adj(lm(y ~ x, data = simdata),
                             design = des2))

  expect_error(as.lmitt(mod, design = des), "Multiple differing")

  expect_error(lmitt(y ~ z, data = simdata,
                     offset = cov_adj(lm(y ~ x, data = simdata),
                                      design = des2),
                     design = des),
               "Multiple differing")

  # Checking that things work when passing multiple of the same designs
  mod1 <- lmitt(y ~ x, data = simdata, weights = ate(des))
  mod2 <- lmitt(y ~ x, data = simdata, weights = ate(des), design = des)
  expect_identical(mod1, mod2)

  mod3 <- lmitt(y ~ adopters(), data = simdata, weights = ate(des))
  mod4 <- lmitt(y ~ adopters(), data = simdata, weights = ate(des), design = des)
  expect_identical(mod3, mod4)


})

test_that("DirectAdjusted with SandwichLayer offset summary uses vcovDA SE's", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  damod <- lmitt(lm(y ~ z, data = simdata, weights = ate(des),
                    offset = cov_adj(cmod)))
  
  s <- summary(damod)
  expect_equal(s$coefficients[, 2L], sqrt(diag(vcovDA(damod))))
  expect_equal(s$coefficients[, 3L],
               damod$coefficients / sqrt(diag(vcovDA(damod))))
  expect_equal(s$coefficients[, 4L],
               2 * pt(abs(damod$coefficients / sqrt(diag(vcovDA(damod)))),
                      damod$df.residual,
                      lower.tail = FALSE))
  
  cmod <- lm(y ~ x, simdata, subset = z == 0)
  damod <- lmitt(lm(y ~ z, data = simdata, weights = ate(des),
                    offset = cov_adj(cmod)))
  s <- summary(damod)
  expect_equal(s$coefficients[, 2L], sqrt(diag(vcovDA(damod))))
  expect_equal(s$coefficients[, 3L],
               damod$coefficients / sqrt(diag(vcovDA(damod))))
  expect_equal(s$coefficients[, 4L],
               2 * pt(abs(damod$coefficients / sqrt(diag(vcovDA(damod)))),
                      damod$df.residual,
                      lower.tail = FALSE))
})

test_that("DirectAdjusted w/o SandwichLayer offset summary uses OLS SE's", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  damod <- lmitt(lm(y ~ z, data = simdata, weights = ate(des)))
  
  s <- summary(damod)
  expect_identical(s, do.call(getS3method("summary", "lm"), list(object = damod)))
})

test_that("vcov.DirectAdjusted handles vcovDA `type` arguments and non-SL offsets", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  damod1 <- lmitt(lm(y ~ z, data = simdata, weights = ate(des),
                     offset = cov_adj(cmod)))
  damod2 <- lmitt(lm(y ~ z, data = simdata, weights = ate(des)))
  
  vmat1 <- vcov(damod1)
  vmat2 <- vcov(damod1, type = "CR1")
  
  expect_error(vcov(damod1, type = "not_a_type"), "should be one of")
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
  damod1 <- lmitt(lm(y ~ z, data = simdata, weights = ate(des),
                     offset = cov_adj(cmod)))
  damod2 <- lmitt(lm(y ~ z, data = simdata, weights = ate(des)))
  
  expect_error(confint(damod1, type = "not_a_type"), "should be one of")
  
  vcovDA_ci.95 <- damod1$coefficients + sqrt(diag(vcovDA(damod1))) %o%
    qt(c(0.025, 0.975), damod1$df.residual)
  dimnames(vcovDA_ci.95) <- list(names(damod1$coefficients), c("2.5 %", "97.5 %"))
  ci1 <- confint(damod1, type = "CR1")
  ci2 <- confint(damod1)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcovDA_ci.95)

  vcovDA_ci.9 <- damod1$coefficients + sqrt(diag(vcovDA(damod1))) %o%
    qt(c(0.05, 0.95), damod1$df.residual)
  dimnames(vcovDA_ci.9) <- list(names(damod1$coefficients), c("5 %", "95 %"))
  ci1 <- confint(damod1, level = 0.9)
  expect_equal(ci1, vcovDA_ci.9)
  
  vcovDA_z.95 <- matrix(damod1$coefficients["z"] + sqrt(vcovDA(damod1)["z", "z"]) *
    qt(c(0.025, 0.975), damod1$df.residual), nrow = 1)
  dimnames(vcovDA_z.95) <- list(c("z"), c("2.5 %", "97.5 %"))
  ci1 <- confint(damod1, "z")
  ci2 <- confint(damod1, 2)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcovDA_z.95)

  vcovDA_z.9 <- matrix(damod1$coefficients["z"] + sqrt(vcovDA(damod1)["z", "z"]) *
    qt(c(0.05, 0.95), damod1$df.residual), nrow = 1)
  dimnames(vcovDA_z.9) <- list(c("z"), c("5 %", "95 %"))
  ci1 <- confint(damod1, "z", level = 0.9)
  ci2 <- confint(damod1, 2, level = 0.9)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcovDA_z.9)
  
  vcovlm.9 <- damod2$coefficients + sqrt(diag(do.call(getS3method("vcov", "lm"),
                                                      list(object = damod2)))) %o%
    qt(c(0.05, 0.95), damod2$df.residual)
  dimnames(vcovlm.9) <- list(names(damod2$coefficients), c("5 %", "95 %"))
  ci1 <- confint(damod2, level = 0.9)
  expect_equal(ci1, vcovlm.9)
  
  vcovlm_z.95 <- matrix(damod2$coefficients["z"] + sqrt(
    do.call(getS3method("vcov", "lm"), list(object = damod2))["z", "z"]) *
    qt(c(0.025, 0.975), damod2$df.residual), nrow = 1)
  dimnames(vcovlm_z.95) <- list(c("z"), c("2.5 %", "97.5 %"))
  ci1 <- confint(damod2, "z")
  ci2 <- confint(damod2, 2)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcovlm_z.95)
})

