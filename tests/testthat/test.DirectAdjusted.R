test_that(paste("DirectAdjusted object created correctly with weights and no",
                "SandwichLayer in the lm call"), {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- new("DirectAdjusted",
              lm(y ~ assigned(), data = simdata, weights = ate(des)),
              Design = des,
              lmitt_fitted = TRUE)

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
              Design = des,
              lmitt_fitted = FALSE)

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
              Design = des,
              lmitt_fitted = TRUE)

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

  expect_true(all(round(mod_lmitt$coef, 4) %in%  round(mod_da$coef, 4)))
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
  des2 <- obs_design(o ~ cluster(cid1, cid2) + block(bid), data = simdata,
                     dichotomy = o >= 3 ~ .)

  # This should error, since there's no `design` argument to `lm`, `assigned()`
  # finds both the `des` inside `ate()` and `des2` inside `cov_adj()`
  ## expect_error(lm(y ~ assigned(),
  ##                 data = simdata,
  ##                 weights = ate(des),
  ##                 offset = cov_adj(lm(y ~ x, data = simdata),
  ##                                  design = des2)),
  ##              "Multiple differing")

  expect_error(lmitt(y ~ 1,
                        data = simdata,
                        design = des,
                        weights = ate(),
                        offset = cov_adj(lm(y ~ x, data = simdata),
                                         design = des2)),
               "Multiple differing")

  expect_error(lmitt(y ~ 1,
                        data = simdata,
                        design = des,
                        weights = ate(),
                        offset = cov_adj(lm(y ~ x, data = simdata),
                                         design = des2)),
               "Multiple differing")


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
  mod2 <- lmitt(lm(y ~ assigned(), simdata, weights = ate(des)))

  expect_false(identical(environment(), environment(formula(mod1))))
  expect_false(identical(environment(), environment(formula(mod2))))

  # `lmitt.formula` builds a data that includes the passed-in data, plus updated
  # RHS and LHS for centering and absorb.
  expect_true(all(simdata %in% environment(formula(mod1))$data))
  expect_true(all(environment(formula(mod2))$data %in% environment(formula(mod1))$data))
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

  uoas <- apply(simdata[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_"))
  vmat3 <- vcov(damod2)
  expect_true(all.equal(
    vmat3,
    sandwich::sandwich(damod2, meat. = sandwich::meatCL, cluster = uoas),
    check.attributes = FALSE))
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

  vcovDA_z.95 <- matrix(damod1$coefficients[1] +
                          sqrt(vcovDA(damod1)[1, 1]) *
                          qt(c(0.025, 0.975), damod1$df.residual), nrow = 1)
  ci1 <- confint(damod1)
  expect_true(all.equal(ci1, vcovDA_z.95, check.attributes = FALSE))

  vcovDA_z.9 <- matrix(damod1$coefficients[1] +
                         sqrt(vcovDA(damod1)[1, 1]) *
                         qt(c(0.05, 0.95), damod1$df.residual), nrow = 1)
  ci1 <- confint(damod1, level = 0.9)
  expect_true(all.equal(ci1, vcovDA_z.9, check.attributes = FALSE))

  uoas <- apply(simdata[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_"))
  vcovlm.9 <- damod2$coefficients + sqrt(diag(sandwich::sandwich(damod2,
                                                                 meat. = sandwich::meatCL,
                                                                 cluster = uoas))) %o%
    qt(c(0.05, 0.95), damod2$df.residual)
  ci1 <- confint(damod2, level = 0.9)
  expect_true(all.equal(ci1, vcovlm.9, check.attributes = FALSE))

  vcovlm_z.95 <- matrix(damod2$coefficients[1] +
                          sqrt(sandwich::sandwich(damod2,
                                                  meat. = sandwich::meatCL,
                                                  cluster = uoas)[1, 1]) *
                          qt(c(0.025, 0.975), damod2$df.residual), nrow = 1)
  ci1 <- confint(damod2)
  expect_true(all.equal(ci1, vcovlm_z.95, check.attributes = FALSE))
})

test_that("absorbed_intercepts", {
  data(simdata)

  blockeddes <- rct_design(z ~ block(bid) + cluster(cid1, cid2), data = simdata)
  noblocksdes <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  blocked_lmitt_fitted_absorbed <- lmitt(y ~ 1, data = simdata,
                                         design = blockeddes, absorb = TRUE)
  expect_message(blocked_lmitt_fitted_not_absorbed <-
                   lmitt(y ~ 1, data = simdata, design = blockeddes,
                         absorb = FALSE))
  blocked_not_lmitt_fitted <- as.lmitt(lm(y ~ assigned(blockeddes), data = simdata),
                                       design = blockeddes)

  expect_error(lmitt(y ~ 1, data = simdata, design = noblocksdes, absorb = TRUE))
  noblocks_lmitt_fitted_not_absorbed <- lmitt(y ~ 1, data = simdata,
                                              design = noblocksdes, absorb = FALSE)

  expect_true(blocked_lmitt_fitted_absorbed@absorbed_intercepts)
  expect_false(blocked_lmitt_fitted_not_absorbed@absorbed_intercepts)
  expect_false(blocked_not_lmitt_fitted@absorbed_intercepts)
  expect_false(noblocks_lmitt_fitted_not_absorbed@absorbed_intercepts)
})

test_that("absorbed_moderators", {
  data(simdata)

  blockeddes <- rct_design(z ~ block(bid) + cluster(cid1, cid2), data = simdata)
  noblocksdes <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  noblocks_lmitt_fittedsbgrp <- lmitt(y ~ force, data = simdata,
                                      design = noblocksdes)
  expect_message(blocked_lmitt_fittedsbgrp <- lmitt(y ~ force, data = simdata,
                                                    design = blockeddes))
  lmitt_fitted_nosbgrp <- lmitt(y ~ 1, data = simdata, design = noblocksdes)
  expect_message(blocked_lmitt_fitted_nosbgrp <- lmitt(y ~ 1, data = simdata,
                                                       design = blockeddes))
  not_lmitt_fitted <- as.lmitt(lm(y ~ assigned(noblocksdes), data = simdata),
                               design = noblocksdes)

  expect_equal(noblocks_lmitt_fittedsbgrp@absorbed_moderators, "force")
  expect_equal(blocked_lmitt_fittedsbgrp@absorbed_moderators, "force")
  expect_equal(lmitt_fitted_nosbgrp@absorbed_moderators, vector("character"))
  expect_equal(blocked_lmitt_fitted_nosbgrp@absorbed_moderators, vector("character"))
  expect_equal(not_lmitt_fitted@absorbed_moderators, vector("character"))
})

test_that("estfun.DirectAdjusted requires a certain model class", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des)
  mod@.S3Class <- "not_a_mod"

  expect_error(estfun(mod), "must have been fitted using")
})

test_that(paste("estfun.DirectAdjusted returns original psi if no offset or no",
                "SandwichLayer offset"), {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  ca <- cov_adj(cmod, newdata = simdata, design = des)

  nolmittmod1 <- lm(scale(y, scale = FALSE) ~ scale(z, scale = FALSE) + 0,
                    simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des)
  nolmittmod2 <- lm(scale(y, scale = FALSE) ~ scale(z, scale = FALSE) + 0,
                    data = simdata, offset = ca)
  mod2 <- lmitt(y ~ 1, data = simdata, design = des, offset = stats::predict(cmod))

  ef_expected1 <- estfun(nolmittmod1)
  ef_expected2 <- estfun(nolmittmod2)

  expect_true(all.equal(estfun(mod1), ef_expected1, check.attributes = FALSE))
  expect_true(all.equal(estfun(mod2), ef_expected2, check.attributes = FALSE))
})

test_that(paste("estfun.DirectAdjusted returns correct dimensions for full",
                "overlap of C and Q"), {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))

  expect_equal(dim(estfun(mod)), c(nrow(simdata), 1))
})

test_that(paste("estfun.DirectAdjusted returns correct dimensions for partial",
                "overlap of C and Q"), {
  data(simdata)
  set.seed(401)
  cmod_data <- rbind(
    simdata[simdata$cid1 == 1, c("x", "y", "cid1", "cid2")],
    data.frame("x" = rnorm(100), "y" = rnorm(100), "cid1" = NA, "cid2" = NA)
  )

  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))

  expect_equal(dim(estfun(mod)),
               c(nrow(simdata) + nrow(cmod_data[is.na(cmod_data$cid1),]), 1))
})

test_that(paste("estfun.DirectAdjusted returns correct dimensions for no",
                "overlap of C and Q"), {
  data(simdata)
  set.seed(401)
  cmod_data <- data.frame(
    "x" = rnorm(100), "y" = rnorm(100), "cid1" = NA, "cid2" = NA
  )

  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))

  expect_equal(dim(estfun(mod)), c(nrow(simdata) + nrow(cmod_data), 1))
})

test_that(paste("bread.DirectAdjusted returns bread.lm for DirectAdjusted objects",
                "with non-SandwichLayer or NULL offsets"), {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  ca <- predict(cmod)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  m1 <- lmitt(y ~ 1, data = simdata, design = des)
  m2 <- lmitt(y ~ 1, data = simdata, design = des, offset = ca)

  expected_out1 <- nrow(simdata) * chol2inv(m1$qr$qr)
  coef_names <- names(m1$coefficients)
  dimnames(expected_out1) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(m1), expected_out1)

  expected_out2 <- nrow(simdata) * chol2inv(m2$qr$qr)
  dimnames(expected_out2) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(m2), expected_out2)
})

test_that("bread.DirectAdjusted fails without a `qr` element", {
  data(simdata)
  simdata$uid <- seq_len(nrow(simdata))

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ uoa(uid), simdata)
  m <- lmitt(y ~ 1, data = simdata, design = des, weights = ate(des),
             offset = cov_adj(cmod))
  m$qr <- NULL

  expect_error(sandwich::bread(m), "Cannot compute")
})

test_that(paste("bread.DirectAdjusted returns expected output for full overlap",
                "of C and Q"), {
  data(simdata)
  simdata$uid <- seq_len(nrow(simdata))

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ uoa(uid), simdata)
  m <- lmitt(y ~ 1, data = simdata, design = des, weights = ate(des),
             offset = cov_adj(cmod))

  expected_out <- nrow(simdata) * chol2inv(m$qr$qr)
  coef_names <- names(m$coefficients)
  dimnames(expected_out) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(m), expected_out)
})

test_that(paste("bread.DirectAdjusted returns expected output for partial overlap",
                "of C and Q"), {
  set.seed(879)
  data(simdata)
  simdata$uid <- seq_len(nrow(simdata))

  cmod_data <- rbind(simdata[, c("uid", "x", "y")],
                     data.frame("uid" = seq(nrow(simdata) + 1,
                                            nrow(simdata) + 25),
                                "x" = rnorm(25),
                                "y" = rnorm(25)))

  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ uoa(uid), simdata)
  m <- lmitt(y ~ 1, data = simdata, design = des, weights = ate(des),
             offset = cov_adj(cmod))

  expected_out <- (nrow(simdata) + 25) * chol2inv(m$qr$qr)
  coef_names <- names(m$coefficients)
  dimnames(expected_out) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(m), expected_out)
})

test_that(paste("bread.DirectAdjusted returns expected output for no overlap",
                "of C and Q"), {
  set.seed(879)
  data(simdata)
  simdata$uid <- seq_len(nrow(simdata))

  cmod_data <- data.frame("uid" = seq(nrow(simdata) + 1, nrow(simdata) + 25),
                          "x" = rnorm(25),
                          "y" = rnorm(25))

  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ uoa(uid), simdata)
  m <- lmitt(y ~ 1, data = simdata, design = des, weights = ate(des),
             offset = cov_adj(cmod))

  expected_out <- (nrow(simdata) + 25) * chol2inv(m$qr$qr)
  coef_names <- names(m$coefficients)
  dimnames(expected_out) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(m), expected_out)
})

test_that(".make_uoa_ids fails without cluster argument or DirectAdjusted model", {
  data(simdata)
  mod <- lm(y ~ z, data = simdata)
  expect_error(.make_uoa_ids(mod), "Cannot deduce")
})

test_that(".make_uoa_ids returns correct ID's for non-DirectAdjusted model", {
  data(simdata)
  
  mod <- lm(y ~ z, data = simdata)
  expected_out <- factor(
    apply(simdata[, "cid1", drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  )
  
  expect_equal(.make_uoa_ids(mod, cluster = "cid1"), expected_out)
})

test_that(".make_uoa_ids returns correct ID's for non-SandwichLayer offset", {
  data(simdata)
  
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des)
  
  expected_out <- factor(
    apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  )
  expect_equal(.make_uoa_ids(mod), expected_out)
})

test_that(".make_uoa_ids returns correct ID's for full overlap of C and Q", {
  data(simdata)
  
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))
  
  expected_out <- factor(
    apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  )
  out <- .make_uoa_ids(dmod)
  
  expect_equal(out, expected_out)
})

test_that(".make_uoa_ids returns correct ID's for no overlap of C and Q", {
  data(simdata)
  set.seed(300)
  cmod_data1 <- data.frame("y" = rnorm(10), "x" = rnorm(10), "cid1" = NA, "cid2" = NA)
  cmod_data2 <- data.frame("y" = rnorm(10), "x" = rnorm(10),
                           "cid1" = rep(c(1, 2), each = 5), "cid2" = NA)
  
  cmod1 <- lm(y ~ x, cmod_data1)
  cmod2 <- lm(y ~ x, cmod_data2)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod1 <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod1))
  dmod2 <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod2))
  
  Q_uoas <- apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1,
                  function(...) paste(..., collapse = "_"))
  
  expect_warning(ids1 <- .make_uoa_ids(dmod1), "treated as IID")
  expect_warning(ids2 <- .make_uoa_ids(dmod2), "NA's for some but not all")
  
  expect_true(is.factor(ids1))
  expect_true(is.factor(ids2))

  expect_equal(length(ids1), nrow(simdata) + 10)
  expect_equal(length(ids2), nrow(simdata) + 10)

  expect_true(all.equal(ids1[1:nrow(simdata)], factor(Q_uoas), check.attributes = FALSE))
  expect_true(all.equal(ids2[1:nrow(simdata)], factor(Q_uoas), check.attributes = FALSE))

  expect_equal(length(unique(ids1)), length(unique(Q_uoas)) + 10)
  expect_equal(length(unique(ids2)), length(unique(Q_uoas)) + 2)
})

test_that(".make_uoa_ids returns correct ID's for partial overlap of C and Q", {
  data(simdata)
  set.seed(300)
  C_not_Q <- data.frame("y" = rnorm(20), "x" = rnorm(20), "cid1" = NA, "cid2" = NA)
  cmod_data <- rbind(simdata[, colnames(C_not_Q)], C_not_Q)
  
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))

  Q_uoas <- apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1,
                  function(...) paste(..., collapse = "_"))

  expect_warning(ids <- .make_uoa_ids(dmod), "treated as IID")
  
  expect_true(is.factor(ids))

  expect_equal(length(ids), nrow(simdata) + 20)

  expect_true(all.equal(ids[1:nrow(simdata)], factor(Q_uoas), check.attributes = FALSE))

  expect_equal(length(unique(ids)), length(unique(Q_uoas)) + 20)
})

test_that(paste(".sanitize_uoas fails without a DirectAdjusted object or",
                "SandwichLayer offset"), {
  data(simdata)
  
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lm(y ~ z, simdata, offset = predict(cmod))
  mod2 <- lmitt(y ~ 1, data = simdata, design = des)
  
  expect_error(.sanitize_uoas(mod1), "must be a DirectAdjusted object")
  expect_error(.sanitize_uoas(mod2), "must be a DirectAdjusted object")
})

test_that(paste(".sanitize_uoas succeeds with valid custom `cluster` argument",
                "(both full and no overlap cases)"), {
  set.seed(300)
  data(simdata)
  
  simdata$uid <- seq_len(nrow(simdata))
  cmod_data1 <- simdata
  cmod_data2 <- data.frame(x = rnorm(30), y = rnorm(30), cid1 = NA, cid2 = NA, uid = NA)
  cmod1 <- lm(y ~ x, cmod_data1)
  cmod2 <- lm(y ~ x, cmod_data2)
  des <- rct_design(z ~ unitid(cid1, cid2, uid), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod2))
  
  out1 <- .sanitize_uoas(mod1, cluster = c("cid1", "cid2"))
  Q_uoas <- apply(simdata[, c("cid1", "cid2")], 1,
                  function(...) paste(..., collapse = "_"))
  expect_warning(out2 <- .sanitize_uoas(mod2, cluster = c("cid1", "cid2")),
                 "Some or all rows")
  
  expect_equal(length(out1), nrow(simdata))
  expect_equal(length(out2), nrow(simdata) + 30)
  
  expect_equal(names(out1)[1:50], Q_uoas)
  expect_equal(names(out2)[1:50], Q_uoas)
  
  expect_equal(sum(out1 == "Q_C"), nrow(simdata))
  expect_equal(sum(out2 == "Q_C"), 0)
  
  expect_equal(sum(out1 == "C"), 0)
  expect_equal(sum(out2 == "C"), 30)
  
  expect_equal(sum(out1 == "Q"), 0)
  expect_equal(sum(out2 == "Q"), nrow(simdata))
})

test_that(paste(".sanitize_uoas generates correct sample assignments with partial",
                "overlap of Q and C"), {
  data(simdata)
  set.seed(300)
  C_not_Q <- data.frame("y" = rnorm(20), "x" = rnorm(20), "cid1" = NA, "cid2" = NA)
  cmod_data <- rbind(simdata[, colnames(C_not_Q)], C_not_Q)
  
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))
  
  expect_warning(out <- .sanitize_uoas(dmod), "Some or all rows")
  expect_equal(sum(out == "Q_C"), nrow(simdata))
  expect_equal(sum(out == "C"), 20)
  expect_equal(sum(out == "Q"), 0)
})

test_that(paste(".sanitize_uoas generates correct sample assignments with no",
                "overlap of Q and C"), {
  data(simdata)
  set.seed(300)
  cmod_data <- data.frame("y" = rnorm(10), "x" = rnorm(10), "cid1" = NA, "cid2" = NA)

  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))
  
  expect_warning(out <- .sanitize_uoas(dmod), "Some or all rows")
  expect_equal(sum(out == "Q"), nrow(simdata))
  expect_equal(sum(out == "C"), 10)
  expect_equal(sum(out == "Q_C"), 0)
})

test_that("sanitize_Q_uoas fails with invalid cluster argument", {
  data(simdata)
  mod <- lm(y ~ z, data = simdata)
  expect_error(.sanitize_Q_uoas(mod, cluster = "not_uoas"),
               "columns not_uoas in ITT effect model data")
  
  invalid_ids <- apply(simdata[, c("cid1", "cid2")], 1,
                       function(...) paste(..., collapse = "_"))
  expect_error(.sanitize_Q_uoas(mod, cluster = invalid_ids),
               ", 5_2 in ITT effect model data")
})

test_that("sanitize_Q_uoas succeeds with valid cluster argument", {
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  cmod <- lm(y ~ z, data = simdata)
  des <- rct_design(z ~ unitid(cid1, cid2, uid), simdata)
  dmod <- lmitt(y ~ 1, design = des, data = simdata, offset = cov_adj(cmod))
  out <- .sanitize_Q_uoas(dmod, cluster = c("cid1", "cid2"))

  expected_out <- apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1,
                        function(...) paste(..., collapse = "_"))
  expect_equal(out, expected_out)
})

test_that("checking proper errors in conversion from lm to DA", {
  data(simdata)
  des <- rct_design(z ~ unitid(cid1, cid2), simdata)

  expect_error(as.lmitt(lm(y ~ x, data = simdata), design = des),
               "aliases are not found")


  # Take in either Design or WeightedDesign
  mod1 <- flexida:::.convert_to_lmitt(lm(y ~ assigned(des), data = simdata),
                                     des,
                                     FALSE,
                                     TRUE,
                                     "a")
  mod2 <- flexida:::.convert_to_lmitt(lm(y ~ assigned(des), data = simdata),
                                     ate(des, data = simdata),
                                     FALSE,
                                     TRUE,
                                     "a")
  expect_identical(mod1@Design, mod2@Design)

  expect_error(flexida:::.convert_to_lmitt(lm(y ~ assigned(des), data = simdata),
                                     1,
                                     FALSE,
                                     TRUE,
                                     "a"), "must be a")


})
