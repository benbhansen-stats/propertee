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
  expect_identical(dalm$model$`(offset)`@keys,
                   cbind(simdata[,var_names(des, "u")], in_Q = rep(TRUE, nrow(simdata))))
})

test_that("DirectAdjusted print/show", {

  data(simdata)
  des <- obs_design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)
  cmod <- lm(y ~ z, data = simdata)

  ######## Pretend it comes from `as.lmitt`
  dalm <- new("DirectAdjusted",
              lm(y ~ assigned(), data = simdata, weights = ate(des),
                 offset = cov_adj(cmod)),
              Design = des,
              lmitt_fitted = FALSE)

  aslm <- as(dalm, "lm")

  expect_silent(invisible(capture.output(expect_identical(print(dalm), dalm))))
  expect_silent(invisible(capture.output(expect_identical(show(dalm), dalm))))

  # Expect "assigned()"
  expect_output(print(dalm), "assigned()")
  expect_output(show(dalm), "assigned()")

  ######## Now pretend it comes from `lmitt.formula`
  # Need `txt_` to match what comes out of a call like `lmitt(y ~ 1...)`
  simdata$z. <- simdata$z
  dalm <- new("DirectAdjusted",
              lm(y ~ z., data = simdata, weights = ate(des),
                 offset = cov_adj(cmod)),
              Design = des,
              lmitt_fitted = TRUE)

  aslm <- as(dalm, "lm")

  expect_silent(invisible(capture.output(expect_identical(print(dalm), dalm))))
  expect_silent(invisible(capture.output(expect_identical(show(dalm), dalm))))

  # Expect "assigned()"
  expect_output(print(dalm), "z\\.")
  expect_output(show(dalm), "z\\.")
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
  expect_identical(mod_da$model$`(offset)`@keys,
                   cbind(simdata[,var_names(des, "u")], in_Q = rep(TRUE, nrow(simdata))))

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


  des3 <- obs_design(o ~ uoa(cid1, cid2), dat = simdata)

  # No error with a dichotomy
  expect_silent(mod <- lmitt(y ~ 1, data = simdata, design = des3,
                             weights = ate(des3, dichotomy = o > 2 ~ .)))

  des4 <- obs_design(o ~ uoa(cid1, cid2), data = simdata,
                     dichotomy = o > 2 ~ .)

  expect_silent(mod <- lmitt(y ~ 1, data = simdata, design = des3,
                              weights = ett(des4)))

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

  expect_error(vcov(damod1, type = "not_a_type"), "not defined")
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

  expect_error(confint(damod1, type = "not_a_type"), "not defined")

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

  vcovDA_z.95 <- matrix(rep(damod1$coef, 2), nrow = 2) +
                          sqrt(diag(vcovDA(damod1))) %o%
                          qt(c(0.025, 0.975), damod1$df.residual)
  ci1 <- confint(damod1)
  expect_true(all.equal(ci1, vcovDA_z.95, check.attributes = FALSE))

  vcovDA_z.9 <- matrix(rep(damod1$coef, 2), nrow = 2) +
                          sqrt(diag(vcovDA(damod1))) %o%
                          qt(c(0.05, 0.95), damod1$df.residual)
  ci1 <- confint(damod1, level = 0.9)
  expect_true(all.equal(ci1, vcovDA_z.9, check.attributes = FALSE))

  uoas <- apply(simdata[, c("cid1", "cid2")], 1,
                function(...) {
                  paste(..., collapse = "_")
                })
  vcovlm.9 <- damod2$coefficients +
    sqrt(diag(sandwich::sandwich(damod2, meat. = sandwich::meatCL,
                                 cluster = uoas))) %o%
    qt(c(0.05, 0.95), damod2$df.residual)
  ci1 <- confint(damod2, level = 0.9)
  expect_true(all.equal(ci1, vcovlm.9, check.attributes = FALSE))

  vcovlm_z.95 <- damod2$coefficients +
    sqrt(diag(sandwich::sandwich(damod2, meat. = sandwich::meatCL,
                            cluster = uoas))) %o%
    qt(c(0.025, 0.975), damod2$df.residual)
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

test_that("@moderator slot", {
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

  expect_equal(noblocks_lmitt_fittedsbgrp@moderator, "force")
  expect_equal(blocked_lmitt_fittedsbgrp@moderator, "force")
  expect_equal(lmitt_fitted_nosbgrp@moderator, vector("character"))
  expect_equal(blocked_lmitt_fitted_nosbgrp@moderator, vector("character"))
  expect_equal(not_lmitt_fitted@moderator, vector("character"))
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

  nolmittmod1 <- lm(y ~ z, simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des)
  nolmittmod2 <- lm(y ~ z, data = simdata, offset = ca)
  mod2 <- lmitt(y ~ 1, data = simdata, design = des, offset = stats::predict(cmod))

  ef_expected1 <- estfun(nolmittmod1)
  ef_expected2 <- estfun(nolmittmod2)

  expect_true(all.equal(estfun(mod1), ef_expected1, check.attributes = FALSE))
  expect_true(all.equal(estfun(mod2), ef_expected2, check.attributes = FALSE))
})

test_that(paste("estfun.DirectAdjusted returns correct dimensions and alignment",
                "when exact alignment between C and Q is possible"), {
  set.seed(438)
  data(simdata)
  simdata$uid <- rownames(simdata)

  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod, by = "uid"))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, design = des, offset = cov_adj(cmod, by = "uid"))

  expect_equal(dim(estfun(mod1)), c(nrow(simdata), 2))
  expect_equal(estfun(mod1), estfun(mod2))
})

test_that(paste("estfun.DirectAdjusted returns correct dimensions when only",
                "inexact alignment between C and Q is possible"), {
  set.seed(438)
  data(simdata)

  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, design = des, offset = cov_adj(cmod))

  expect_equal(dim(estfun(mod1)), c(nrow(simdata), 2))
  expect_equal(dim(estfun(mod2)), c(nrow(simdata), 2))
})

test_that(paste("estfun.DirectAdjusted returns correct dimensions and alignment",
                "when exact alignment between C and Q is possible but they only",
                "partially overlap"), {
  set.seed(438)
  data(simdata)
  simdata$uid <- rownames(simdata)

  Q_data <- simdata[simdata$uid %in% seq_len(20),]
  shuffled_Q_data <- simdata[sample(rownames(Q_data)),]
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, design = des, offset = cov_adj(cmod, by = "uid"))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, design = des, offset = cov_adj(cmod, by = "uid"))

  expect_equal(dim(estfun(mod1)), c(nrow(simdata), 2))
  expect_equal(estfun(mod1), estfun(mod2))
})

test_that(paste("estfun.DirectAdjusted returns correct dimensions for partial",
                "overlap of C and Q when only inexact alignment is possible"), {
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
               c(nrow(simdata) + nrow(cmod_data[is.na(cmod_data$cid1),]), 2))
})

test_that(paste("estfun.DirectAdjusted returns correct dimensions for no",
                "overlap of C and Q when only inexact alignment is possible"), {
  data(simdata)
  set.seed(401)
  cmod_data <- data.frame(
    "x" = rnorm(100), "y" = rnorm(100), "cid1" = NA, "cid2" = NA
  )

  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))

  expect_equal(dim(estfun(mod)), c(nrow(simdata) + nrow(cmod_data), 2))
})

test_that("estfun.DirectAdjusted returns correct dimensions for rectangular A11_inv", {
  data(simdata)
  cmod <- robustbase::lmrob(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  dmod <- lmitt(lm(y ~ z.(des), simdata, offset = cov_adj(cmod)), design = des)

  out <- estfun(dmod)
  expect_equal(dim(out), c(50, 2))
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

test_that(paste(".align_and_extend_estfuns fails if not a DirectAdjusted object",
                "with a SandwichLayer offset"), {
  data(simdata)

  mod1 <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  mod2 <- lmitt(y ~ 1, design = des, data = simdata)
  expect_error(.align_and_extend_estfuns(mod1), "must be a fitted")
  expect_error(.align_and_extend_estfuns(mod2), "must be a fitted")
})

test_that(paste(".align_and_extend_estfuns when exact alignment of C and Q is",
                "possible and the samples fully overlap"), {
  set.seed(438)
  data(simdata)

  simdata$obs_id <- seq_len(nrow(simdata))
  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  cmod1 <- lm(y ~ x, simdata)
  cmod2 <- lm(y ~ x, shuffled_simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod1, by = "obs_id"))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, design = des, offset = cov_adj(cmod2, by = "obs_id"))

  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_equal(ef1$phi, ef1$phi[sort(simdata$obs_id, index.return = TRUE)$ix,])
  expect_equal(ef1$phi, ef2$phi)
  expect_equal(ef1$psi, ef1$psi[sort(simdata$obs_id, index.return = TRUE)$ix,, drop = FALSE])
  expect_equal(ef1$psi, ef2$psi)
})

test_that(paste(".align_and_extend_estfuns when exact alignment of C and Q is",
                "possible and Q is a subset of C"), {
  set.seed(438)
  data(simdata)

  simdata$obs_id <- seq_len(nrow(simdata))
  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  Q_data <- simdata[1:20,]
  shuffled_Q_data <- Q_data[sample(rownames(Q_data)),]
  cmod1 <- lm(y ~ x, simdata)
  cmod2 <- lm(y ~ x, shuffled_simdata)
  des1 <- rct_design(z ~ cluster(cid1, cid2), data = Q_data)
  des2 <- rct_design(z ~ cluster(cid1, cid2), data = shuffled_Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, design = des1, offset = cov_adj(cmod1, by = "obs_id"))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, design = des2, offset = cov_adj(cmod2, by = "obs_id"))

  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_equal(ef1$phi, ef1$phi[sort(simdata$obs_id, index.return = TRUE)$ix,])
  expect_equal(ef1$phi, ef2$phi)
  expect_equal(ef1$psi, ef1$psi[sort(simdata$obs_id, index.return = TRUE)$ix,, drop = FALSE])
  expect_equal(ef1$psi[21:50], rep(0, 30))
  expect_equal(ef1$psi, ef2$psi)
})

test_that(paste(".align_and_extend_estfuns when exact alignment of C and Q is",
                "possible and C is a subset of Q"), {
  set.seed(438)
  data(simdata)

  simdata$obs_id <- seq_len(nrow(simdata))
  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  C_data <- simdata[1:20,]
  shuffled_C_data <- C_data[sample(rownames(C_data)),]
  cmod1 <- lm(y ~ x, C_data)
  cmod2 <- lm(y ~ x, shuffled_C_data)
  des1 <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  des2 <- rct_design(z ~ cluster(cid1, cid2), data = shuffled_simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des1, offset = cov_adj(cmod1, by = "obs_id"))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, design = des2, offset = cov_adj(cmod2, by = "obs_id"))

  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_equal(ef1$phi, ef1$phi[sort(simdata$obs_id, index.return = TRUE)$ix,])
  expect_true(all(ef1$phi[30:50] == 0))
  expect_equal(ef1$phi, ef2$phi)
  expect_equal(ef1$psi, ef1$psi[sort(simdata$obs_id, index.return = TRUE)$ix,, drop = FALSE])
  expect_equal(ef1$psi, ef2$psi)
})

test_that(paste(".align_and_extend_estfuns when exact alignment of C and Q is",
                "possible and C and Q have no overlap"), {
  set.seed(438)
  data(simdata)

  simdata$obs_id <- seq_len(nrow(simdata))
  Q_data <- simdata[21:50,]
  C_data <- simdata[1:20,]
  shuffled_Q_data <- Q_data[sample(rownames(Q_data)),]
  shuffled_C_data <- C_data[sample(rownames(C_data)),]
  cmod1 <- lm(y ~ x, C_data)
  cmod2 <- lm(y ~ x, shuffled_C_data)
  des1 <- rct_design(z ~ cluster(cid1, cid2), data = Q_data)
  des2 <- rct_design(z ~ cluster(cid1, cid2), data = shuffled_Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, design = des1, offset = cov_adj(cmod1, by = "obs_id"))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, design = des2, offset = cov_adj(cmod2, by = "obs_id"))

  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_true(all(ef1$phi[21:50,] == 0))
  expect_equal(ef1$phi, ef2$phi)
  expect_equal(ef1$psi, ef1$psi[sort(simdata$obs_id, index.return = TRUE)$ix,, drop = FALSE])
  expect_true(all(ef1$psi[1:20,] == 0))
  expect_equal(ef1$psi, ef2$psi)
})

test_that(paste(".align_and_extend_estfuns when exact alignment of C and Q isn't",
                "possible and the samples fully overlap"), {
  set.seed(438)
  data(simdata)

  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  cmod1 <- lm(y ~ x, simdata)
  cmod2 <- lm(y ~ x, shuffled_simdata)
  des1 <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  des2 <- rct_design(z ~ cluster(cid1, cid2), data = shuffled_simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des1, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, design = des2, offset = cov_adj(cmod2))

  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)
  by_ix <- sort(apply(simdata[, c("cid1", "cid2")], 1,
                      function(...) paste(..., collapse = "_")))

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_equal(dim(ef2$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef2$psi), c(nrow(simdata), 2))

  phi1_sorted <- lapply(split(ef1$phi, by_ix), sort)
  phi2_sorted <- lapply(split(ef2$phi, by_ix), sort)
  expect_true(all(sapply(unique(by_ix),
                         function(id) all.equal(phi1_sorted[[id]], phi2_sorted[[id]]))))

  psi1_sorted <- lapply(split(ef1$psi, by_ix), sort)
  psi2_sorted <- lapply(split(ef2$psi, by_ix), sort)
  expect_true(all(sapply(unique(by_ix),
                         function(id) all.equal(psi1_sorted[[id]], psi2_sorted[[id]]))))
})

test_that(paste(".align_and_extend_estfuns when exact alignment of C and Q isn't",
                "possible and Q is a subset of C"), {
  set.seed(438)
  data(simdata)

  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  Q_data <- simdata[1:20,]
  shuffled_Q_data <- Q_data[sample(rownames(Q_data)),]
  cmod1 <- lm(y ~ x, simdata)
  cmod2 <- lm(y ~ x, shuffled_simdata)
  des1 <- rct_design(z ~ cluster(cid1, cid2), data = Q_data)
  des2 <- rct_design(z ~ cluster(cid1, cid2), data = shuffled_Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, design = des1, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, design = des2, offset = cov_adj(cmod2))

  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)
  by_ix <- sort(apply(simdata[, c("cid1", "cid2")], 1,
                      function(...) paste(..., collapse = "_")))

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_equal(dim(ef2$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef2$psi), c(nrow(simdata), 2))

  phi1_sorted <- lapply(split(ef1$phi, by_ix), sort)
  phi2_sorted <- lapply(split(ef2$phi, by_ix), sort)
  expect_true(all(sapply(unique(by_ix),
                         function(id) all.equal(phi1_sorted[[id]], phi2_sorted[[id]]))))

  psi1_sorted <- lapply(split(ef1$psi, by_ix), sort)
  psi2_sorted <- lapply(split(ef2$psi, by_ix), sort)
  expect_true(all(sapply(unique(by_ix),
                         function(id) all.equal(psi1_sorted[[id]], psi2_sorted[[id]]))))
  expect_true(all(ef1$psi[21:nrow(simdata),] == 0))
  expect_true(all(ef2$psi[21:nrow(simdata),] == 0))
})

test_that(paste(".align_and_extend_estfuns when exact alignment of C and Q isn't",
                "possible and C is a subset of Q"), {
  set.seed(438)
  data(simdata)

  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  C_data <- simdata[1:20,]
  shuffled_C_data <- C_data[sample(rownames(C_data)),]
  cmod1 <- lm(y ~ x, C_data)
  cmod2 <- lm(y ~ x, shuffled_C_data)
  des1 <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  des2 <- rct_design(z ~ cluster(cid1, cid2), data = shuffled_simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des1, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, design = des2, offset = cov_adj(cmod2))

  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)
  by_ix <- sort(apply(simdata[, c("cid1", "cid2")], 1,
                      function(...) paste(..., collapse = "_")))

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_equal(dim(ef2$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef2$psi), c(nrow(simdata), 2))

  phi1_sorted <- lapply(split(ef1$phi, by_ix), sort)
  phi2_sorted <- lapply(split(ef2$phi, by_ix), sort)
  expect_true(all(sapply(unique(by_ix),
                         function(id) all.equal(phi1_sorted[[id]], phi2_sorted[[id]]))))
  expect_true(all(ef1$phi[21:nrow(simdata),] == 0))
  expect_true(all(ef2$phi[21:nrow(simdata),] == 0))

  psi1_sorted <- lapply(split(ef1$psi, by_ix), sort)
  psi2_sorted <- lapply(split(ef2$psi, by_ix), sort)
  expect_true(all(sapply(unique(by_ix),
                         function(id) all.equal(psi1_sorted[[id]], psi2_sorted[[id]]))))
})

test_that(paste(".align_and_extend_estfuns when exact alignment of C and Q isn't",
                "possible and the samples have no overlap"), {
  set.seed(438)
  data(simdata)

  Q_data <- simdata[21:50,]
  C_data <- simdata[1:20,]
  shuffled_Q_data <- Q_data[sample(rownames(Q_data)),]
  shuffled_C_data <- C_data[sample(rownames(C_data)),]
  cmod1 <- lm(y ~ x, C_data)
  cmod2 <- lm(y ~ x, shuffled_C_data)
  des1 <- rct_design(z ~ cluster(cid1, cid2), data = Q_data)
  des2 <- rct_design(z ~ cluster(cid1, cid2), data = shuffled_Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, design = des1, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, design = des2, offset = cov_adj(cmod2))

  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)
  by_ix <- sort(apply(simdata[, c("cid1", "cid2")], 1,
                      function(...) paste(..., collapse = "_")))

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_equal(dim(ef2$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef2$psi), c(nrow(simdata), 2))

  phi1_sorted <- lapply(split(ef1$phi, by_ix), sort)
  phi2_sorted <- lapply(split(ef2$phi, by_ix), sort)
  expect_true(all(sapply(unique(by_ix),
                         function(id) all.equal(phi1_sorted[[id]], phi2_sorted[[id]]))))
  expect_true(all(ef1$phi[21:50,] == 0))
  expect_true(all(ef2$phi[21:50,] == 0))

  psi1_sorted <- lapply(split(ef1$psi, by_ix), sort)
  psi2_sorted <- lapply(split(ef2$psi, by_ix), sort)
  expect_true(all(sapply(unique(by_ix),
                         function(id) all.equal(psi1_sorted[[id]], psi2_sorted[[id]]))))
  expect_true(all(ef1$psi[1:20,] == 0))
  expect_true(all(ef2$psi[1:20,] == 0))
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

  expect_warning(ids1 <- .make_uoa_ids(dmod1), "treated as independent")
  # expect_warning(ids2 <- .make_uoa_ids(dmod2), "NA's for some but not all")
  expect_warning(ids2 <- .make_uoa_ids(dmod2), "treated as independent")

  expect_true(is.factor(ids1))
  expect_true(is.factor(ids2))

  expect_equal(length(ids1), nrow(simdata) + 10)
  expect_equal(length(ids2), nrow(simdata) + 10)

  expect_true(all.equal(ids1[1:nrow(simdata)], factor(Q_uoas), check.attributes = FALSE))
  expect_true(all.equal(ids2[1:nrow(simdata)], factor(Q_uoas), check.attributes = FALSE))

  expect_equal(length(unique(ids1)), length(unique(Q_uoas)) + 10)
  expect_equal(length(unique(ids2)), length(unique(Q_uoas)) + 10)
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

  expect_warning(ids <- .make_uoa_ids(dmod), "treated as independent")

  expect_true(is.factor(ids))

  expect_equal(length(ids), nrow(simdata) + 20)

  expect_true(all.equal(ids[1:nrow(simdata)], factor(Q_uoas), check.attributes = FALSE))

  expect_equal(length(unique(ids)), length(unique(Q_uoas)) + 20)
})

test_that(paste(".make_uoa_ids returns correct ID's when cov_adj's 'by' argument",
                "provides a different ordering"), {
  data(simdata)
  simdata$id <- sample(seq_len(nrow(simdata)))
  
  set.seed(300)
  C_not_Q <- data.frame("y" = rnorm(20), "x" = rnorm(20), "cid1" = NA, "cid2" = NA,
                        "id" = seq(51, 70))
  cmod_data <- rbind(simdata[, colnames(C_not_Q)], C_not_Q)
  
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod, by = "id"))
  
  Q_uoas <- apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1,
                  function(...) paste(..., collapse = "_"))
  
  ids <- .make_uoa_ids(dmod)
  
  expect_true(is.factor(ids))
  
  expect_equal(length(ids), nrow(simdata) + 20)
  
  expect_true(all.equal(ids[1:nrow(simdata)],
                        factor(Q_uoas)[sort(simdata$id, index.return = TRUE)$ix],
                        check.attributes = FALSE))
  
  expect_equal(length(unique(ids)), length(unique(Q_uoas)) + 20)
})

test_that(paste(".order_samples fails without a DirectAdjusted object or",
                "SandwichLayer offset"), {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lm(y ~ z, simdata, offset = predict(cmod))
  mod2 <- lmitt(y ~ 1, data = simdata, design = des)

  expect_error(.order_samples(mod1), "must be a DirectAdjusted object")
  expect_error(.order_samples(mod2), "must be a DirectAdjusted object")
})

test_that(".order_samples when the samples fully overlap", {
  set.seed(300)
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod, by = "uid"))

  out <- .order_samples(mod)

  expect_equal(length(out$Q_order), nrow(simdata))
  expect_equal(length(out$C_order), nrow(simdata))
  expect_equal(length(out$Q_union_C_order), nrow(simdata))

  expect_equal(names(out$Q_order), as.character(seq_len(nrow(simdata))))
  expect_equal(names(out$C_order), as.character(seq_len(nrow(simdata))))
  expect_equal(names(out$Q_union_C_order), as.character(seq_len(nrow(simdata))))
  
  expect_true(all.equal(out$Q_order, seq_len(nrow(simdata)), check.attributes = FALSE))
  expect_true(all.equal(out$C_order, seq_len(nrow(simdata)), check.attributes = FALSE))
  expect_true(all.equal(out$Q_union_C_order, seq_len(nrow(simdata)), check.attributes = FALSE))
})

test_that(".order_samples when Q is a subset of C", {
  set.seed(300)
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  cmod_data <- rbind(
    simdata[, c("x", "y", "cid1", "cid2", "uid")],
    data.frame(x = rnorm(30), y = rnorm(30), cid1 = NA, cid2 = NA,
               uid = seq_len(30) + nrow(simdata))
  )
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod, by = "uid"))

  out <- .order_samples(mod)

  expect_equal(length(out$Q_order), nrow(simdata))
  expect_equal(length(out$C_order), nrow(simdata) + 30)
  expect_equal(length(out$Q_union_C_order), nrow(simdata) + 30)

  expect_equal(names(out$Q_order), as.character(seq_len(nrow(simdata))))
  expect_equal(names(out$C_order), as.character(seq_len(nrow(simdata) + 30)))
  expect_equal(names(out$Q_union_C_order), as.character(seq_len(nrow(simdata) + 30)))
  
  expect_true(all.equal(out$Q_order, seq_len(nrow(simdata)), check.attributes = FALSE))
  expect_true(all.equal(out$C_order, seq_len(nrow(simdata) + 30), check.attributes = FALSE))
  expect_true(all.equal(out$Q_union_C_order, seq_len(nrow(simdata) + 30), check.attributes = FALSE))
})


test_that(".order_samples when C is a subset of Q", {
  set.seed(300)
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  cmod_data <- simdata[1:20,]
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod, by = "uid"))

  out <- .order_samples(mod)

  expect_equal(length(out$Q_order), nrow(simdata))
  expect_equal(length(out$C_order), 20)
  expect_equal(length(out$Q_union_C_order), nrow(simdata))

  expect_equal(names(out$Q_order), as.character(seq_len(50)))
  expect_equal(names(out$C_order), as.character(seq_len(20)))
  expect_equal(names(out$Q_union_C_order), as.character(seq_len(50)))

  expect_true(all.equal(out$Q_order, seq_len(nrow(simdata)), check.attributes = FALSE))
  expect_true(all.equal(out$C_order, seq_len(20), check.attributes = FALSE))
  expect_true(all.equal(out$Q_union_C_order, seq_len(nrow(simdata)), check.attributes = FALSE))
})

test_that(".order_samples when the samples do not overlap", {
  set.seed(300)
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  cmod_data <- simdata[1:20,]
  Q_data <- simdata[21:50,]
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_data)
  mod <- lmitt(y ~ 1, data = Q_data, design = des, offset = cov_adj(cmod, by = "uid"))

  out <- .order_samples(mod)

  expect_equal(length(out$Q_order), 30)
  expect_equal(length(out$C_order), 20)
  expect_equal(length(out$Q_union_C_order), nrow(simdata))

  expect_equal(names(out$Q_order), as.character(seq(21, 50)))
  expect_equal(names(out$C_order), as.character(seq_len(20)))
  expect_equal(names(out$Q_union_C_order), as.character(seq_len(nrow(simdata))))

  expect_true(all.equal(out$Q_order, seq_len(30), check.attributes = FALSE))
  expect_true(all.equal(out$C_order, seq_len(20), check.attributes = FALSE))
  expect_true(all.equal(out$Q_union_C_order, c(seq(31, 50), seq_len(30)), check.attributes = FALSE))
})

test_that(".order_samples when no `by` argument provided", {
  set.seed(300)
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  C_data <- simdata[1:20,]
  Q_data <- simdata[21:50,]
  cmod <- lm(y ~ x, C_data)
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_data)
  mod <- lmitt(y ~ 1, data = Q_data, design = des, offset = cov_adj(cmod))

  out <- .order_samples(mod)

  Q_uoas <- sort(apply(Q_data[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_")))
  C_uoas <- sort(apply(C_data[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_")))
  names(Q_uoas) <- NULL
  names(C_uoas) <- NULL

  expect_equal(length(out$Q_order), 30)
  expect_equal(length(out$C_order), 20)
  expect_equal(length(out$Q_union_C_order), nrow(simdata))

  expect_equal(names(out$Q_order), Q_uoas)
  expect_equal(names(out$C_order), C_uoas)
  expect_equal(names(out$Q_union_C_order), sort(c(Q_uoas, C_uoas)))

  expect_true(all.equal(out$Q_order, seq_len(30), check.attributes = FALSE))
  expect_true(all.equal(out$C_order, seq_len(20), check.attributes = FALSE))
  expect_true(all.equal(out$Q_union_C_order, c(seq(31, 50), seq_len(30)), check.attributes = FALSE))
})

test_that("sanitize_Q_ids fails with invalid cluster argument", {
  data(simdata)
  mod <- lm(y ~ z, data = simdata)
  expect_error(.sanitize_Q_ids(mod, by = "not_uoas"),
               "columns not_uoas in ITT effect model data")

  invalid_ids <- apply(simdata[, c("cid1", "cid2")], 1,
                       function(...) paste(..., collapse = "_"))
  expect_error(.sanitize_Q_ids(mod, by = invalid_ids),
               ", 5_2 in ITT effect model data")
})

test_that("sanitize_Q_ids succeeds with valid `by` argument", {
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  cmod <- lm(y ~ z, data = simdata)
  des <- rct_design(z ~ unitid(cid1, cid2, uid), simdata)
  dmod <- lmitt(y ~ 1, design = des, data = simdata, offset = cov_adj(cmod))
  out <- .sanitize_Q_ids(dmod, by = c("cid1", "cid2"))

  expected_out <- apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1,
                        function(...) paste(..., collapse = "_"))
  expect_equal(out, expected_out)
})

test_that(".base_S3class_estfun fails with invalid base S3 class", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des)
  mod@.S3Class <- "invalid_class"

  expect_error(.base_S3class_estfun(mod), "must have been fitted")
})

test_that("checking proper errors in conversion from lm to DA", {
  data(simdata)
  des <- rct_design(z ~ unitid(cid1, cid2), simdata)

  expect_error(as.lmitt(lm(y ~ x, data = simdata), design = des),
               "aliases are not found")


  # Take in either Design or WeightedDesign
  mod1 <- propertee:::.convert_to_lmitt(lm(y ~ assigned(des), data = simdata),
                                     des,
                                     FALSE,
                                     TRUE,
                                     "a",
                                     call("quote", call("ls")))
  mod2 <- propertee:::.convert_to_lmitt(lm(y ~ assigned(des), data = simdata),
                                     ate(des, data = simdata),
                                     FALSE,
                                     TRUE,
                                     "a",
                                     call("quote", call("ls")))
  expect_identical(mod1@Design, mod2@Design)

  expect_error(propertee:::.convert_to_lmitt(lm(y ~ assigned(des), data = simdata),
                                     1,
                                     FALSE,
                                     TRUE,
                                     "a",
                                     call("quote", call("ls"))), "must be a")


})


test_that("disable update.DA", {
  data(simdata)
  des <- rct_design(z ~ unitid(cid1, cid2), simdata)

  mod <- lmitt(y ~ 1, data = simdata, design = des)

  expect_error(update(mod),
               "DirectAdjusted objects do not support")
})

test_that("lmitt_call", {
  data(simdata)
  des <- rct_design(z ~ unitid(cid1, cid2), simdata)

  # Make sure slot is actual call
  mod <- lmitt(y ~ 1, data = simdata, design = des)
  expect_true(is.call(mod@lmitt_call))

  # Call via `lmitt()` should match
  lmittcall <- str2lang("lmitt(y ~ 1, data = simdata, design = des)")
  mod <- eval(lmittcall)
  expect_equal(lmittcall, mod@lmitt_call)

  # Call via `lmitt.formula()` should match
  lmittformcall <- str2lang("lmitt.formula(y ~ 1, data = simdata, design = des)")
  mod <- eval(lmittformcall)
  expect_equal(lmittformcall, mod@lmitt_call)

  # as.lmitt, make sure is actual call
  lmmod <- lm(y ~ a.(des), data = simdata)
  mod <- as.lmitt(lmmod, design = des)
  expect_true(is.call(mod@lmitt_call))

  # Call via `as.lmitt()` should match
  aslmittcall <- str2lang("as.lmitt(lmmod, design = des)")
  mod <- eval(aslmittcall)
  expect_equal(aslmittcall, mod@lmitt_call)

  # Call via `lmitt(lm` should match
  lmittlmcall <- str2lang("lmitt(lmmod, design = des)")
  mod <- eval(lmittlmcall)
  expect_equal(lmittlmcall, mod@lmitt_call)

  # Call via `lmitt.lm` should match
  lmittlmcall <- str2lang("lmitt.lm(lmmod, design = des)")
  mod <- eval(lmittlmcall)
  expect_equal(lmittlmcall, mod@lmitt_call)

})

test_that("lmitt_fitted object", {
  data(simdata)
  des <- rct_design(z ~ unitid(cid1, cid2), simdata)

  mod <- lmitt(y ~ as.factor(o), data = simdata, design = des)
  expect_true(mod@lmitt_fitted)

  mod <- as.lmitt(lm(y ~ adopters(des), data = simdata), design = des)
  expect_false(mod@lmitt_fitted)

  mod <- lmitt(lm(y ~ adopters(des), data = simdata), design = des)
  expect_false(mod@lmitt_fitted)

})

test_that("printed effects aren't confused by bad naming", {
  data(simdata)
  simdata$abz.c <- as.factor(simdata$o)
  des <- rct_design(z ~ unitid(cid1, cid2), simdata)

  mod <- lmitt(y ~ abz.c, data = simdata, design = des)
  co <- capture.output(show(mod))
  cos <- strsplit(trimws(co), "[[:space:]]+")

  expect_true(all(vapply(cos, length, numeric(1)) == 4))
  expect_true(all(!isTRUE(grepl("^abz", cos[[1]]))))
  expect_true(all(grepl("^z", cos[[1]])))

  # to force ` in variable names via as.factor
  mod <- lmitt(y ~ as.factor(abz.c), data = simdata, design = des)
  co <- capture.output(show(mod))
  cos <- strsplit(trimws(co[c(1, 3)]), "[[:space:]]+")
  cos <- Reduce("c", cos)
  expect_length(cos, 4)
  expect_true(all(!isTRUE(grepl("^as\\.factor", cos))))
  expect_true(all(grepl("^`z", cos)))
})


test_that("as.lmitt call inside mapply", {
  min_df <- data.frame(schoolid = seq_len(20),
                       response_col = rnorm(20),
                       grade = rep(3:4, each = 10),
                       adsy = rep(rep(c(0, 1), each = 5), 2))

  min_des <- rct_design(adsy ~ unitid(schoolid), min_df)

  res <- mapply(function(fmla) {
    as.lmitt(lm(fmla, min_df, weights = ett(min_des)))
  },
  list(response_col ~ z.(), response_col ~ grade + grade:z.()),
  SIMPLIFY = FALSE)


  expect_true(all(sapply(res, is, "DirectAdjusted")))
  expect_true(all(sapply(sapply(sapply(res,
                                       slot, "lmitt_call"),
                                "[[", 1),
                         deparse) == "as.lmitt"))

  res <- lapply(list(response_col ~ z.(), response_col ~ grade + grade:z.()),
                function(fmla) {
                  as.lmitt(lm(fmla, min_df, weights = ett(min_des)))
                })

  expect_true(all(sapply(res, is, "DirectAdjusted")))
  expect_true(all(sapply(sapply(sapply(res,
                                       slot, "lmitt_call"),
                                "[[", 1),
                         deparse) == "as.lmitt"))

})

test_that("#81 continuous moderator shows appropriate coefficients", {
  data(simdata)
  des <- obs_design(z ~ uoa(cid1, cid2) , data = simdata)

  mod <- lmitt(y ~ x, data = simdata, design = des)
  expect_length(mod$coeff, 4)
  coefnames <- strsplit(capture.output(show(mod))[1], " +")[[1]]
  coefnames <- coefnames[nchar(coefnames) > 0]
  expect_true(all(grepl("z.", coefnames, fixed = TRUE)))
})
