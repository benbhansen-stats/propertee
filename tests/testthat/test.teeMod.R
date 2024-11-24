test_that(paste("teeMod object created correctly with weights and no",
                "SandwichLayer in the lm call"), {

  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa2, uoa1) + block(bid), data = simdata)

  sslm <- new("teeMod",
              lm(y ~ assigned(), data = simdata, weights = ate(spec)),
              StudySpecification = spec,
              lmitt_fitted = TRUE)

  expect_s4_class(sslm, "teeMod")
  expect_true(inherits(sslm, "lm"))

  expect_identical(sslm$model$"(weights)"@StudySpecification, spec)
  expect_identical(sslm$model$"(weights)"@StudySpecification, sslm@StudySpecification)
})

test_that(paste("teeMod object created correctly with weights and ",
                "SandwichLayer in the lm call"), {
  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa2, uoa1) + block(bid), data = simdata)
  cmod <- lm(y ~ x, data = simdata)
  sslm <- new("teeMod",
              lm(y ~ assigned(), data = simdata, weights = ate(spec),
                 offset = cov_adj(cmod)),
              StudySpecification = spec,
              lmitt_fitted = FALSE)

  expect_s4_class(sslm, "teeMod")
  expect_true(inherits(sslm, "lm"))

  expect_equal(sslm$model$`(offset)`@.Data, as.numeric(cmod$fitted.values))
  expect_identical(sslm$model$`(offset)`@fitted_covariance_model, cmod)
  expect_equal(sslm$model$`(offset)`@prediction_gradient,
               stats::model.matrix(cmod))
  expect_identical(sslm$model$`(offset)`@keys,
                   cbind(simdata[,var_names(spec, "u")], in_Q = rep(TRUE, nrow(simdata))))
})

test_that("teeMod print/show", {

  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa2, uoa1) + block(bid), data = simdata)
  cmod <- lm(y ~ z, data = simdata)

  ######## Pretend it comes from `as.lmitt`
  sslm <- new("teeMod",
              lm(y ~ assigned(), data = simdata, weights = ate(spec),
                 offset = cov_adj(cmod)),
              StudySpecification = spec,
              lmitt_fitted = FALSE)

  aslm <- as(sslm, "lm")

  expect_silent(invisible(capture.output(expect_identical(print(sslm), sslm))))
  expect_silent(invisible(capture.output(expect_identical(show(sslm), sslm))))

  # Expect "assigned()"
  expect_output(print(sslm), "assigned()")
  expect_output(show(sslm), "assigned()")

  ######## Now pretend it comes from `lmitt.formula`
  # Need `txt_` to match what comes out of a call like `lmitt(y ~ 1...)`
  simdata$z. <- simdata$z
  sslm <- new("teeMod",
              lm(y ~ z., data = simdata, weights = ate(spec),
                 offset = cov_adj(cmod)),
              StudySpecification = spec,
              lmitt_fitted = TRUE)

  aslm <- as(sslm, "lm")

  expect_silent(invisible(capture.output(expect_identical(print(sslm), sslm))))
  expect_silent(invisible(capture.output(expect_identical(show(sslm), sslm))))

  # Expect "assigned()"
  expect_output(print(sslm), "z\\.")
  expect_output(show(sslm), "z\\.")
})

test_that("lm to teeMod succeeds with weights and no SandwichLayer", {

  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  mod <- lm(y ~ assigned(), data = simdata, weights = ate(spec))

  mod_ss <- as.lmitt(mod)

  expect_s4_class(mod_ss, "teeMod")
  expect_true(inherits(mod_ss, "lm"))

  expect_identical(mod_ss$model$"(weights)"@StudySpecification, spec)

  mod_lmitt <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)

  expect_true(all(round(mod_lmitt$coef, 4) %in% round(mod_ss$coef, 4)))
  expect_identical(mod_ss@StudySpecification, mod_lmitt@StudySpecification)
})

test_that("lm to teeMod with weights and SandwichLayer", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  cmod <- lm(y ~ x, data = simdata)

  mod <- lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))

  mod_ss <- as.lmitt(mod)

  expect_s4_class(mod_ss, "teeMod")
  expect_true(inherits(mod_ss, "lm"))

  expect_equal(mod_ss$model$`(offset)`@.Data, as.numeric(cmod$fitted.values))
  expect_identical(mod_ss$model$`(offset)`@fitted_covariance_model, cmod)
  expect_equal(mod_ss$model$`(offset)`@prediction_gradient,
               stats::model.matrix(cmod))
  expect_identical(mod_ss$model$`(offset)`@keys,
                   cbind(simdata[,var_names(spec, "u")], in_Q = rep(TRUE, nrow(simdata))))

  mod_lmitt <- lmitt(y ~ 1, data = simdata, weights = ate(),
                     offset = cov_adj(cmod), specification = spec)

  expect_true(all(round(mod_lmitt$coef, 4) %in%  round(mod_ss$coef, 4)))
  expect_identical(mod_ss@StudySpecification, mod_lmitt@StudySpecification)
})

test_that("Conversion from lm to teeMod fails without an lm object", {
  expect_error(as.lmitt(1),
               "lm object")
})

test_that("lm to teeMod fails without a StudySpecification object", {
  data(simdata)

  expect_error(as.lmitt(lm(y ~ assigned(), data = simdata,
                                    weights = seq_len(nrow(simdata)))),
               "Unable to locate StudySpecification")
})

test_that("vcov, confint, etc", {
  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa2, uoa1) + block(bid), data = simdata)

  sslm <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec)))

  expect_true(is.matrix(vcov(sslm)))
  expect_equal(dim(vcov(sslm)), c(2, 2))

  expect_true(is.matrix(confint(sslm)))
  expect_equal(dim(confint(sslm)), c(2, 2))


})

test_that("subsetting model with weights and/or cov_adj", {
  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  mod1 <- lm(y ~ assigned(),
             data = simdata,
             weights = ate(spec),
             offset = cov_adj(lm(y ~ x, data = simdata)))

  expect_true(inherits(mod1$model$`(weights)`, "WeightedStudySpecification"))
  expect_true(inherits(mod1$model$`(offset)`, "SandwichLayer"))

  # add subsetting
  mod2 <- lm(y ~ assigned(),
             data = simdata,
             weights = ate(spec),
             offset = cov_adj(lm(y ~ x, data = simdata)),
             subset = simdata$dose < 300)

  expect_true(inherits(mod2$model$`(weights)`, "WeightedStudySpecification"))
  expect_true(inherits(mod2$model$`(offset)`, "SandwichLayer"))

})

test_that("Differing specifications found", {

  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  spec2 <- obs_spec(o ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  # This should error, since there's no `specification` argument to `lm`, `assigned()`
  # finds both the `spec` inside `ate()` and `spec2` inside `cov_adj()`
  ## expect_error(lm(y ~ assigned(),
  ##                 data = simdata,
  ##                 weights = ate(spec),
  ##                 offset = cov_adj(lm(y ~ x, data = simdata),
  ##                                  specification = spec2)),
  ##              "Multiple differing")

  expect_error(lmitt(y ~ 1,
                        data = simdata,
                        specification = spec,
                        weights = ate(),
                        offset = cov_adj(lm(y ~ x, data = simdata),
                                         specification = spec2)),
               "Multiple differing")

  expect_error(lmitt(y ~ 1,
                        data = simdata,
                        specification = spec,
                        weights = ate(),
                        offset = cov_adj(lm(y ~ x, data = simdata),
                                         specification = spec2)),
               "Multiple differing")


  mod <- lm(y ~ assigned(), data = simdata,
            offset = cov_adj(lm(y ~ x, data = simdata),
                             specification = spec2))

  expect_error(as.lmitt(mod, specification = spec), "Multiple differing")

  expect_error(lmitt(y ~ 1, data = simdata,
                     offset = cov_adj(lm(y ~ x, data = simdata),
                                      specification = spec2),
                     specification = spec),
               "Multiple differing")

  # Checking that things work when passing multiple of the same specifications
  mod1 <- lmitt(y ~ x, data = simdata, weights = ate(), specification = spec)
  mod2 <- lmitt(y ~ x, data = simdata, weights = ate(spec), specification = spec)
  expect_equal(mod1$coefficients, mod2$coefficients)

  mod3 <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)
  mod4 <- lmitt(y ~ 1, data = simdata, weights = ate(spec), specification = spec)
  expect_equal(mod3$coefficients, mod4$coefficients)


  spec3 <- obs_spec(o ~ uoa(uoa1, uoa2), dat = simdata)

  # No error with a dichotomy
  expect_silent(mod <- lmitt(y ~ 1, data = simdata, specification = spec3,
                             weights = ate(spec3, dichotomy = o > 2 ~ .)))

})

test_that("teeMod object has its own evaluation environment", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec)
  mod2 <- lmitt(lm(y ~ assigned(), simdata, weights = ate(spec)))

  expect_false(identical(environment(), environment(formula(mod1))))
  expect_false(identical(environment(), environment(formula(mod2))))

  # `lmitt.formula` builds a data that incluspec the passed-in data, plus updated
  # RHS and LHS for centering and absorb.
  expect_true(all(simdata %in% environment(formula(mod1))$data))
  expect_true(all(environment(formula(mod2))$data %in% environment(formula(mod1))$data))
  expect_equal(environment(formula(mod1))$specification, spec)
  expect_equal(environment(formula(mod1))$specification, environment(formula(mod2))$specification)
})

test_that("vcov.teeMod handles vcov_tee arguments and non-SL offsets", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  ssmod1 <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec),
                     offset = cov_adj(cmod)))
  ssmod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)

  vmat1 <- vcov(ssmod1)
  vmat2 <- vcov(ssmod1, type = "CR0")
  vmat3 <- vcov(ssmod1, cadjust = FALSE)

  expect_error(vcov(ssmod1, type = "not_a_type"), "not defined")
  expect_identical(vmat1, vmat2)
  expect_identical(vmat1, vcov_tee(ssmod1))
  expect_identical(vmat3, vcov_tee(ssmod1, cadjust = FALSE))

  uoas <- apply(simdata[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))
  vmat3 <- vcov(ssmod2)
  expect_true(all.equal(
    vmat3,
    sandwich::sandwich(ssmod2, meat. = sandwich::meatCL, cluster = uoas),
    check.attributes = FALSE))
})

test_that("confint.teeMod handles vcov_tee `type` arguments and non-SL offsets", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  ssmod1 <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec,
                     offset = cov_adj(cmod))
  ssmod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)

  expect_error(confint(ssmod1, type = "not_a_type"), "not defined")

  # default CI
  vcov_tee_ci.95 <- ssmod1$coefficients + sqrt(diag(vcov_tee(ssmod1))) %o%
    qt(c(0.025, 0.975), ssmod1$df.residual)
  dimnames(vcov_tee_ci.95) <- list(names(ssmod1$coefficients), c("2.5 %", "97.5 %"))
  ci1 <- confint(ssmod1, type = "CR0")
  ci2 <- confint(ssmod1)
  expect_equal(ci1, ci2)
  expect_equal(ci1, vcov_tee_ci.95)

  # HC1 CI
  vcov_tee_HC1_ci.95 <- ssmod1$coefficients + sqrt(diag(vcov_tee(ssmod1, type = "HC1"))) %o%
    qt(c(0.025, 0.975), ssmod1$df.residual)
  dimnames(vcov_tee_HC1_ci.95) <- list(names(ssmod1$coefficients), c("2.5 %", "97.5 %"))
  ci_HC1 <- confint(ssmod1, type = "HC1")
  expect_equal(ci_HC1, vcov_tee_HC1_ci.95)

  # CI with different level
  vcov_tee_ci.9 <- ssmod1$coefficients + sqrt(diag(vcov_tee(ssmod1))) %o%
    qt(c(0.05, 0.95), ssmod1$df.residual)
  dimnames(vcov_tee_ci.9) <- list(names(ssmod1$coefficients), c("5 %", "95 %"))
  ci1 <- confint(ssmod1, level = 0.9)
  expect_equal(ci1, vcov_tee_ci.9)

  # CI with lmitt.lm
  uoas <- apply(simdata[, c("uoa1", "uoa2")], 1,
                function(...) {
                  paste(..., collapse = "_")
                })

  vcovlm_z.95 <- ssmod2$coefficients +
    sqrt(diag(sandwich::sandwich(ssmod2, meat. = sandwich::meatCL,
                            cluster = uoas))) %o%
    qt(c(0.025, 0.975), ssmod2$df.residual)
  ci1 <- confint(ssmod2)
  expect_true(all.equal(ci1, vcovlm_z.95, check.attributes = FALSE))
})

test_that("absorbed_intercepts", {
  data(simdata)

  blockedspec <- rct_spec(z ~ block(bid) + cluster(uoa1, uoa2), data = simdata)
  noblocksspec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  blocked_lmitt_fitted_absorbed <- lmitt(y ~ 1, data = simdata,
                                         specification = blockedspec, absorb = TRUE)
  expect_message(blocked_lmitt_fitted_not_absorbed <-
                   lmitt(y ~ 1, data = simdata, specification = blockedspec,
                         absorb = FALSE))
  blocked_not_lmitt_fitted <- as.lmitt(lm(y ~ assigned(blockedspec), data = simdata),
                                       specification = blockedspec)

  noblocks_lmitt_fitted_absorbed <- lmitt(y ~ 1, data = simdata,
                                          specification = noblocksspec, absorb = TRUE)
  noblocks_lmitt_fitted_not_absorbed <- lmitt(y ~ 1, data = simdata,
                                              specification = noblocksspec, absorb = FALSE)

  expect_true(blocked_lmitt_fitted_absorbed@absorbed_intercepts)
  expect_false(blocked_lmitt_fitted_not_absorbed@absorbed_intercepts)
  expect_false(blocked_not_lmitt_fitted@absorbed_intercepts)
  expect_false(noblocks_lmitt_fitted_not_absorbed@absorbed_intercepts)
  expect_true(noblocks_lmitt_fitted_absorbed@absorbed_intercepts)
})

test_that("@moderator slot", {
  data(simdata)

  blockedspec <- rct_spec(z ~ block(bid) + cluster(uoa1, uoa2), data = simdata)
  noblocksspec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  noblocks_lmitt_fittedsbgrp <- lmitt(y ~ force, data = simdata,
                                      specification = noblocksspec)
  expect_message(blocked_lmitt_fittedsbgrp <- lmitt(y ~ force, data = simdata,
                                                    specification = blockedspec))
  lmitt_fitted_nosbgrp <- lmitt(y ~ 1, data = simdata, specification = noblocksspec)
  expect_message(blocked_lmitt_fitted_nosbgrp <- lmitt(y ~ 1, data = simdata,
                                                       specification = blockedspec))
  not_lmitt_fitted <- as.lmitt(lm(y ~ assigned(noblocksspec), data = simdata),
                               specification = noblocksspec)

  expect_equal(noblocks_lmitt_fittedsbgrp@moderator, "force")
  expect_equal(blocked_lmitt_fittedsbgrp@moderator, "force")
  expect_equal(lmitt_fitted_nosbgrp@moderator, vector("character"))
  expect_equal(blocked_lmitt_fitted_nosbgrp@moderator, vector("character"))
  expect_equal(not_lmitt_fitted@moderator, vector("character"))
})

test_that("estfun.teeMod requires a certain model class", {
  data(simdata)

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification = spec)
  mod@.S3Class <- "not_a_mod"

  expect_error(estfun(mod), "must have been fitted using")
})

test_that(paste("estfun.teeMod returns original psi if no offset or no",
                "SandwichLayer offset"), {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  ca <- cov_adj(cmod, newdata = simdata, specification = spec)

  nolmittmod1 <- lm(y ~ z, simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec)
  nolmittmod2 <- lm(y ~ z, data = simdata, offset = ca)
  mod2 <- lmitt(y ~ 1, data = simdata, specification = spec, offset = stats::predict(cmod))

  ef_expected1 <- estfun(nolmittmod1)
  ef_expected2 <- estfun(nolmittmod2)

  expect_true(all.equal(estfun(mod1), ef_expected1, check.attributes = FALSE))
  expect_true(all.equal(estfun(mod2), ef_expected2, check.attributes = FALSE))
})

test_that("estfun.teeMod gets scaling constants right with no overlap", {
  set.seed(438)
  data(simdata)
  simdata_copy <- simdata
  nc <- 30
  nq <- nrow(simdata_copy)
  n <- nc + nq
  cmod_ssta <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_ssta)

  simdata_copy$id <- seq(nq)
  spec <- rct_spec(z ~ unitid(id), simdata_copy)
  dmod <- lmitt(y ~ 1, specification = spec, data = simdata_copy, offset = cov_adj(cmod))

  ef_pieces <- .align_and_extend_estfuns(dmod, by = "id")
  a11_inv <- .get_a11_inverse(dmod)
  a21 <- .get_a21(dmod)
  expect_equal(
    estfun(dmod),
    ef_pieces[["psi"]] - nq / nc * ef_pieces[["phi"]] %*% t(a11_inv) %*% t(a21)
  )
})

test_that("estfun.teeMod gets scaling constants right with partial overlap", {
  set.seed(438)
  data(simdata)
  simdata_copy <- simdata
  nq <- nrow(simdata_copy)
  simdata_copy$id <- seq(nq)

  aux_nc <- 30
  n <- nc <- nq + aux_nc
  aux_cmod_ssta <- data.frame(y = rnorm(aux_nc), x = rnorm(aux_nc), id = nq + seq(aux_nc))
  cmod <- lm(y ~ x, rbind(simdata_copy[, c("y", "x", "id")], aux_cmod_ssta))

  spec <- rct_spec(z ~ unitid(id), simdata_copy)
  dmod <- lmitt(y ~ 1, specification = spec, data = simdata_copy, offset = cov_adj(cmod))

  ef_pieces <- .align_and_extend_estfuns(dmod, by = "id")
  a11_inv <- .get_a11_inverse(dmod)
  a21 <- .get_a21(dmod)
  expect_equal(
    estfun(dmod),
    ef_pieces[["psi"]] - nq / nc * ef_pieces[["phi"]] %*% t(a11_inv) %*% t(a21)
  )
})

test_that("estfun.teeMod with missing values", {
  set.seed(438)
  data(simdata)
  simdata_copy <- simdata
  simdata_copy$z[1:3] <- NA_real_
  nc <- 30
  nq <- nrow(simdata_copy)
  n <- nc + nq
  cmod_ssta <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_ssta)
  
  simdata_copy$id <- seq(nq)
  spec <- rct_spec(z ~ unitid(id), simdata_copy)
  dmod <- lmitt(y ~ 1, specification = spec, data = simdata_copy, offset = cov_adj(cmod))
  
  ef <- estfun(dmod)
  expect_equal(nrow(ef), nc + nq)
  expect_true(all.equal(ef[1:3,], matrix(0, nrow = 3, ncol = 2), check.attributes = FALSE))
  class(dmod$na.action) <- "exclude"
  ef_pieces <- .align_and_extend_estfuns(dmod, by = "id")
  a11_inv <- .get_a11_inverse(dmod)
  a21 <- .get_a21(dmod)
  expect_equal(
    ef,
    ef_pieces[["psi"]] - nq / nc * ef_pieces[["phi"]] %*% t(a11_inv) %*% t(a21)
  )
})

test_that(paste("estfun.teeMod returns correct dimensions and alignment",
                "when exact alignment between C and Q is possible"), {
  set.seed(438)
  data(simdata)
  simdata_copy <- simdata
  simdata_copy$uid <- rownames(simdata_copy)

  shuffled_simdata <- simdata_copy[sample(rownames(simdata_copy)),]
  cmod <- lm(y ~ x, simdata_copy)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata_copy)
  mod1 <- lmitt(y ~ 1, data = simdata_copy, specification = spec, offset = cov_adj(cmod, by = "uid"))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, specification = spec, offset = cov_adj(cmod, by = "uid"))

  expect_equal(dim(estfun(mod1)), c(nrow(simdata_copy), 2))
  expect_equal(estfun(mod1), estfun(mod2))
})

test_that(paste("estfun.teeMod returns correct dimensions when only",
                "inexact alignment between C and Q is possible"), {
  set.seed(438)
  data(simdata)

  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, specification = spec, offset = cov_adj(cmod))

  expect_equal(dim(estfun(mod1)), c(nrow(simdata), 2))
  expect_equal(dim(estfun(mod2)), c(nrow(simdata), 2))
})

test_that(paste("estfun.teeMod returns correct dimensions and alignment",
                "when exact alignment between C and Q is possible but they only",
                "partially overlap"), {
  set.seed(438)
  data(simdata)
  simdata$uid <- rownames(simdata)

  Q_data <- simdata[simdata$uid %in% seq_len(20),]
  shuffled_Q_data <- simdata[sample(rownames(Q_data)),]
  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, specification = spec, offset = cov_adj(cmod, by = "uid"))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, specification = spec, offset = cov_adj(cmod, by = "uid"))

  expect_equal(dim(estfun(mod1)), c(nrow(simdata), 2))
  expect_equal(estfun(mod1), estfun(mod2))
})

test_that(paste("estfun.teeMod returns correct dimensions for partial",
                "overlap of C and Q when only inexact alignment is possible"), {
  data(simdata)
  set.seed(401)
  cmod_ssta <- rbind(
    simdata[simdata$uoa1 == 1, c("x", "y", "uoa1", "uoa2")],
    data.frame("x" = rnorm(100), "y" = rnorm(100), "uoa1" = NA, "uoa2" = NA)
  )

  cmod <- lm(y ~ x, cmod_ssta)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))

  expect_equal(dim(estfun(mod)),
               c(nrow(simdata) + nrow(cmod_ssta[is.na(cmod_ssta$uoa1),]), 2))
})

test_that(paste("estfun.teeMod returns correct dimensions for no",
                "overlap of C and Q when only inexact alignment is possible"), {
  data(simdata)
  set.seed(401)
  cmod_ssta <- data.frame(
    "x" = rnorm(100), "y" = rnorm(100), "uoa1" = NA, "uoa2" = NA
  )

  cmod <- lm(y ~ x, cmod_ssta)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))

  expect_equal(dim(estfun(mod)), c(nrow(simdata) + nrow(cmod_ssta), 2))
})

if (requireNamespace("robustbase", quietly = TRUE)) {
  test_that("estfun.teeMod returns correct dimensions for rectangular A11_inv", {
    data(simdata)
    cmod <- robustbase::lmrob(y ~ x, simdata)
    spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
    dmod <- lmitt(lm(y ~ z.(spec), simdata, offset = cov_adj(cmod)), specification = spec)

    out <- estfun(dmod)
    expect_equal(dim(out), c(50, 2))
  })
}

test_that(paste("bread.teeMod returns bread.lm for teeMod objects",
                "with non-SandwichLayer or NULL offsets"), {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  ca <- predict(cmod)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  m1 <- lmitt(y ~ 1, data = simdata, specification = spec)
  m2 <- lmitt(y ~ 1, data = simdata, specification = spec, offset = ca)

  expected_out1 <- nrow(simdata) * chol2inv(m1$qr$qr)
  coef_names <- names(m1$coefficients)
  dimnames(expected_out1) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(m1), expected_out1)

  expected_out2 <- nrow(simdata) * chol2inv(m2$qr$qr)
  dimnames(expected_out2) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(m2), expected_out2)
})

test_that(paste("bread.teeMod returns expected output for full overlap",
                "of C and Q"), {
  data(simdata)
  simdata$uid <- seq_len(nrow(simdata))

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ uoa(uid), simdata)
  m <- lmitt(y ~ 1, data = simdata, specification = spec, weights = ate(spec),
             offset = cov_adj(cmod))

  expected_out <- nrow(simdata) * chol2inv(m$qr$qr)
  coef_names <- names(m$coefficients)
  dimnames(expected_out) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(m), expected_out)
})

test_that(paste("bread.teeMod returns expected output for partial overlap",
                "of C and Q"), {
  set.seed(438)
  data(simdata)
  nq <- nrow(simdata)
  simdata$id <- seq(nq)

  aux_nc <- 30
  n <- nc <- nq + aux_nc
  aux_cmod_ssta <- data.frame(y = rnorm(aux_nc), x = rnorm(aux_nc), id = nq + seq(aux_nc))
  cmod <- lm(y ~ x, rbind(simdata[, c("y", "x", "id")], aux_cmod_ssta))

  spec <- rct_spec(z ~ unitid(id), simdata)
  dmod <- lmitt(y ~ 1, specification = spec, data = simdata, offset = cov_adj(cmod))

  expected_out <- n * chol2inv(dmod$qr$qr)
  coef_names <- names(dmod$coefficients)
  dimnames(expected_out) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(dmod), expected_out)
})

test_that(paste("bread.teeMod returns expected output for no overlap",
                "of C and Q"), {
  set.seed(438)
  data(simdata)
  nc <- 30
  nq <- nrow(simdata)
  n <- nc + nq
  cmod_ssta <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_ssta)

  simdata$id <- seq(nq)
  spec <- rct_spec(z ~ unitid(id), simdata)
  dmod <- lmitt(y ~ 1, specification = spec, data = simdata, offset = cov_adj(cmod))

  expected_out <- n * chol2inv(dmod$qr$qr)
  coef_names <- names(dmod$coefficients)
  dimnames(expected_out) <- list(coef_names, coef_names)
  expect_equal(sandwich::bread(dmod), expected_out)
})

test_that("bread.teeMod handles model with less than full rank", {
  data(simdata)
  copy_simdata <- simdata
  copy_simdata$o_fac <- as.factor(copy_simdata$o)
  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), copy_simdata)

  ### lmitt.formula
  ssmod <- lmitt(y ~ o_fac, data = copy_simdata, specification = spec, offset = cov_adj(cmod))
  keep_ix <- ssmod$qr$pivot[1L:ssmod$rank]
  expect_true(
    all.equal(
      bread(ssmod),
      nrow(simdata) * solve(crossprod(model.matrix(ssmod))[keep_ix,keep_ix])
    )
  )
})

test_that(paste(".align_and_extend_estfuns fails if not a teeMod object",
                "with a SandwichLayer offset"), {
  data(simdata)

  mod1 <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  mod2 <- lmitt(y ~ 1, specification = spec, data = simdata)
  expect_error(.align_and_extend_estfuns(mod1), "must be a fitted")
  expect_error(.align_and_extend_estfuns(mod2), "must be a fitted")
})

test_that(paste(".align_and_extend_estfuns with `by` and the samples fully overlap"), {
  set.seed(438)
  data(simdata)
  simdata_copy <- simdata

  simdata_copy$obs_id <- seq_len(nrow(simdata_copy))
  shuffle_ix <- sample(rownames(simdata_copy))
  shuffled_simdata <- simdata_copy[shuffle_ix,]
  cmod1 <- lm(y ~ x, simdata_copy)
  cmod2 <- lm(y ~ x, shuffled_simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata_copy)
  mod1 <- lmitt(y ~ 1, data = simdata_copy, specification = spec, offset = cov_adj(cmod1, by = "obs_id"))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, specification = spec, offset = cov_adj(cmod2, by = "obs_id"))

  cmod_ef <- estfun(cmod1)
  mod_lm_ef <- estfun(as(mod1, "lm"))
  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)

  # tests to run (for each .align_and_extend_estfuns() test):
  # 1) do we get a matrix of estimating equations of dimension n?
  # 2) are rows' contributions to the matrix of estimating equations aligned in
  #     phi and psi?
  # 3) are variance estimates invariant to the original order of the data,
  #     even when we change the clustering level?
  expect_equal(dim(ef1$phi), c(nrow(simdata_copy), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata_copy), 2))
  expect_true(all.equal(
    ef1$phi,
    cmod_ef[sort(as.character(seq(nrow(simdata_copy)))),],
    check.attributes = FALSE))
  expect_true(all.equal(
    ef1$psi,
    mod_lm_ef[sort(as.character(seq(nrow(simdata_copy)))),],
    check.attributes = FALSE))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1, cluster = "bid"), vcov_tee(mod2, cluster = "bid"))
})

test_that(paste(".align_and_extend_estfuns with `by` and Q is a subset of C"), {
  set.seed(438)
  data(simdata)

  simdata$obs_id <- seq_len(nrow(simdata))
  simdata_shuffle_ix <- sample(rownames(simdata))
  shuffled_simdata <- simdata[simdata_shuffle_ix,]
  Q_data <- simdata[1L:20L,]
  Q_shuffle_ix <- sample(rownames(Q_data))
  shuffled_Q_data <- Q_data[Q_shuffle_ix,]
  cmod1 <- lm(y ~ x, simdata)
  cmod2 <- lm(y ~ x, shuffled_simdata)
  spec1 <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_data)
  spec2 <- rct_spec(z ~ cluster(uoa1, uoa2), data = shuffled_Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, specification = spec1, offset = cov_adj(cmod1, by = "obs_id"))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, specification = spec2, offset = cov_adj(cmod2, by = "obs_id"))

  cmod_ef <- estfun(cmod1)
  mod_lm_ef <- estfun(as(mod1, "lm"))
  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  nonzero_ix <- as.numeric(sort(as.character(seq(20L))))
  zero_ix <- seq(21L, 50L)
  expect_true(all.equal(ef1$phi, cmod_ef[c(nonzero_ix, zero_ix),], check.attributes = FALSE))
  expect_true(all.equal(ef1$psi[1L:20L,], mod_lm_ef[nonzero_ix,], check.attributes = FALSE))
  expect_true(all(ef1$psi[zero_ix] == 0))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1, cluster = "bid"), vcov_tee(mod2, cluster = "bid"))
})

test_that(paste(".align_and_extend_estfuns with `by` and C is a subset of Q"), {
  set.seed(438)
  data(simdata)
  simdata_copy <- simdata

  simdata_copy$obs_id <- seq_len(nrow(simdata_copy))
  simdata_shuffle_ix <- sample(rownames(simdata_copy))
  shuffled_simdata <- simdata_copy[simdata_shuffle_ix,]
  C_data <- simdata_copy[1:20,]
  C_shuffle_ix <- sample(rownames(C_data))
  shuffled_C_data <- C_data[C_shuffle_ix,]
  cmod1 <- lm(y ~ x, C_data)
  cmod2 <- lm(y ~ x, shuffled_C_data)
  spec1 <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata_copy)
  spec2 <- rct_spec(z ~ cluster(uoa1, uoa2), data = shuffled_simdata)
  mod1 <- lmitt(y ~ 1, data = simdata_copy, specification = spec1, offset = cov_adj(cmod1, by = "obs_id"))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, specification = spec2, offset = cov_adj(cmod2, by = "obs_id"))

  cmod_ef <- estfun(cmod1)
  mod_lm_ef <- estfun(as(mod1, "lm"))
  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)

  expect_equal(dim(ef1$phi), c(nrow(simdata_copy), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata_copy), 2))
  nonzero_ix <- 31L:50L
  zero_ix <- 1L:30L
  expect_true(all(ef1$phi[1L:30L] == 0))
  expect_true(all.equal(
    ef1$phi[31L:50L,],
    cmod_ef[sort(as.character(seq(20L))),],
    check.attributes = FALSE))
  expect_true(all.equal(
    ef1$psi,
    mod_lm_ef[c(21L:50L, sort(as.character(seq(20L)))),], check.attributes = FALSE))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1, cluster = "bid"), vcov_tee(mod2, cluster = "bid"))
})

test_that(paste(".align_and_extend_estfuns with `by` and C and Q have no overlap"), {
  set.seed(438)
  data(simdata)

  simdata$obs_id <- seq_len(nrow(simdata))
  Q_data <- simdata[21:50,]
  C_data <- simdata[1:20,]
  Q_shuffle_ix <- sample(rownames(Q_data))
  shuffled_Q_data <- Q_data[Q_shuffle_ix,]
  C_shuffle_ix <- sample(rownames(C_data))
  shuffled_C_data <- C_data[C_shuffle_ix,]
  cmod1 <- lm(y ~ x, C_data)
  cmod2 <- lm(y ~ x, shuffled_C_data)
  spec1 <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_data)
  spec2 <- rct_spec(z ~ cluster(uoa1, uoa2), data = shuffled_Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, specification = spec1, offset = cov_adj(cmod1, by = "obs_id"))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, specification = spec2, offset = cov_adj(cmod2, by = "obs_id"))

  cmod_ef <- estfun(cmod1)
  mod_lm_ef <- estfun(as(mod1, "lm"))
  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_true(all(ef1$phi[1L:30L,] == 0))
  expect_true(all.equal(ef1$phi[31L:50L,], cmod_ef[1L:20L,], check.attributes = FALSE))
  expect_true(all(ef1$psi[31L:50L,] == 0))
  expect_true(all.equal(ef1$psi[1L:30L,], mod_lm_ef[1L:30L,], check.attributes = FALSE))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1, cluster = "bid"), vcov_tee(mod2, cluster = "bid"))
})

test_that(paste(".align_and_extend_estfuns when the samples fully overlap (no `by`)"), {
  set.seed(438)
  data(simdata)

  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  cmod1 <- lm(y ~ x, simdata)
  cmod2 <- lm(y ~ x, shuffled_simdata)
  spec1 <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  spec2 <- rct_spec(z ~ cluster(uoa1, uoa2), data = shuffled_simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec1, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, specification = spec2, offset = cov_adj(cmod2))

  cmod_ef <- estfun(cmod1)
  mod_lm_ef <- estfun(as(mod1, "lm"))
  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)
  by_ix <- sort(apply(simdata[, c("uoa1", "uoa2")], 1,
                      function(...) paste(..., collapse = "_")))

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_true(all.equal(ef1$phi, cmod_ef[sort(seq(50L)),], check.attributes = FALSE))
  expect_true(all.equal(ef1$psi, mod_lm_ef[sort(seq(50L)),], check.attributes = FALSE))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1, cluster = "bid"), vcov_tee(mod2, cluster = "bid"))
})

test_that(paste(".align_and_extend_estfuns when Q is a subset of C (no `by`)"), {
  set.seed(438)
  data(simdata)

  shuffled_simdata <- simdata[sample(rownames(simdata)),]
  Q_data <- simdata[1:20,]
  shuffled_Q_data <- Q_data[sample(rownames(Q_data)),]
  cmod1 <- lm(y ~ x, simdata)
  cmod2 <- lm(y ~ x, shuffled_simdata)
  spec1 <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_data)
  spec2 <- rct_spec(z ~ cluster(uoa1, uoa2), data = shuffled_Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, specification = spec1, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, specification = spec2, offset = cov_adj(cmod2))

  cmod_ef <- estfun(cmod1)
  mod_lm_ef <- estfun(as(mod1, "lm"))
  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)
  by_ix <- sort(apply(simdata[, c("uoa1", "uoa2")], 1,
                      function(...) paste(..., collapse = "_")))

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_true(all.equal(ef1$phi, cmod_ef, check.attributes = FALSE))
  expect_true(all.equal(ef1$psi[1L:20L,], mod_lm_ef[1L:20L,], check.attributes = FALSE))
  expect_true(all(ef1$psi[21L:50L] == 0))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1, cluster = "bid"), vcov_tee(mod2, cluster = "bid"))
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
  spec1 <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  spec2 <- rct_spec(z ~ cluster(uoa1, uoa2), data = shuffled_simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec1, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = shuffled_simdata, specification = spec2, offset = cov_adj(cmod2))

  cmod_ef <- estfun(cmod1)
  mod_lm_ef <- estfun(as(mod1, "lm"))
  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)
  by_ix <- sort(apply(simdata[, c("uoa1", "uoa2")], 1,
                      function(...) paste(..., collapse = "_")))

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_true(all.equal(ef1$phi[31L:50L,], cmod_ef[1L:20L,], check.attributes = FALSE))
  expect_true(all(ef1$phi[1L:30L] == 0))
  expect_true(all.equal(ef1$psi, mod_lm_ef[c(21L:50L, 1L:20L),], check.attributes = FALSE))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1, cluster = "bid"), vcov_tee(mod2, cluster = "bid"))
})

test_that(paste(".align_and_extend_estfuns when the samples have no overlap (no `by`)"), {
  set.seed(438)
  data(simdata)

  Q_data <- simdata[21:50,]
  C_data <- simdata[1:20,]
  shuffled_Q_data <- Q_data[sample(rownames(Q_data)),]
  shuffled_C_data <- C_data[sample(rownames(C_data)),]
  cmod1 <- lm(y ~ x, C_data)
  cmod2 <- lm(y ~ x, shuffled_C_data)
  spec1 <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_data)
  spec2 <- rct_spec(z ~ cluster(uoa1, uoa2), data = shuffled_Q_data)
  mod1 <- lmitt(y ~ 1, data = Q_data, specification = spec1, offset = cov_adj(cmod1))
  mod2 <- lmitt(y ~ 1, data = shuffled_Q_data, specification = spec2, offset = cov_adj(cmod2))

  cmod_ef <- estfun(cmod1)
  mod_lm_ef <- estfun(as(mod1, "lm"))
  ef1 <- .align_and_extend_estfuns(mod1)
  ef2 <- .align_and_extend_estfuns(mod2)
  by_ix <- sort(apply(simdata[, c("uoa1", "uoa2")], 1,
                      function(...) paste(..., collapse = "_")))

  expect_equal(dim(ef1$phi), c(nrow(simdata), 2))
  expect_equal(dim(ef1$psi), c(nrow(simdata), 2))
  expect_true(all.equal(ef1$phi[31L:50L,], cmod_ef[1L:20L,], check.attributes = FALSE))
  expect_true(all(ef1$phi[1L:30L,] == 0))
  expect_true(all.equal(ef1$psi[1L:30L,], mod_lm_ef[1L:30L,], check.attributes = FALSE))
  expect_true(all(ef1$psi[31L:50L,] == 0))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1, cluster = "bid"), vcov_tee(mod2, cluster = "bid"))
})

test_that(".make_uoa_ids fails without cluster argument or teeMod model", {
  data(simdata)
  mod <- lm(y ~ z, data = simdata)
  expect_error(.make_uoa_ids(mod), "Must be provided")
})

test_that(".make_uoa_ids returns correct ID's for non-SandwichLayer offset", {
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification = spec)

  expected_out <- factor(
    apply(simdata[, c("uoa1", "uoa2"), drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  )
  expect_equal(.make_uoa_ids(mod, vcov_type = "CR"), expected_out)
})

test_that(".make_uoa_ids returns correct ID's for full overlap of C and Q", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  dmod <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))

  expected_out <- factor(
    apply(simdata[, c("uoa1", "uoa2"), drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  )
  out <- .make_uoa_ids(dmod, vcov_type = "CR")

  expect_true(all.equal(out, expected_out, check.attributes = FALSE))
})

test_that(".make_uoa_ids for no overlap of C and Q", {
  data(simdata)
  set.seed(300)
  cmod_ssta1 <- data.frame("y" = rnorm(10), "x" = rnorm(10), "uoa1" = NA, "uoa2" = NA)
  cmod_ssta2 <- data.frame("y" = rnorm(10), "x" = rnorm(10),
                           "uoa1" = rep(c(1, 2), each = 5), "uoa2" = NA)

  cmod1 <- lm(y ~ x, cmod_ssta1)
  cmod2 <- lm(y ~ x, cmod_ssta2)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  dmod1 <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod1))
  dmod2 <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod2))

  Q_uoas <- apply(simdata[, c("uoa1", "uoa2"), drop = FALSE], 1,
                  function(...) paste(..., collapse = "_"))

  ids1 <- .make_uoa_ids(dmod1, vcov_type = "CR")
  ids2 <- .make_uoa_ids(dmod2, vcov_type = "CR")

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
  C_not_Q <- data.frame("y" = rnorm(20), "x" = rnorm(20), "uoa1" = NA, "uoa2" = NA)
  cmod_ssta <- rbind(simdata[, colnames(C_not_Q)], C_not_Q)

  cmod <- lm(y ~ x, cmod_ssta)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  dmod <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))

  Q_uoas <- apply(simdata[, c("uoa1", "uoa2"), drop = FALSE], 1,
                  function(...) paste(..., collapse = "_"))

  ids <- .make_uoa_ids(dmod, vcov_type = "CR")

  expect_true(is.factor(ids))

  expect_equal(length(ids), nrow(simdata) + 20)

  expect_true(all.equal(ids[1:nrow(simdata)], factor(Q_uoas), check.attributes = FALSE))

  expect_equal(length(unique(ids)), length(unique(Q_uoas)) + 20)
})

test_that(paste(".make_uoa_ids returns correct ID's when cov_adj's 'by' argument",
                "provided a different ordering"), {
  set.seed(300)
  data(simdata)
  simdata_copy <- simdata
  simdata_copy$id <- sample(seq_len(nrow(simdata_copy)))

  C_not_Q <- data.frame("y" = rnorm(20), "x" = rnorm(20), "uoa1" = NA, "uoa2" = NA,
                        "id" = seq(51, 70))
  cmod_ssta <- rbind(simdata_copy[, colnames(C_not_Q)], C_not_Q)

  cmod <- lm(y ~ x, cmod_ssta)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata_copy)
  dmod <- lmitt(y ~ 1, data = simdata_copy, specification = spec, offset = cov_adj(cmod, by = "id"))

  Q_uoas <- apply(simdata_copy[, c("uoa1", "uoa2"), drop = FALSE], 1,
                  function(...) paste(..., collapse = "_"))

  ids <- .make_uoa_ids(dmod, vcov_type = "CR")

  expect_true(is.factor(ids))
  expect_equal(length(ids), nrow(simdata_copy) + 20)
  expect_true(all.equal(
    ids[1L:nrow(simdata_copy)],
    factor(Q_uoas)[match(sort(as.character(simdata_copy$id)), simdata_copy$id)],
    check.attributes = FALSE))
  expect_equal(length(unique(ids)), length(unique(Q_uoas)) + 20)
})

test_that(".make_uoa_ids with `by` when Q ID's aren't unique", {
  set.seed(23)
  cmod_ssta <- data.frame(id = letters[1:15],
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_ssta)

  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)

  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
  mod_w_by <- lmitt(y ~ yr, specification = newspec, data = analysis_dat,
                    offset = cov_adj(cmod, by = "by_col"))

  by_overlap <- sort(as.character(setdiff(seq_len(15), seq(1, 15, 3))))
  nonoverlap <- as.character(seq(1, 15, 3))
  overlap_Q_ix <- match(by_overlap, analysis_dat$by_col)
  overlap_C_ix <- match(by_overlap, cmod_ssta$by_col)

  out <- c(analysis_dat$id[overlap_Q_ix],
           cmod_ssta$id[match(nonoverlap, cmod_ssta$by_col)])
  expect_equal(length(ids <- .make_uoa_ids(mod, vcov_type = "HC0", cluster = "id",
                                           by = "by_col")),
               15)
  expect_equal(ids, factor(out, levels = unique(out)))
  expect_equal(ids, .make_uoa_ids(mod_w_by, vcov_type = "HC0", cluster = "id",
                                  by = "by_col"))
})

test_that(".make_uoa_ids with `by` when C ID's aren't unique", {
  set.seed(23)
  cmod_ssta <- data.frame(yr = rep(c("00", "01", "02"), 5),
                          id = rep(letters[1:5], each = 3),
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_ssta)

  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)

  analysis_dat <- data.frame(id = letters[1:5],
                             x = rnorm(5),
                             a = c(rep(1, 3), rep(0, 2)),
                             y = rnorm(5),
                             by_col = seq(3, 15, 3))
  mod <- lmitt(y ~ 1, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
  mod_w_by <- lmitt(y ~ 1, specification = newspec, data = analysis_dat,
                    offset = cov_adj(cmod, by = "by_col"))

  by_overlap <- sort(as.character(seq(3, 15, 3)))
  nonoverlap <- as.character(setdiff(seq_len(15), seq(3, 15, 3)))
  overlap_Q_ix <- match(by_overlap, analysis_dat$by_col)
  overlap_C_ix <- match(by_overlap, cmod_ssta$by_col)

  out <- c(analysis_dat$id[overlap_Q_ix],
           cmod_ssta$id[match(nonoverlap, cmod_ssta$by_col)])
  expect_equal(length(ids <- .make_uoa_ids(mod, vcov_type = "HC0", cluster = "id",
                                           by = "by_col")),
               15)
  expect_equal(ids, factor(out, levels = unique(out)))
  expect_equal(ids, .make_uoa_ids(mod_w_by, vcov_type = "HC0", cluster = "id",
                                  by = "by_col"))
})

test_that(".make_uoa_ids with `by` when Q and C ID's aren't unique", {
  set.seed(23)
  cmod_ssta <- data.frame(yr = rep(c("00", "01", "02"), 5),
                          id = rep(letters[1:5], each = 3),
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_ssta)

  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)

  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
  mod_w_by <- lmitt(y ~ yr, specification = newspec, data = analysis_dat,
                    offset = cov_adj(cmod, by = "by_col"))

  by_overlap <- sort(as.character(setdiff(seq_len(15), seq(1, 15, 3))))
  nonoverlap <- as.character(seq(1, 15, 3))
  overlap_Q_ix <- match(by_overlap, analysis_dat$by_col)
  overlap_C_ix <- match(by_overlap, cmod_ssta$by_col)

  out <- c(analysis_dat$id[overlap_Q_ix],
           cmod_ssta$id[match(nonoverlap, cmod_ssta$by_col)])
  expect_equal(length(ids <- .make_uoa_ids(mod, vcov_type = "HC0", cluster = "id",
                                           by = "by_col")),
               15)
  expect_equal(ids, factor(out, levels = unique(out)))
  expect_equal(ids, .make_uoa_ids(mod_w_by, vcov_type = "HC0", cluster = "id",
                                  by = "by_col"))
})

test_that(".make_uoa_ids without `by` when Q and C ID's aren't unique", {
  set.seed(23)
  cmod_ssta <- data.frame(yr = rep(c("00", "01", "02"), 5),
                          id = rep(letters[1:5], each = 3),
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_ssta)

  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)

  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
  expect_error(.make_uoa_ids(mod, vcov_type = "HC0", cluster = "id"),
               "not uniquely specified. Provide a `by` argument")
})

test_that("model-based .make_uoa_ids replaces uoa ID's in small blocks with block ID's", {
  data(simdata)

  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  suppressMessages(dmod <- lmitt(y ~ 1, specification = spec, data = simdata))

  expect_equal(
    out <- .make_uoa_ids(dmod, "CR"),
    factor(paste0("bid",
                  c(rep("1", sum(simdata$bid == 1)), rep("2", sum(simdata$bid == 2)),
                    rep("3", sum(simdata$bid == 3)))),
                  levels = paste0("bid", c("1", "2", "3")))
  )
  expect_equal(out, .make_uoa_ids(dmod, "MB"))
  expect_equal(out, .make_uoa_ids(dmod, "HC"))
  expect_equal(
    .make_uoa_ids(dmod, "DB"),
    factor(apply(simdata[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_")))
  )
})

test_that(paste(".order_samples fails without a teeMod object or",
                "SandwichLayer offset"), {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, simdata, offset = predict(cmod))
  mod2 <- lmitt(y ~ 1, data = simdata, specification = spec)

  expect_error(.order_samples(mod1), "must be a teeMod object")
  expect_error(.order_samples(mod2), "must be a teeMod object")
})

test_that(".order_samples when the samples fully overlap", {
  set.seed(300)
  data(simdata)
  simdata_copy <- simdata

  simdata_copy$uid <- seq_len(nrow(simdata_copy))
  cmod <- lm(y ~ x, simdata_copy)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata_copy)
  mod <- lmitt(y ~ 1, data = simdata_copy, specification = spec, offset = cov_adj(cmod, by = "uid"))

  out <- .order_samples(mod)

  expect_equal(length(out$Q_not_C), 0)
  expect_equal(length(out$Q_in_C), nrow(simdata_copy))
  expect_equal(length(out$C_in_Q), nrow(simdata_copy))
  expect_equal(length(out$C_not_Q), 0)

  expect_equal(names(out$Q_in_C), ord <- sort(as.character(seq(nrow(simdata_copy)))))
  expect_equal(names(out$C_in_Q), ord)

  expect_true(all.equal(out$Q_in_C, ord, check.attributes = FALSE))
  expect_true(all.equal(out$C_in_Q, ord, check.attributes = FALSE))
})

test_that(".order_samples when Q is a subset of C", {
  set.seed(300)
  data(simdata)
  simdata_copy <- simdata

  simdata_copy$uid <- seq_len(nrow(simdata_copy))
  cmod_ssta <- rbind(
    simdata_copy[, c("x", "y", "uoa1", "uoa2", "uid")],
    data.frame(x = rnorm(30), y = rnorm(30), uoa1 = NA, uoa2 = NA,
               uid = seq_len(30) + nrow(simdata_copy))
  )
  cmod <- lm(y ~ x, cmod_ssta)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata_copy)
  mod <- lmitt(y ~ 1, data = simdata_copy, specification = spec, offset = cov_adj(cmod, by = "uid"))

  out <- .order_samples(mod)

  expect_equal(length(out$Q_not_C), 0)
  expect_equal(length(out$Q_in_C), nrow(simdata_copy))
  expect_equal(length(out$C_in_Q), nrow(simdata_copy))
  expect_equal(length(out$C_not_Q), 30)

  expect_equal(names(out$Q_in_C), overlap_ord <- sort(as.character(seq(nrow(simdata_copy)))))
  expect_equal(names(out$C_in_Q), overlap_ord)
  expect_equal(names(out$C_not_Q), as.character(seq(51, 80)))

  expect_true(all.equal(out$Q_in_C, overlap_ord, check.attributes = FALSE))
  expect_true(all.equal(out$C_in_Q, overlap_ord, check.attributes = FALSE))
  expect_true(all.equal(out$C_not_Q, as.character(seq(51, 80)), check.attributes = FALSE))
})


test_that(".order_samples when C is a subset of Q", {
  set.seed(300)
  data(simdata)
  simdata_copy <- simdata

  simdata_copy$uid <- seq_len(nrow(simdata_copy))
  cmod_ssta <- simdata_copy[1:20,]
  cmod <- lm(y ~ x, cmod_ssta)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata_copy)
  mod <- lmitt(y ~ 1, data = simdata_copy, specification = spec, offset = cov_adj(cmod, by = "uid"))

  out <- .order_samples(mod)

  expect_equal(length(out$Q_not_C), 30)
  expect_equal(length(out$Q_in_C), 20)
  expect_equal(length(out$C_in_Q), 20)
  expect_equal(length(out$C_not_Q), 0)

  expect_equal(names(out$Q_not_C), as.character(seq(21, 50)))
  expect_equal(names(out$Q_in_C), overlap_ord <- sort(as.character(seq(20))))
  expect_equal(names(out$C_in_Q), overlap_ord)

  expect_true(all.equal(out$Q_not_C, as.character(seq(21, nrow(simdata_copy))), check.attributes = FALSE))
  expect_true(all.equal(out$Q_in_C, overlap_ord, check.attributes = FALSE))
  expect_true(all.equal(out$C_in_Q, overlap_ord, check.attributes = FALSE))
})

test_that(".order_samples when the samples do not overlap", {
  set.seed(300)
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  cmod_ssta <- simdata[1:20,]
  Q_data <- simdata[21:50,]
  cmod <- lm(y ~ x, cmod_ssta)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_data)
  mod <- lmitt(y ~ 1, data = Q_data, specification = spec, offset = cov_adj(cmod, by = "uid"))

  out <- .order_samples(mod)

  expect_equal(length(out$Q_not_C), 30)
  expect_equal(length(out$Q_in_C), 0)
  expect_equal(length(out$C_in_Q), 0)
  expect_equal(length(out$C_not_Q), 20)

  expect_equal(names(out$Q_not_C), as.character(seq(30)))
  expect_equal(names(out$C_not_Q), as.character(seq(20)))

  expect_true(all.equal(out$Q_not_C, as.character(seq(21, 50)), check.attributes = FALSE))
  expect_true(all.equal(out$C_not_Q, as.character(seq(20)), check.attributes = FALSE))
})

test_that(".order_samples when no `by` argument provided", {
  set.seed(300)
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  C_data <- simdata[1:20,]
  Q_data <- simdata[21:50,]
  cmod <- lm(y ~ x, C_data)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_data)
  mod <- lmitt(y ~ 1, data = Q_data, specification = spec, offset = cov_adj(cmod))

  out <- .order_samples(mod)

  Q_uoas <- apply(Q_data[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))
  C_uoas <- apply(C_data[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))

  expect_equal(length(out$Q_not_C), 30)
  expect_equal(length(out$Q_in_C), 0)
  expect_equal(length(out$C_in_Q), 0)
  expect_equal(length(out$C_not_Q), 20)

  expect_true(all.equal(names(out$Q_not_C), as.character(seq(30)), check.attributes = FALSE))
  expect_true(all.equal(names(out$C_not_Q), as.character(seq(20)), check.attributes = FALSE))

  expect_true(all.equal(out$Q_not_C, Q_uoas, check.attributes = FALSE))
  expect_true(all.equal(out$C_not_Q, C_uoas, check.attributes = FALSE))
})

test_that(".order_samples with `by` when Q ID's aren't unique", {
  set.seed(23)
  cmod_ssta <- data.frame(id = letters[1:15],
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_ssta)

  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)

  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
  mod_w_by <- lmitt(y ~ yr, specification = newspec, data = analysis_dat,
                    offset = cov_adj(cmod, by = "by_col"))

  overlap <- sort(as.character(setdiff(seq_len(15), seq(1, 15, 3))))
  nonoverlap <- as.character(seq(1, 15, 3))
  out <- list("Q_not_C" = setNames(character(0L), character(0L)),
              "Q_in_C" = setNames(overlap, match(overlap, analysis_dat$by_col)),
              "C_in_Q" = setNames(overlap, match(overlap, cmod_ssta$by_col)),
              "C_not_Q" = setNames(nonoverlap, which(cmod_ssta$by_col %in% nonoverlap)))
  expect_equal(ord <- .order_samples(mod, by = "by_col"), out)
  expect_equal(ord, .order_samples(mod_w_by))
})

test_that(".order_samples with `by` when C ID's aren't unique", {
  set.seed(23)
  cmod_ssta <- data.frame(yr = rep(c("00", "01", "02"), 5),
                          id = rep(letters[1:5], each = 3),
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_ssta)

  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)

  analysis_dat <- data.frame(id = letters[1:5],
                             x = rnorm(5),
                             a = c(rep(1, 3), rep(0, 2)),
                             y = rnorm(5),
                             by_col = seq(3, 15, 3))
  mod <- lmitt(y ~ 1, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
  mod_w_by <- lmitt(y ~ 1, specification = newspec, data = analysis_dat,
                    offset = cov_adj(cmod, by = "by_col"))

  overlap <- sort(as.character(seq(3, 15, 3)))
  nonoverlap <- as.character(setdiff(seq_len(15), seq(3, 15, 3)))
  out <- list("Q_not_C" = setNames(character(0L), character(0L)),
              "Q_in_C" = setNames(overlap, match(overlap, analysis_dat$by_col)),
              "C_in_Q" = setNames(overlap, match(overlap, cmod_ssta$by_col)),
              "C_not_Q" = setNames(nonoverlap, which(cmod_ssta$by_col %in% nonoverlap)))
  expect_equal(ord <- .order_samples(mod, by = "by_col"), out)
  expect_equal(.order_samples(mod_w_by), ord)
})

test_that(".order_samples with `by` when both Q and C ID's aren't unique", {
  set.seed(23)
  cmod_ssta <- data.frame(yr = rep(c("00", "01", "02"), 5),
                          id = rep(letters[1:5], each = 3),
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_ssta)

  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)

  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
  mod_w_by <- lmitt(y ~ yr, specification = newspec, data = analysis_dat,
                    offset = cov_adj(cmod, by = "by_col"))

  overlap <- sort(as.character(setdiff(seq_len(15), seq(1, 15, 3))))
  nonoverlap <- as.character(seq(1, 15, 3))
  out <- list("Q_not_C" = setNames(character(0L), character(0L)),
              "Q_in_C" = setNames(overlap, match(overlap, analysis_dat$by_col)),
              "C_in_Q" = setNames(overlap, match(overlap, cmod_ssta$by_col)),
              "C_not_Q" = setNames(nonoverlap, which(cmod_ssta$by_col %in% nonoverlap)))
  expect_equal(ord <- .order_samples(mod, by = "by_col"), out)
  expect_equal(.order_samples(mod_w_by), ord)
})

test_that(".order_samples without `by` when Q or C ID's aren't unique", {
  set.seed(23)
  cmod_ssta <- data.frame(yr = rep(c("00", "01", "02"), 5),
                          id = rep(letters[1:5], each = 3),
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_ssta)

  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)

  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
  expect_error(.order_samples(mod),
               "not uniquely specified. Provide a `by` argument")
})

test_that("sanitize_Q_ids succeeds with valid `by` argument", {
  data(simdata)

  simdata_copy <- simdata
  simdata_copy$uid <- seq_len(nrow(simdata_copy))
  cmod <- lm(y ~ z, data = simdata_copy)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2, uid), simdata_copy)
  dmod <- lmitt(y ~ 1, specification = spec, data = simdata_copy, offset = cov_adj(cmod))
  out <- .sanitize_Q_ids(dmod, c("uoa1", "uoa2"))$cluster

  expected_out <- apply(simdata_copy[, c("uoa1", "uoa2"), drop = FALSE], 1,
                        function(...) paste(..., collapse = "_"))
  expect_equal(out, expected_out)
})

test_that(".base_S3class_estfun fails with invalid base S3 class", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification = spec)
  mod@.S3Class <- "invalid_class"

  expect_error(.base_S3class_estfun(mod), "must have been fitted")
})

test_that("checking proper errors in conversion from lm to teeMod", {
  data(simdata)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)

  expect_error(as.lmitt(lm(y ~ x, data = simdata), specification = spec),
               "aliases are not found")


  # Take in either StudySpecification or WeightedStudySpecification
  mod1 <- propertee:::.convert_to_lmitt(lm(y ~ assigned(spec), data = simdata),
                                     spec,
                                     FALSE,
                                     TRUE,
                                     "a",
                                     call("quote", call("ls")))
  mod2 <- propertee:::.convert_to_lmitt(lm(y ~ assigned(spec), data = simdata),
                                     ate(spec, data = simdata),
                                     FALSE,
                                     TRUE,
                                     "a",
                                     call("quote", call("ls")))
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)

  expect_error(propertee:::.convert_to_lmitt(lm(y ~ assigned(spec), data = simdata),
                                     1,
                                     FALSE,
                                     TRUE,
                                     "a",
                                     call("quote", call("ls"))), "must be a")


})


test_that("disable update.teeMod", {
  data(simdata)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)

  mod <- lmitt(y ~ 1, data = simdata, specification = spec)

  expect_error(update(mod),
               "teeMod objects do not support")
})

test_that("lmitt_call", {
  data(simdata)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)

  # Make sure slot is actual call
  mod <- lmitt(y ~ 1, data = simdata, specification = spec)
  expect_true(is.call(mod@lmitt_call))

  # Call via `lmitt()` should match
  lmittcall <- str2lang("lmitt(y ~ 1, data = simdata, specification = spec)")
  mod <- eval(lmittcall)
  expect_equal(lmittcall, mod@lmitt_call)

  # Call via `lmitt.formula()` should match
  lmittformcall <- str2lang("lmitt.formula(y ~ 1, data = simdata, specification = spec)")
  mod <- eval(lmittformcall)
  expect_equal(lmittformcall, mod@lmitt_call)

  # as.lmitt, make sure is actual call
  lmmod <- lm(y ~ a.(spec), data = simdata)
  mod <- as.lmitt(lmmod, specification = spec)
  expect_true(is.call(mod@lmitt_call))

  # Call via `as.lmitt()` should match
  aslmittcall <- str2lang("as.lmitt(lmmod, specification = spec)")
  mod <- eval(aslmittcall)
  expect_equal(aslmittcall, mod@lmitt_call)

  # Call via `lmitt(lm` should match
  lmittlmcall <- str2lang("lmitt(lmmod, specification = spec)")
  mod <- eval(lmittlmcall)
  expect_equal(lmittlmcall, mod@lmitt_call)

  # Call via `lmitt.lm` should match
  lmittlmcall <- str2lang("lmitt.lm(lmmod, specification = spec)")
  mod <- eval(lmittlmcall)
  expect_equal(lmittlmcall, mod@lmitt_call)

})

test_that("lmitt_fitted object", {
  data(simdata)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)

  mod <- lmitt(y ~ as.factor(o), data = simdata, specification = spec)
  expect_true(mod@lmitt_fitted)

  mod <- as.lmitt(lm(y ~ adopters(spec), data = simdata), specification = spec)
  expect_false(mod@lmitt_fitted)

  mod <- lmitt(lm(y ~ adopters(spec), data = simdata), specification = spec)
  expect_false(mod@lmitt_fitted)

})

test_that("printed effects aren't confused by bad naming", {
  data(simdata)
  simdata$abz.c <- as.factor(simdata$o)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)

  mod <- lmitt(y ~ abz.c, data = simdata, specification = spec)
  co <- capture.output(show(mod))
  cos <- strsplit(trimws(co), "[[:space:]]+")

  expect_true(all(vapply(cos, length, numeric(1)) == 4))
  expect_true(all(!isTRUE(grepl("^abz", cos[[1]]))))
  expect_true(all(grepl("^z", cos[[1]])))

  # to force ` in variable names via as.factor
  mod <- lmitt(y ~ as.factor(abz.c), data = simdata, specification = spec)
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

  min_spec <- rct_spec(adsy ~ unitid(schoolid), min_df)

  res <- mapply(function(fmla) {
    as.lmitt(lm(fmla, min_df, weights = ett(min_spec)))
  },
  list(response_col ~ z.(), response_col ~ grade + grade:z.()),
  SIMPLIFY = FALSE)


  expect_true(all(sapply(res, is, "teeMod")))
  expect_true(all(sapply(sapply(sapply(res,
                                       slot, "lmitt_call"),
                                "[[", 1),
                         deparse) == "as.lmitt"))

  res <- lapply(list(response_col ~ z.(), response_col ~ grade + grade:z.()),
                function(fmla) {
                  as.lmitt(lm(fmla, min_df, weights = ett(min_spec)))
                })

  expect_true(all(sapply(res, is, "teeMod")))
  expect_true(all(sapply(sapply(sapply(res,
                                       slot, "lmitt_call"),
                                "[[", 1),
                         deparse) == "as.lmitt"))

})

test_that("#81 continuous moderator shows appropriate coefficients", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) , data = simdata)

  mod <- lmitt(y ~ x, data = simdata, specification = spec)
  expect_length(mod$coeff, 4)
  coefnames <- strsplit(capture.output(show(mod))[1], " +")[[1]]
  coefnames <- coefnames[nchar(coefnames) > 0]
  expect_true(all(grepl("z.", coefnames, fixed = TRUE)))
})

test_that("Invalid input to .convert_to_lmitt", {
  data(simdata)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)

  expect_error(.convert_to_lmitt(1, spec, FALSE, TRUE, "a",
                                 call("quote", call("ls"))))
})

test_that(".estfun_DB_blockabsorb returns 0 if not asking for
          specification-based SE or tee model does not absorb intercepts", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  ssmod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec))
  ssmod_abs <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                     absorb = TRUE)

  expect_true(all(.estfun_DB_blockabsorb(ssmod) == 0))
  expect_true(all(.estfun_DB_blockabsorb(ssmod, db = TRUE) == 0))
  expect_true(all(.estfun_DB_blockabsorb(ssmod_abs) == 0))
})

test_that(".estfun_DB_blockabsorb returns correct value", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  ssmod_abs <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                     absorb = TRUE)

  expect_false(all(.estfun_DB_blockabsorb(ssmod_abs, db = TRUE) == 0))

  phi <- .get_phi_tilde(ssmod_abs, db = TRUE)
  aa <- .get_appinv_atp(ssmod_abs, db = TRUE)
  expect_true(all.equal(
    cbind(0, phi %*% aa),
    .estfun_DB_blockabsorb(ssmod_abs, db = TRUE)
  ))
})
