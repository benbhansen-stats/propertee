save_options <- options()
options("propertee_message_on_unused_blocks" = FALSE)

test_that("lmitt", {

  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  mod <- lmitt(y ~ 1, weights = ate(), data = simdata, specification = spec)

  expect_s4_class(mod, "teeMod")
  expect_true(inherits(mod, "lm"))

  mod_ett <- lmitt(y ~ 1, weights = ett(), data = simdata, specification = spec)

  expect_s4_class(mod_ett, "teeMod")
})

test_that("lmitt and lm return the same in simple cases", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid),
                   data = simdata)

  ml <- lm(y ~ assigned(), data = simdata, weights = ate(spec))
  ssmod <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)
  ml2ssmod <- lmitt(ml)

  expect_true(all.equal(ssmod$coef[1:2], ml$coef, check.attributes = FALSE))
  expect_true(all.equal(ssmod$coef, ml2ssmod$coef, check.attributes = FALSE))
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(ssmod$model$`(weights)`, slot)),
                           deparse1(methods::slot(ml$model$`(weights)`, slot)))
               } else {
                 identical(methods::slot(ssmod$model$`(weights)`, slot),
                           methods::slot(ml$model$`(weights)`, slot))
               },
               logical(1L)))
  )
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(ssmod$model$`(weights)`, slot)),
                           deparse1(methods::slot(ml2ssmod$model$`(weights)`, slot)))
               } else {
                 identical(methods::slot(ssmod$model$`(weights)`, slot),
                           methods::slot(ml2ssmod$model$`(weights)`, slot))
               },
               logical(1L)))
  )

  ml <- lm(y ~ z, data = simdata, weights = ett(spec))
  ssmod <- lmitt(y ~ 1, data = simdata, weights = ett(), specification = spec)

  expect_true(all.equal(ssmod$coef[1:2], ml$coef, check.attributes = FALSE))
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(ssmod$model$`(weights)`, slot)),
                           deparse1(methods::slot(ml$model$`(weights)`, slot)))
               } else {
                 identical(methods::slot(ssmod$model$`(weights)`, slot),
                           methods::slot(ml$model$`(weights)`, slot))
               },
               logical(1L)))
  )

})

test_that("covariate adjustment", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  camod <- lm(y ~ x, data = simdata)

  ssmod <- lmitt(y ~ 1, data = simdata, weights = ate(),
              offset = cov_adj(camod), specification = spec)
  expect_true(!is.null(ssmod$model$"(offset)"))
  expect_true(inherits(ssmod$model$"(offset)", "SandwichLayer"))

})


test_that("StudySpecification argument", {
  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)
  mod3 <- lmitt(y ~ 1, data = simdata, specification = z ~ cluster(uoa1, uoa2))
  expect_true(mod3@StudySpecification@type == "Obs")
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_identical(mod1@StudySpecification@structure, mod3@StudySpecification@structure)
  expect_identical(mod1@StudySpecification@unit_of_assignment_type, mod3@StudySpecification@unit_of_assignment_type)
  expect_identical(deparse1(mod1@StudySpecification@call), deparse1(mod3@StudySpecification@call))

  spec2 <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec2)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec2)
  mod3 <- lmitt(y ~ 1, data = simdata,
                specification = z ~ cluster(uoa1, uoa2) + forcing(force))
  expect_true(mod3@StudySpecification@type == "RD")
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_identical(mod1@StudySpecification@structure, mod3@StudySpecification@structure)
  expect_identical(mod1@StudySpecification@unit_of_assignment_type, mod3@StudySpecification@unit_of_assignment_type)
  expect_identical(deparse1(mod1@StudySpecification@call), deparse1(mod3@StudySpecification@call))

  spec3 <- obs_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec3,
                subset = simdata$dose < 300)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(),
                subset = simdata$dose < 300, specification = spec3)
  mod3 <- lmitt(y ~ 1, data = simdata,
                specification = z ~ cluster(uoa1, uoa2) + block(bid),
                subset = simdata$dose < 300)
  expect_true(mod3@StudySpecification@type == "Obs")
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_identical(mod1@StudySpecification@structure, mod3@StudySpecification@structure)
  expect_identical(mod1@StudySpecification@unit_of_assignment_type, mod3@StudySpecification@unit_of_assignment_type)
  expect_identical(deparse1(mod1@StudySpecification@call), deparse1(mod3@StudySpecification@call))

  expect_error(lmitt(y ~ 1, data = simdata, specification = 4),
               "formula specifying such a specification")

})

test_that("Dichotomy argument", {

  data(simdata)
  spec <- obs_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec, dichotomy = dose > 200 ~ .)
  mod2 <- lmitt(y ~ 1, data = simdata,
                specification = dose ~ cluster(uoa1, uoa2),
                dichotomy = dose > 200 ~ .)
  expect_identical(mod1@StudySpecification@structure, mod2@StudySpecification@structure)
  expect_identical(mod1@StudySpecification@unit_of_assignment_type, mod2@StudySpecification@unit_of_assignment_type)
  expect_identical(deparse1(mod1@StudySpecification@call), deparse1(mod2@StudySpecification@call))
  expect_true(all.equal(mod1$coefficients, mod2$coefficients,
                        check.attributes =  FALSE))

})

test_that("Dichotomy from weights", {
  data(simdata)
  spec <- obs_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec,
                weights = ate(dichotomy = dose > 200 ~ .))
  mod2 <- lmitt(y ~ 1, data = simdata, specification = spec, weights = "ate",
                dichotomy = dose > 200 ~ .)
  mod3 <- lmitt(y ~ 1, data = simdata, specification = spec, weights = ate(spec),
                dichotomy = dose > 200 ~ .)
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_identical(mod1@StudySpecification, mod3@StudySpecification)
  expect_true(all.equal(mod1$coefficients, mod2$coefficients,
                        check.attributes =  FALSE))
  expect_true(all.equal(mod1$coefficients, mod3$coefficients,
                        check.attributes =  FALSE))
  expect_true(
    all.equal(mod2$model$`(weights)`, ate(spec, data = simdata, dichotomy = dose > 200 ~ .))
  )
  expect_true(
    all.equal(mod3$model$`(weights)`, ate(spec, data = simdata, dichotomy = dose > 200 ~ .))
  )
})

test_that("Minimal support for continuous treatments", {

  data(simdata)
  spec1 <- obs_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  expect_error(mod1 <- lmitt(y ~ 1, data = simdata, specification = spec1, absorb=FALSE),
               "continuous treatment")
  spec2 <- obs_spec(dose ~ cluster(uoa1, uoa2) +block(bid), data = simdata)
  expect_error(mod2 <- lmitt(y ~ 1, data = simdata, specification = spec2, absorb=TRUE),
               "continuous treatment")

  ## But we don't also allow factor or ordinal treatments
  simdata_ <-
      transform(simdata, dosef=cut(dose, breaks=c(50, 100, 200, 300),
                          labels=c("lo", "med", "hi"),
                          include.lowest=TRUE)
                )
  spec3 <-
      obs_spec(dosef ~ cluster(uoa1, uoa2), data = simdata_)
  expect_error(lmitt(y ~ 1, data = simdata, specification = spec3),
               "not supported")
})

test_that("control means for intercept only", {
  set.seed(93)
  pairs <- data.frame(b = c(rep(seq_len(2), each = 3), rep(3, 2)),
                      new_trt = rep(c(1, 0), 4),
                      x = rnorm(8),
                      dv = rnorm(8),
                      id = seq_len(8))
  spec <- rct_spec(new_trt ~ block(b) + unitid(id), pairs)
  
  covadj_data <- data.frame(dv = rnorm(30), x = rnorm(30), id = NA, b = NA)
  cmod <- lm(dv ~ x, covadj_data)
  cad <- cov_adj(cmod, newdata = pairs, spec = spec)
  
  # lmitt
  suppressMessages(m1 <- lmitt(dv ~ 1, spec, pairs))
  suppressMessages(m2 <- lmitt(lm(dv ~ z.(spec, pairs), pairs), spec))
  expect_equal(length(m1$coef), 3)
  expect_equal(m1$coef[3], c("dv:(Intercept)" = mean(pairs$dv[pairs$new_trt == 0])))
  expect_equal(m2$coef[3], m1$coef[3])

  # cov adj
  suppressMessages(m3 <- lmitt(dv ~ 1, spec, pairs, offset = cad))
  suppressMessages(m4 <- lmitt(lm(dv ~ a.(spec, pairs), pairs, offset = cad), spec))
  expect_equal(length(m3$coef), 4)
  expect_equal(m3$coef[3:4],
               c("dv:(Intercept)" = mean(pairs$dv[pairs$new_trt == 0]),
                 "cov_adj:(Intercept)" = mean(cad[pairs$new_trt == 0])))
  expect_equal(m3$coef[3:4], m4$coef[3:4])
  
  # case weights
  cw <- rep(c(10, 12), each = 4) 
  suppressMessages(m5 <- lmitt(dv ~ 1, spec, pairs, w = cw, offset = cad))
  suppressMessages(m6 <- lmitt(lm(dv ~ adopters(spec, pairs), pairs, w = cw, off = cad), spec))
  expect_equal(m5$coef[3:4],
               c("dv:(Intercept)" = sum(cw * pairs$dv * (1-pairs$new_trt)) / sum(cw * (1-pairs$new_trt)),
                 "cov_adj:(Intercept)" = sum(cw * cad * (1-pairs$new_trt)) / sum(cw * (1-pairs$new_trt))))
  expect_equal(m5$coef[3:4], m6$coef[3:4])
  
  # absorb = TRUE
  m7 <- lmitt(dv ~ 1, spec, pairs, offset = cad, absorb = TRUE)
  pis <- c(rep(2/3, 3), rep(1/3, 3), rep(1/2, 2))
  expect_equal(m7$coef[3:4],
               c("dv:(Intercept)" = sum(pis * pairs$dv * (1-pairs$new_trt)) / sum(pis * (1-pairs$new_trt)),
                 "cov_adj:(Intercept)" = sum(pis * cad * (1-pairs$new_trt)) / sum(pis * (1-pairs$new_trt))))
  
  # case weights and absorb = TRUE
  m8 <- lmitt(dv ~ 1, spec, pairs, w = cw, offset = cad, absorb = TRUE)
  cw_pis <- cw * (rowsum(cw * pairs$new_trt, pairs$b) / rowsum(cw, pairs$b))[pairs$b]
  expect_equal(m8$coef[3:4],
               c("dv:(Intercept)" = sum(cw_pis * pairs$dv * (1-pairs$new_trt)) / sum(cw_pis * (1-pairs$new_trt)),
                 "cov_adj:(Intercept)" = sum(cw_pis * cad * (1-pairs$new_trt)) / sum(cw_pis * (1-pairs$new_trt))))

  # missing data--stratum where trt is missing outcome still contributes to ctrl mean estimation
  pairs <- data.frame(b = c(rep(seq_len(2), each = 3), rep(3:4, each = 2)),
                      new_trt = rep(c(1, 0), 5),
                      x = rnorm(10),
                      dv = c(rnorm(8), NA_real_, rnorm(1)),
                      id = seq_len(10))
  spec <- rct_spec(new_trt ~ block(b) + unitid(id), pairs)
  m9 <- lmitt(dv ~ 1, spec, pairs, absorb = TRUE)
  pis <- c(pis, rep(1/2, 2))
  expect_equal(
    m9$coef[3],
    c("dv:(Intercept)" = sum(pis * pairs$dv * (1-pairs$new_trt), na.rm = TRUE) /
        sum(pis * (1-pairs$new_trt), na.rm = TRUE)))
})

test_that("control means for continuous moderator", {
  set.seed(93)
  pairs <- data.frame(b = c(rep(seq_len(2), each = 3), rep(3, 2)),
                      new_trt = rep(c(1, 0), 4),
                      x = rnorm(8),
                      mvar = sample(3, 8, replace = TRUE),
                      dv = rnorm(8),
                      id = seq_len(8),
                      pis = c(rep(2/3, 3), rep(1/3, 3), rep(1/2, 2)),
                      cw = rep(c(10, 12), each = 4))
  spec <- rct_spec(new_trt ~ block(b) + unitid(id), pairs)
  
  covadj_data <- data.frame(dv = rnorm(30), x = rnorm(30), id = NA, b = NA)
  cmod <- lm(dv ~ x, covadj_data)
  cad <- cov_adj(cmod, newdata = pairs, spec = spec)
  
  unwtd.ctrl.reg <- lm(cbind(dv, cad) ~ mvar, pairs, w = 1 - new_trt)
  cw.ctrl.reg <- lm(cbind(dv, cad) ~ mvar, pairs, w = cw * (1 - new_trt))
  pis.ctrl.reg <- lm(cbind(dv, cad) ~ mvar, pairs, w = pis * (1 - new_trt))
  cw_pis <- (rowsum(pairs$cw * pairs$new_trt, pairs$b) / rowsum(pairs$cw, pairs$b))[pairs$b]
  cw_pis.ctrl.reg <- lm(cbind(dv, cad) ~ mvar, pairs, w = pairs$cw * cw_pis * (1 - new_trt))
  
  # no weights
  suppressMessages(m1 <- lmitt(dv ~ mvar, spec, pairs))
  expect_equal(length(m1$coef), 6)
  expect_equal(
    m1$coef[5:6],
    setNames(unwtd.ctrl.reg$coef[,1], c("dv:(Intercept)", "dv:mvar"))
  )
  
  # cov adj
  suppressMessages(m2 <- lmitt(dv ~ mvar, spec, pairs, offset = cad))
  expect_equal(length(m2$coef), 8)
  expect_equal(
    m2$coef[5:8],
    setNames(c(unwtd.ctrl.reg$coef),
             c("dv:(Intercept)", "dv:mvar", "cov_adj:(Intercept)", "cov_adj:mvar"))
  )
  
  # case weights
  suppressMessages(m3 <- lmitt(dv ~ mvar, spec, pairs, weights = pairs$cw, offset = cad))
  expect_equal(length(m3$coef), 8)
  expect_equal(
    m3$coef[5:8],
    setNames(c(cw.ctrl.reg$coef),
             c("dv:(Intercept)", "dv:mvar", "cov_adj:(Intercept)", "cov_adj:mvar"))
  )
  
  # absorb = TRUE
  m4 <- lmitt(dv ~ mvar, spec, pairs, offset = cad, absorb = TRUE)
  expect_equal(length(m4$coef), 8)
  expect_equal(
    m4$coef[5:8],
    setNames(c(pis.ctrl.reg$coef),
             c("dv:(Intercept)", "dv:mvar", "cov_adj:(Intercept)", "cov_adj:mvar"))
  )
  
  # case weights and absorb = TRUE
  m5 <- lmitt(dv ~ mvar, spec, pairs, weights = pairs$cw, offset = cad, absorb = TRUE)
  expect_equal(length(m5$coef), 8)
  expect_equal(
    m5$coef[5:8],
    setNames(c(cw_pis.ctrl.reg$coef),
             c("dv:(Intercept)", "dv:mvar", "cov_adj:(Intercept)", "cov_adj:mvar"))
  )
})

test_that("control means for categorical moderator", {
  set.seed(93)
  # pairs <- data.frame(b = c(rep(seq_len(2), each = 3), rep(3, 2)),
  #                     new_trt = rep(c(1, 0), 4),
  #                     x = rnorm(8),
  #                     mvar = letters[rep(seq_len(3), 3)[-9]],
  #                     dv = rnorm(8),
  #                     id = seq_len(8),
  #                     pis = c(rep(2/3, 3), rep(1/3, 3), rep(1/2, 2)),
  #                     cw = rep(c(10, 12), each = 4))
  pairs <- data.frame(b = c(rep(seq_len(2), each = 7), rep(3, 10)),
                      new_trt = rep(c(1, 0), 12),
                      x = rnorm(24),
                      mvar = c(rep("a", 3), rep("b", 4), rep("a", 7), rep(c("a", "b"), each = 5)),
                      dv = rnorm(24),
                      id = seq_len(24),
                      pis = c(rep(2/3, 3), rep(1/2, 4), rep(3/7, 7), rep(3/5, 5), rep(2/5, 5)),
                      cw = rep(c(10, 12), each = 12))
  spec <- rct_spec(new_trt ~ block(b) + unitid(id), pairs)
  
  covadj_data <- data.frame(dv = rnorm(30), x = rnorm(30), id = NA, b = NA)
  cmod <- lm(dv ~ x, covadj_data)
  cad <- cov_adj(cmod, newdata = pairs, spec = spec)
  
  unwtd.ctrl.reg <- lm(cbind(dv, cad) ~ 0+mvar, pairs, w = 1 - new_trt)
  cw.ctrl.reg <- lm(cbind(dv, cad) ~ 0+mvar, pairs, w = cw * (1 - new_trt))
  pis.ctrl.reg <- lm(cbind(dv, cad) ~ 0+mvar, pairs, w = pis * (1 - new_trt))
  cw_pis <- (rowsum(pairs$cw * pairs$new_trt, paste(pairs$b, pairs$mvar, sep = "_")) /
               rowsum(pairs$cw, paste(pairs$b, pairs$mvar, sep = "_")))[paste(pairs$b, pairs$mvar, sep = "_"),]
  cw_pis.ctrl.reg <- lm(cbind(dv, cad) ~ 0+mvar, pairs, w = pairs$cw * cw_pis * (1 - new_trt))
  
  # no weights
  suppressMessages(m1 <- lmitt(dv ~ mvar, spec, pairs))
  expect_equal(length(m1$coef), 6)
  expect_equal(
    m1$coef[5:6],
    setNames(unwtd.ctrl.reg$coef[,1], c("dv:mvara", "dv:mvarb"))
  )
  
  # cov adj
  suppressMessages(m2 <- lmitt(dv ~ mvar, spec, pairs, offset = cad))
  expect_equal(length(m2$coef), 8)
  expect_equal(
    m2$coef[5:8],
    setNames(c(unwtd.ctrl.reg$coef),
             c("dv:mvara", "dv:mvarb", "cov_adj:mvara", "cov_adj:mvarb"))
  )
  
  # case weights
  suppressMessages(m3 <- lmitt(dv ~ mvar, spec, pairs, weights = pairs$cw, offset = cad))
  expect_equal(length(m3$coef), 8)
  expect_equal(
    m3$coef[5:8],
    setNames(c(cw.ctrl.reg$coef),
             c("dv:mvara", "dv:mvarb", "cov_adj:mvara", "cov_adj:mvarb"))
  )
  
  # absorb = TRUE
  suppressMessages(m4 <- lmitt(dv ~ mvar, spec, pairs, offset = cad, absorb = TRUE))
  expect_equal(length(m4$coef), 8)
  expect_equal(
    m4$coef[5:8],
    setNames(c(pis.ctrl.reg$coef),
             c("dv:mvara", "dv:mvarb", "cov_adj:mvara", "cov_adj:mvarb"))
  )
  
  # case weights and absorb = TRUE
  suppressMessages(m5 <- lmitt(dv ~ mvar, spec, pairs, w = pairs$cw, offset = cad, absorb = TRUE))
  expect_equal(length(m5$coef), 8)
  expect_equal(
    m5$coef[5:8],
    setNames(c(cw_pis.ctrl.reg$coef),
             c("dv:mvara", "dv:mvarb", "cov_adj:mvara", "cov_adj:mvarb"))
  )
})

test_that("control means for dichotomized treatment", {
  set.seed(93)
  strata <- data.frame(b = rep(seq_len(2), each = 10),
                       trt2dich = letters[rep(rep(seq_len(3), 4)[1:10], 2)],
                       x = rnorm(20),
                       dv = rnorm(20),
                       id = seq_len(20))
  spec <- rct_spec(trt2dich ~ block(b) + unitid(id), strata)
  
  covadj_data <- data.frame(dv = rnorm(30), x = rnorm(30), id = NA, b = NA)
  cmod <- lm(dv ~ x, covadj_data)
  expect_warning(cad <- cov_adj(cmod, newdata = strata, spec = spec),
                 "treatment is an independent variable")
  
  suppressMessages(m1 <- lmitt(dv ~ 1, spec, strata, offset = cad,
                               dichotomy = trt2dich %in% c("a", "b") ~ trt2dich == "c"))
  expect_equal(length(m1$coef), 4)
  expect_equal(
    m1$coef[3:4],
    c("dv:(Intercept)" = mean(strata$dv[strata$trt2dich == "c"]),
      "cov_adj:(Intercept)" = mean(cad[strata$trt2dich == "c"]))
  )
  
  m2 <- lmitt(dv ~ 1, spec, strata, offset = cad, absorb = TRUE,
              dichotomy = trt2dich %in% c("a", "b") ~ trt2dich == "c")
  pis <- tapply(strata$trt2dich %in% c("a", "b"), strata$b, mean)[strata$b]
  expect_equal(length(m2$coef), 4)
  expect_equal(
    m2$coef[3:4],
    c("dv:(Intercept)" = sum(pis * strata$dv * (strata$trt2dich == "c")) / sum(pis * (strata$trt2dich == "c")),
      "cov_adj:(Intercept)" = sum(pis * cad * (strata$trt2dich == "c")) / sum(pis * (strata$trt2dich == "c")))
  )
})

test_that("lmitt finds StudySpecification wherever it's stored", {
  data(simdata)
  spec_form <- z ~ uoa(uoa1, uoa2) + block(bid)
  spec <- obs_spec(spec_form, data = simdata)
  camod <- lm(y ~ x, data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, weights = ate(),
                offset = cov_adj(camod, specification = spec),
                specification = spec)
  mod2 <- lmitt(y ~ 1, specification = spec_form, data = simdata, weights = ate(),
                offset = cov_adj(camod))
  mod3 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(),
                offset = cov_adj(camod))

  expect_equal(coef(mod1), coef(mod2))
  expect_equal(coef(mod1), coef(mod3))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1), vcov_tee(mod3))
})

test_that("Allowed inputs to lmitt #73", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  ### Allowed versions

  expect_no_error(l1 <- lmitt(y ~ 1, specification = spec, data = simdata))
  expect_no_error(l2 <- lmitt(y ~ dose, specification = spec, data = simdata))

  ### Disallowed versions

  expect_error(lmitt(y ~ 0, specification = spec, data = simdata))
  expect_error(lmitt(y ~ dose + 0, specification = spec, data = simdata))
  expect_error(lmitt(y ~ 0 + dose, specification = spec, data = simdata))
  expect_error(lmitt(y ~ dose + bid, specification = spec, data = simdata))
  expect_error(lmitt(y ~ 1 + dose + bid, specification = spec, data = simdata))
  expect_error(lmitt(y ~ 1 + x, specification = spec, data = simdata))
  expect_error(lmitt(y ~ x + 1, specification = spec, data = simdata))

  ### Main + Subgroup effects.
  ### NYI - see #73
  expect_error(lmitt(y ~ 1 + dose, specification = spec, data = simdata),
               "To estimate subgroup effects") # Remove this once #73 is fixed

  #expect_no_error(l3 <- lmitt(y ~ 1 + dose, specification = spec, data = simdata))
  #expect_no_error(l4 <- lmitt(y ~ dose + 1, specification = spec, data = simdata))

  #expect_identical(l3$coeff, l4$coeff)

})

test_that("weights argument can be string", {

  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  l1 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate())
  l2 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate")

  l3 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ett())
  l4 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ett")

  l1@lmitt_call <- call("ls")
  l2@lmitt_call <- call("ls")
  l3@lmitt_call <- call("ls")
  l4@lmitt_call <- call("ls")

  expect_true(all.equal(l1, l2))
  expect_true(all.equal(l3, l4))

  # all.equal returns strings detailing differences
  expect_true(is.character(all.equal(l1, l3)))

  # Passing a different string warns users, but allows `lm()` to error
  # in case user is trying to pass something to someplace else
  expect_error(expect_warning(lmitt(y ~ 1, specification = spec, data =simdata,
                                    weights = "abc"), "other than"))

})

test_that("Regular weights can still be used", {

  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ 1, specification = spec, weights = simdata$dose, data = simdata)
  expect_true(is(mod, "teeMod"))

})

test_that("Aliases for assigned() aren't allowed either", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_error(lmitt(y ~ assigned(), specification = spec, data = simdata),
               "Do not specify")

  expect_error(lmitt(y ~ adopters(), specification = spec, data = simdata),
               "Do not specify")

  expect_error(lmitt(y ~ a.(), specification = spec, data = simdata),
               "Do not specify")


  expect_error(lmitt(y ~ z.(), specification = spec, data = simdata),
               "Do not specify")

  #### because of the `.` in `a.()`, let's make sure the regexp isn't breaking
  ad <- assigned
  # Passing in a different weight to ensure we get a predictable error that
  # ocurs *after* checking the formula
  spec2 <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  expect_error(lmitt(y ~ ad(), specification = spec, data = simdata,
                     weights = ate(spec2)),
               "Multiple differing")

  # Just to ensure that `lmitt.formula()` doesn't get changed in a way such that
  # checking for differing specifications happens prior:
  expect_error(lmitt(y ~ a.(), specification = spec, data = simdata,
                     weights = ate(spec2)),
               "Do not specify")

  rm(ad)

})

test_that("User can pass WeightedStudySpecification", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ 1, data = simdata, specification = ate(spec))

  expect_identical(spec, mod@StudySpecification)

})

test_that("#82 Informative error if specification in weight/offset but not lmitt", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_error(lmitt(y ~ 1, data = simdata, weights = ate(spec)),
               "into the weight function")

  cmod <- lm(y ~ x, data = simdata)
  expect_error(lmitt(y ~ 1, data = simdata, offset = cov_adj(cmod, specification = spec)),
               "into the offset function")

})

test_that("#81 continuous moderator", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) , data = simdata)

  mod <- lmitt(y ~ x, data = simdata, specification = spec)
  expect_length(mod$coeff, 6)

  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ x, data = simdata, specification = spec, absorb = TRUE)
  expect_length(mod$coeff, 6)

})

test_that("non-integer units of assignment", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, specification = spec, absorb = TRUE))

  simdata$uoa1 <- as.character(simdata$uoa1)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, specification = spec, absorb = TRUE))

  simdata$uoa1 <- as.factor(simdata$uoa1)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, specification = spec, absorb = TRUE))

  expect_true(TRUE) # Avoid an empty test warning

})

test_that("#137 ensure absorb is using the correct moderator", {
  data(simdata)

  mod1 <- lmitt(y ~ o, data = simdata,
                specification = z ~ cluster(uoa1, uoa2) + block(bid))
  mod2 <- lmitt(y ~ o, data = simdata,
                specification = z ~ cluster(uoa1, uoa2) + block(bid),
                absorb = TRUE)

  o1 <- model.matrix(mod1)[, "o"]
  o2 <- model.matrix(mod2)[, "o"]

  # o1 and o2 should be different, since the latter should have been
  # group-centered due to `absorb = TRUE`. Due to a bug discovered in #137, the
  # un-centered version of the continuous moderator was being using previously,
  # leading o1 and o2 to be identical (incorrectly).

  expect_true(!(all.equal(o1, o2, check.attributes = FALSE) == TRUE))


})

test_that("#140 handling 0 weights", {

  # A single 0 weight observation in a larger block - that observation's
  # conbribution to the estimating equation should be 0
  data(simdata)
  simdata$weight <- 1
  simdata$weight[1] <- 0

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  mod <- lmitt(y ~ x, data = simdata, specification = spec,
               weights = simdata$weight)
  ee <- propertee:::estfun.teeMod(mod)
  expect_true(all(ee[1, ] == 0))

  # A block with only a single non-zero weight should contribute nothing to the
  # estimating equation
  data(simdata)
  simdata$weight <- 1
  simdata$weight[simdata$bid == 1] <- 0
  simdata$weight[simdata$bid == 1][1] <- 1

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  mod <- lmitt(y ~ x, data = simdata, specification = spec,
               absorb = TRUE, weights = simdata$weight)
  ee <- propertee:::estfun.teeMod(mod)
  expect_true(all(ee[simdata$bid == 1, 2:ncol(ee)] == 0))
  expect_true(all(ee[simdata$bid == 1 & simdata$weight == 0, 1] == 0))

  data(simdata)
  simdata$z[simdata$bid == 1] <- 1
  simdata$weight <- 1
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  # Block 1 has 0 variance in treatment
  mod <- lmitt(y ~ x, data = simdata, specification = spec,
               absorb = TRUE, weights = simdata$weight)
  ee <- estfun.teeMod(mod)
  expect_true(all(ee[simdata$bid == 1, "z."] == 0))
  expect_true(all(
    sapply(ee[simdata$bid == 1, c("(Intercept)", "x", "z._x")],
           function(vals) any(vals != 0))
  ))

})


options(save_options)
#### !!!!!!!!!!!NOTE!!!!!!!!!!!!!
# This test below should NOT have `options()$propertee_message_on_unused_blocks`
# set to FALSE. So it needs to stay below the restoration of options line above.
# Other tests should probably go above the restoration of options line.

test_that("Message if specification has block info but isn't used in lmitt", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_message(lmitt(y ~ 1, data = simdata, specification = spec),
                 "is not used in this model")
  expect_message(lmitt(y ~ 1, data = simdata, specification = ate(spec)),
                 "is not used in this model")
  expect_message(lmitt(y ~ 1, data = simdata, specification = ate(spec),
                       weights = abs(simdata$x)),
                 "is not used in this model")

  expect_silent(x <- lmitt(y ~ 1, data = simdata, absorb = TRUE, specification = spec))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = "ate", specification = spec))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = ett()*3, specification = spec))

})

### READ COMMENT ABOUT LAST TEST!

test_that(paste("absorb=TRUE doesn't drop rows when some strata have weights=0",
                "(some strata only have treated clusters)"), {
  data(simdata)
  simdata[simdata$uoa1 == 5 & simdata$uoa2 == 2, "bid"] <- 4

  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  absorb_mod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate", absorb = TRUE)
  no_absorb_mod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate", absorb = FALSE)

  expect_equal(dim(model.matrix(absorb_mod)), dim(model.matrix(no_absorb_mod)))
})

test_that(paste("absorb=TRUE doesn't drop rows when some strata have weights=0",
                "(some strata only have untreated clusters)"), {
  data(simdata)
  simdata[simdata$uoa1 == 4 & simdata$uoa2 == 2, "bid"] <- 4

  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  absorb_mod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate", absorb = TRUE)
  no_absorb_mod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate", absorb = FALSE)

  expect_equal(dim(model.matrix(absorb_mod)), dim(model.matrix(no_absorb_mod)))
})

test_that("#147 lmitt.formula accepts references to formula objects", {
  data(simdata)
  lmitt_form <- y ~ 1
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  expect_equal(
    co <- capture.output(lmitt(lmitt_form, data = simdata, specification = spec)),
    capture.output(lmitt(y ~ 1, data = simdata, specification = spec))
  )

  make_lmitt <- function(data, dv_col) {
    lmitt_form <- as.formula(paste0(dv_col, "~1"))
    spec <- rct_spec(z ~ unitid(uoa1, uoa2), data)
    lmitt(lmitt_form, specification = spec, data = data)
  }
  expect_equal(capture.output(make_lmitt(simdata, "y")), co)

})
