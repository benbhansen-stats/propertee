test_that("equivalent of lm and lmitt calls", {
  data(simdata)
  simdata$o <- as.factor(simdata$o)
  des <- obs_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)


  ### No weights

  # Just treatment
  mod1_lm <- lm(y ~ assigned(des), data = simdata)
  expect_equal(length(mod1_lm$coeff), 2)
  mod1_lmitt <- lmitt(y ~ 1, data = simdata, design = des)
  expect_equal(length(mod1_lmitt$coeff), 2)
  expect_true(all.equal(mod1_lm$coef[2], mod1_lmitt$coef[2],
                        check.attributes = FALSE))

  # Absorb = TRUE
  mod2_lm <- lm(y ~ assigned(des) + as.factor(bid), data = simdata)
  expect_equal(length(mod2_lm$coeff), 4)
  mod2_lmitt <- lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE)
  expect_equal(length(mod2_lmitt$coeff), 2)
  expect_true(all.equal(mod2_lm$coef[2], mod2_lmitt$coef[2],
                        check.attributes = FALSE))

  # Subgroup
  mod3_lm <- lm(y ~ assigned(des):o + o, data = simdata)
  expect_equal(length(mod3_lm$coeff), 8)
  mod3_lmitt <- lmitt(y ~ o, data = simdata, design = des)
  expect_equal(length(mod3_lmitt$coeff), 4)
  expect_true(all.equal(mod3_lm$coef[5:8], mod3_lmitt$coef,
                        check.attributes = FALSE))

  # Subgroup, absorb = TRUE
  mod4_lm <- lm(y ~ assigned(des):o + o  + as.factor(bid), data = simdata)
  expect_equal(length(mod4_lm$coeff), 10)
  mod4_lmitt <- lmitt(y ~ o, data = simdata, design = des, absorb = TRUE)
  expect_equal(length(mod4_lmitt$coeff), 4)
  expect_true(all.equal(mod4_lm$coef[7:10], mod4_lmitt$coef,
                        check.attributes = FALSE))


  ### Weights

  # Just treatment
  wmod1_lm <- lm(y ~ assigned(des), data = simdata, weights = ate(des))
  expect_equal(length(wmod1_lm$coeff), 2)
  wmod1_lmitt <- lmitt(y ~ 1, data = simdata, design = des, weights = "ate")
  expect_equal(length(wmod1_lmitt$coeff), 2)
  expect_true(all.equal(wmod1_lm$coef[2], wmod1_lmitt$coef[2],
                        check.attributes = FALSE))

  # Absorb = TRUE
  wmod2_lm <- lm(y ~ assigned(des) + as.factor(bid), data = simdata, weights = ett(des))
  expect_equal(length(wmod2_lm$coeff), 4)
  wmod2_lmitt <- lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE, weights = "ett")
  expect_equal(length(wmod2_lmitt$coeff), 2)
  expect_true(all.equal(wmod2_lm$coef[2], wmod2_lmitt$coef[2],
                        check.attributes = FALSE))

  # Subgroup
  wmod3_lm <- lm(y ~ assigned(des):o + o, data = simdata, weights = ett(des))
  expect_equal(length(wmod3_lm$coeff), 8)
  wmod3_lmitt <- lmitt(y ~ o, data = simdata, design = des, weights = "ett")
  expect_equal(length(wmod3_lmitt$coeff), 4)
  expect_true(all.equal(wmod3_lm$coef[5:8], wmod3_lmitt$coef,
                        check.attributes = FALSE))

  # Subgroup, absorb = TRUE
  wmod4_lm <- lm(y ~ assigned(des):o + o  + as.factor(bid), data = simdata, weights = ate(des))
  expect_equal(length(wmod4_lm$coeff), 10)
  wmod4_lmitt <- lmitt(y ~ o, data = simdata, design = des, absorb = TRUE, weights = "ate")
  expect_equal(length(wmod4_lmitt$coeff), 4)
  expect_true(all.equal(wmod4_lm$coef[7:10], wmod4_lmitt$coef,
                        check.attributes = FALSE))

})
