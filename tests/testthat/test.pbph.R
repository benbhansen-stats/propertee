test_that("Basics of pbph", {
  data(simdata)

  cmod <- lm(y ~ x, data = simdata, subset = (z == 0))

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)

  expect_warning(mod <- lmitt_pbph(cmod, des, simdata))

  expect_true(is(mod, "DirectAdjusted"))
  expect_true(is(mod, "DirectAdjustedPBPH"))

  expect_length(mod$coefficients, 2)

  expect_true(all(grepl("assigned()", names(mod$coefficients), fixed = TRUE)))

  smod <- summary(mod)

  expect_false(is(smod, "summary.DirectAdjusted"))
  expect_true(is(smod, "summary.DirectAdjustedPBPH"))


})
