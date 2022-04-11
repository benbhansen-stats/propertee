test_that("Treatment inputs", {
  data(simdata)

  # numeric 0/1 treatment
  d <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  expect_true(is.numeric(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))

  # numeric non-0/1
  d <- rct_design(o ~ cluster(cid1, cid2), data = simdata)
  expect_true(is.numeric(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))

  # factor 0/1
  simdata$foo <- as.factor(simdata$z)
  expect_warning(d <- rct_design(foo ~ cluster(cid1, cid2), data = simdata),
                 "conversion")
  expect_true(is.numeric(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))

  # factor numeric non 0/1
  simdata$foo <- as.factor(simdata$o)
  expect_warning(d <- rct_design(foo ~ cluster(cid1, cid2), data = simdata),
                 "conversion")
  expect_true(is.numeric(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))

  # factor string levels
  simdata$foo <- factor(letters[simdata$o])
  expect_warning(d <- rct_design(foo ~ cluster(cid1, cid2), data = simdata),
                 "conversion")
  expect_true(is.character(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))

  # ordinal 0/1
  simdata$foo <- as.ordered(simdata$z)
  expect_warning(d <- rct_design(foo ~ cluster(cid1, cid2), data = simdata),
                 "conversion")
  expect_true(is.numeric(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))

  # ordinal numeric non 0/1
  simdata$foo <- as.ordered(simdata$o)
  expect_warning(d <- rct_design(foo ~ cluster(cid1, cid2), data = simdata),
                 "conversion")
  expect_true(is.numeric(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))

  # ordinal string levels
  simdata$foo <- ordered(letters[simdata$o])
  expect_warning(d <- rct_design(foo ~ cluster(cid1, cid2), data = simdata),
                 "conversion")
  expect_true(is.character(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))

  # character
  simdata$foo <- letters[simdata$o]
  d <- rct_design(foo ~ cluster(cid1, cid2), data = simdata)
  expect_true(is.character(treatment(d)[, 1]))
  expect_false(has_binary_treatment(d))

  # logical
  simdata$foo <- as.logical(simdata$z)
  d <- rct_design(foo ~ cluster(cid1, cid2), data = simdata)
  expect_true(is.logical(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))

  # conditional
  expect_warning(d <- rct_design(o > 2 ~ uoa(cid1, cid2), data = simdata),
                 "conditional logic")
  expect_true(is.logical(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))

  old_opt <- options()
  on.exit(options(old_opt))
  options(flexida_warn_on_conditional_treatment = FALSE)
  d <- rct_design(o > 2 ~ uoa(cid1, cid2), data = simdata)
  expect_true(is.logical(treatment(d)[, 1]))
  expect_true(has_binary_treatment(d))
  options(old_opt)

  # tests for other things...?
})
