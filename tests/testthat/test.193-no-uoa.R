old_opt <- options()
on.exit(options(old_opt))
options(propertee_warn_on_no_unit_of_assignment = FALSE)

test_that("Support lack of UOA", {
  data(simdata)
  # I was having weird issues with testthat when I worked with simdata directly;
  # so that's why I'm making a copy here - JE
  sd_ided <- simdata
  sd_ided$id <- rownames(sd_ided)

  spec1 <- rct_spec(z ~ uoa(id), data = sd_ided)
  spec2 <- rct_spec(z ~ 1, data = sd_ided)
  expect_true(all.equal(spec1@structure,
                        spec2@structure,
                        check.attributes = FALSE))
  expect_identical(spec1@column_index,
                   spec2@column_index)
  expect_equal(spec1@unit_of_assignment_type, "unit_of_assignment")
  expect_equal(spec2@unit_of_assignment_type, "none")

  mod1 <- lmitt(y ~ 1, spec = spec1, data = sd_ided)
  mod2 <- lmitt(y ~ 1, spec = spec2, data = sd_ided)
  expect_true(isTRUE(all.equal(coefficients(mod1),
                               coefficients(mod2))))

  smod1 <- summary(mod1)
  smod2 <- summary(mod2)
  expect_true(isTRUE(all.equal(coefficients(smod1),
                               coefficients(smod2))))

  ## Shuffle
  sd_ided_shuffled <- sd_ided[sample(1:50), ]
  spec3 <- rct_spec(z ~ 1, data = sd_ided_shuffled)
  mod3 <- lmitt(y ~ 1, spec = spec3, data = sd_ided_shuffled)
  expect_true(isTRUE(all.equal(coefficients(mod1),
                               coefficients(mod3))))
  smod3 <- summary(mod3)
  expect_true(isTRUE(all.equal(coefficients(smod1),
                               coefficients(smod3))))


})

test_that("covadj", {
  data(simdata)
  sd_ided <- simdata
  sd_ided$id <- rownames(sd_ided)

  spec1 <- rct_spec(z ~ uoa(id), data = sd_ided)
  spec2 <- rct_spec(z ~ 1, data = sd_ided)

  camod <- lm(y ~ x, data = sd_ided, subset = z == 0)

  mod1 <- lmitt(y ~ 1, spec = spec1, offset = cov_adj(camod), data = sd_ided)
  mod2 <- lmitt(y ~ 1, spec = spec2, offset = cov_adj(camod), data = sd_ided)
  expect_true(isTRUE(all.equal(coefficients(mod1),
                               coefficients(mod2))))

  smod1 <- summary(mod1)
  smod2 <- summary(mod2)
  expect_true(isTRUE(all.equal(coefficients(smod1),
                               coefficients(smod2))))



  data(STARplus)
  st <- STARplus
  st$cond_small <- st$cond_at_entry == "small"
  spec1 <- obs_spec(cond_small ~ unit_of_assignment(stdntid) +
                      block(school_at_entry),
                   data = st, na.fail = FALSE)
  ca <- lm(read_yr1 ~ dob + gender, data = st, subset = !(cond_small))

  mod1 <- lmitt(read_yr1 ~ 1, spec1, data = st, absorb = TRUE,
                                 offset = cov_adj(ca))

  spec2 <- obs_spec(cond_small ~ block(school_at_entry),
                   data = st, na.fail = FALSE)
  suppressWarnings(mod2 <- lmitt(read_yr1 ~ 1, spec2, absorb = TRUE, data = st,
                                 offset = cov_adj(ca)))

  expect_true(isTRUE(all.equal(coefficients(mod1),
                               coefficients(mod2))))
  smod1 <- summary(mod1, vcov.type = "CR0")
  smod2 <- summary(mod2, vcov.type = "CR0")
  expect_true(isTRUE(all.equal(coefficients(smod1),
                               coefficients(smod2))))
})


options(old_opt)

test_that("warning message", {
  data(simdata)
  # Running each a few times to ensure options are reset appropriately
  expect_warning(rct_spec(z ~ 1, data = simdata))
  expect_no_warning(lmitt(y ~ 1, spec = z ~ 1, data = simdata))
  expect_warning(rct_spec(z ~ 1, data = simdata))
  expect_no_warning(lmitt(y ~ 1, spec = z ~ 1, data = simdata))
  expect_warning(rct_spec(z ~ 1, data = simdata))
})
