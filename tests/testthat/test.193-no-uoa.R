test_that("Support lack of UOA", {
  data(simdata)
  sd_idd <- simdata
  sd_idd$id <- rownames(sd_idd)

  spec1 <- rct_spec(z ~ uoa(id), data = sd_idd)
  spec2 <- rct_spec(z ~ 1, data = sd_idd)
  expect_true(all.equal(spec1@structure,
                        spec2@structure,
                        check.attributes = FALSE))
  expect_identical(spec1@column_index,
                   spec2@column_index)
  expect_equal(spec1@unit_of_assignment_type, "unit_of_assignment")
  expect_equal(spec2@unit_of_assignment_type, "none")

  mod1 <- lmitt(y ~ 1, spec = spec1, data = sd_idd)
  mod2 <- lmitt(y ~ 1, spec = spec2, data = sd_idd)
  expect_identical(coefficients(mod1),
                   coefficients(mod2))

  smod1 <- summary(mod1)
  smod2 <- summary(mod2)
  expect_identical(coefficients(smod2), coefficients(smod1))
})
