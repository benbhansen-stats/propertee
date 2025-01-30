test_that("Support lack of UOA", {
  data(simdata)
  simdata$id <- rownames(simdata)

  spec1 <- rct_spec(z ~ uoa(id), data = simdata)
  spec2 <- rct_spec(z ~ 1, data = simdata)
  expect_true(all.equal(spec1@structure,
                        spec2@structure,
                        check.attributes = FALSE))
  expect_identical(spec1@column_index,
                   spec2@column_index)
  expect_equal(spec1@unit_of_assignment_type, "unit_of_assignment")
  expect_equal(spec2@unit_of_assignment_type, "none")

  mod1 <- lmitt(y ~ 1, spec = spec1, data = simdata)
  mod2 <- lmitt(y ~ 1, spec = spec2, data = simdata)
  expect_identical(coefficients(mod1),
                   coefficients(mod2))

  summary(mod1)

})
