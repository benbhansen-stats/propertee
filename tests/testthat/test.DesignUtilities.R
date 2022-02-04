test_that("Design formula checking", {
  expect_true(.check_design_formula(y ~ cluster(x)))
  expect_true(.check_design_formula(y ~ cluster(x,z,q,r)))
  expect_true(.check_design_formula(y ~ unitid(x)))
  expect_true(.check_design_formula(y ~ unitid(x,z,q,r)))
  expect_true(.check_design_formula(y ~ unitOfAssignment(x)))
  expect_true(.check_design_formula(y ~ unitOfAssignment(x,z,q,r)))
  expect_true(.check_design_formula(y ~ uoa(x)))
  expect_true(.check_design_formula(y ~ uoa(x,z,q,r)))

  expect_error(.check_design_formula(~ cluster(x)),
               "treatment")

  expect_error(.check_design_formula(y ~ x),
               "cluster")

  expect_error(.check_design_formula(y ~ cluster(x) + unitid(z)),
               "Only one of")

  expect_error(.check_design_formula(y ~ cluster(x) + cluster(z)),
               "Only one instance of `cluster")

  expect_error(.check_design_formula(y ~ unitid(x) + unitid(z)),
               "Only one instance of `unitid")

  expect_error(.check_design_formula(y ~ unitOfAssignment(x) + unitOfAssignment(z)),
               "Only one instance of `unitOf")

  expect_error(.check_design_formula(y ~ uoa(x) + uoa(z)),
               "Only one instance of `unitOf")

  expect_true(.check_design_formula(y ~ cluster(x) + block(z)))
  expect_true(.check_design_formula(y ~ cluster(x) + block(z, a, b, c)))

  expect_error(.check_design_formula(y ~ cluster(x) + block(z) + block(q)),
               "only one block")

  expect_error(.check_design_formula(y ~ cluster(x) + forcing(z)),
               "only allowed")

  expect_true(.check_design_formula(y ~ cluster(x) + forcing(z),
                                 allowForcing = TRUE))

  expect_error(.check_design_formula(y ~ cluster(x) + forcing(z) + forcing(q),
                                  allowForcing = TRUE),
               "only one forcing")

})
