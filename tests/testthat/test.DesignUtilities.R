test_that("Design formula checking", {
  expect_true(.check_design_formula(y ~ cluster(x)))
  expect_true(.check_design_formula(y ~ cluster(x, z, q, r)))
  expect_true(.check_design_formula(y ~ unitid(x)))
  expect_true(.check_design_formula(y ~ unitid(x, z, q, r)))
  expect_true(.check_design_formula(y ~ unit_of_assignment(x)))
  expect_true(.check_design_formula(y ~ unit_of_assignment(x, z, q, r)))
  expect_true(.check_design_formula(y ~ uoa(x)))
  expect_true(.check_design_formula(y ~ uoa(x, z, q, r)))

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

  expect_error(.check_design_formula(y ~ unit_of_assignment(x) + unit_of_assignment(z)),
               "Only one instance of `unit_of")

  expect_error(.check_design_formula(y ~ uoa(x) + uoa(z)),
               "Only one instance of `unit_of")

  expect_true(.check_design_formula(y ~ cluster(x) + block(z)))
  expect_true(.check_design_formula(y ~ cluster(x) + block(z, a, b, c)))

  expect_error(.check_design_formula(y ~ cluster(x) + block(z) + block(q)),
               "only one block")

  expect_error(.check_design_formula(y ~ cluster(x) + forcing(z)),
               "only allowed")

  expect_true(.check_design_formula(y ~ cluster(x) + forcing(z),
                                 allow_forcing = TRUE))

  expect_error(.check_design_formula(y ~ cluster(x) + forcing(z) + forcing(q),
                                  allow_forcing = TRUE),
               "only one forcing")

})

test_that("binary treatment and dichotomization", {
  expect_error(is_dichotomized(1),
               "must be")
  expect_error(has_binary_treatment(1),
               "must be")

  des1 <- obs_design(z ~ uoa(cid1, cid2), data = simdata)
  des2 <- obs_design(o ~ uoa(cid1, cid2), data = simdata)
  des3 <- obs_design(o ~ uoa(cid1, cid2), data = simdata, dichotomize = o > 3 ~ o == 1)

  expect_false(is_dichotomized(des1))
  expect_true(has_binary_treatment(des1))
  expect_false(is_dichotomized(des2))
  expect_false(has_binary_treatment(des2))
  expect_true(is_dichotomized(des3))
  expect_true(has_binary_treatment(des3))

  # Adding afterwards
  dichotomization(des2) <- o >= 2 ~ .
  expect_true(is_dichotomized(des2))
  expect_true(has_binary_treatment(des2))

  # Removing
  dichotomization(des3) <- NULL
  expect_false(is_dichotomized(des3))
  expect_false(has_binary_treatment(des3))
})
