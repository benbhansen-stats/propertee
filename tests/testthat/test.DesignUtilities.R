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

  expect_error(.check_design_formula(y ~ unit_of_assignment(x) +
                                       unit_of_assignment(z)),
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

test_that("binary treatment and dichotomy", {
  expect_error(is_dichotomized(1),
               "must be")
  expect_error(has_binary_treatment(1),
               "must be")
  expect_error(is_binary_or_dichotomized(1),
               "must be")

  des1 <- obs_design(z ~ uoa(cid1, cid2), data = simdata)
  des2 <- obs_design(o ~ uoa(cid1, cid2), data = simdata)
  des3 <- obs_design(o ~ uoa(cid1, cid2), data = simdata,
                     dichotomy = o > 3 ~ o == 1)
  des4 <- obs_design(z ~ uoa(cid1, cid2), data = simdata,
                     dichotomy = z == 1 ~ z == 0)

  expect_false(is_dichotomized(des1))
  expect_true(has_binary_treatment(des1))
  expect_true(is_binary_or_dichotomized(des1))

  expect_false(is_dichotomized(des2))
  expect_false(has_binary_treatment(des2))
  expect_false(is_binary_or_dichotomized(des2))

  expect_true(is_dichotomized(des3))
  expect_false(has_binary_treatment(des3))
  expect_true(is_binary_or_dichotomized(des3))

  expect_true(is_dichotomized(des4))
  expect_true(has_binary_treatment(des4))
  expect_true(is_binary_or_dichotomized(des4))

  # Adding afterwards
  dichotomy(des2) <- o >= 2 ~ .
  expect_true(is_dichotomized(des2))
  expect_false(has_binary_treatment(des2))
  expect_true(is_binary_or_dichotomized(des2))

  # Removing
  dichotomy(des3) <- NULL
  expect_false(is_dichotomized(des3))
  expect_false(has_binary_treatment(des3))
  expect_false(is_binary_or_dichotomized(des3))
})

test_that("identical_Designs function", {
  data(simdata)

  des1 <- rct_design(dose ~ cluster(cid1, cid2), data = simdata)
  des2 <- rct_design(dose ~ cluster(cid1, cid2), data = simdata,
                     dichotomy = dose > 200 ~ dose <= 200)
  des3 <- rct_design(dose ~ cluster(cid1, cid2) + block(bid), data = simdata)

  expect_true(identical_Designs(des1, des2))
  expect_false(identical_Designs(des1, des3))
  expect_false(identical_Designs(des2, des3))

  expect_false(identical_Designs(des1, des2, TRUE))
  expect_false(identical_Designs(des1, des3, TRUE))
  expect_false(identical_Designs(des2, des3, TRUE))

  # Was hitting an error if the Design passed in via `design=` contained a
  # dichotomy. `des2` above has a dichotomy.
  expect_no_error(lmitt(y ~ 1, design = des2, data = simdata, weights = "ate"))

})

test_that("identify_small_blocks", {
  des_data <- data.frame(uid = letters[1:12],
                         bid = rep(LETTERS[1:3], each = 4),
                         a = c(rep(rep(c(0, 1), each = 2), 2), 0, rep(1, 3)))
  des <- rct_design(a ~ unitid(uid) + block(bid), des_data)
  
  small_blocks <- identify_small_blocks(des)
  expect_equal(length(small_blocks), 3)
  expect_equal(names(small_blocks), LETTERS[1:3])
  expect_true(all.equal(small_blocks, c(FALSE, FALSE, TRUE), check.attributes = FALSE))
})
