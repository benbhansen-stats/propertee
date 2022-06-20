test_that("apply dichotomization function", {
  zz <- data.frame(z = 1:10)

  bt <- .apply_dichotomy(zz, z > 5 ~ z == 2)
  expect_true(is.numeric(bt))
  expect_length(bt, nrow(zz))
  expect_true(all(bt %in% c(0, 1, NA)))


  expect_error(.apply_dichotomy(zz, 1),
               "must be formula")

  expect_error(.apply_dichotomy(1, z > 5 ~ .),
               "to be a named")

  expect_error(.apply_dichotomy(zz, . ~ .),
               "least one side")

  expect_error(.apply_dichotomy(zz, z < 6 ~ z > 4),
               "overlaps")
})
