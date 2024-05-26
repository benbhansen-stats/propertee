test_that("apply dichotomization function", {
  zz <- data.frame(z = 1:10)

  bt <- .apply_dichotomy(zz, z > 5 ~ z == 2)
  expect_true(is.numeric(bt))
  expect_length(bt, nrow(zz))
  expect_true(all(bt %in% c(0, 1, NA)))

  expect_error(.apply_dichotomy(1, z > 5 ~ .),
               "to be a named")

  expect_error(.apply_dichotomy(zz, . ~ .),
               "least one side")

  expect_error(.apply_dichotomy(zz, z < 6 ~ z > 4),
               "overlaps")
})

test_that("apply dichotomy on factor", {
  zz <- data.frame(z = factor(sample(2:4, 20, TRUE)))

  expect_true(0 < mean(.apply_dichotomy(zz, z == "3" ~ .)))
  expect_true(mean(.apply_dichotomy(zz, z == "3" ~ .)) < 1)
  expect_identical(.apply_dichotomy(zz, z == 3 ~ .),
                   .apply_dichotomy(zz, z == "3" ~ .))

  # Ensure that we're comparing against the levels not the values - "2" is the
  # first level, but there is no level 1, so this should flag NO treated units
  expect_true(sum(.apply_dichotomy(zz, z == 1 ~ .)) == 0)

  expect_warning(.apply_dichotomy(zz, z > 3 ~ .),
                 "not meaningful")


  zz <- data.frame(z = factor(sample(c("a", "b", "c"), 20, TRUE)))

  expect_length(table(.apply_dichotomy(zz, z == "b" ~ z == "c"),
                      useNA = 'ifany'), 3)

  zz <- data.frame(z = ordered(sample(3:5, 20, TRUE)))

  expect_true(0 < mean(.apply_dichotomy(zz, z == "3" ~ .)))
  expect_true(mean(.apply_dichotomy(zz, z == "3" ~ .)) < 1)
  expect_identical(.apply_dichotomy(zz, z == 3 ~ .),
                   .apply_dichotomy(zz, z == "3" ~ .))

  expect_true(0 < mean(.apply_dichotomy(zz, z > "4" ~ .)))
  expect_true(mean(.apply_dichotomy(zz, z > "4" ~ .)) < 1)
  expect_identical(.apply_dichotomy(zz, z > 4 ~ .),
                   .apply_dichotomy(zz, z > "4" ~ .))

})
