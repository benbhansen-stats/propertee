test_that(".get_dichotomy with no lmitt.formula() and no arugment", {
  expect_true(is.null(.get_dichotomy()))
})

test_that(".get_dichotomy with no lmitt.formula()", {
  expect_true(all.equal(.get_dichotomy(z ~ 1), z ~ 1))
})