test_that("Design formula checking", {
  expect_true(checkDesignFormula(y ~ cluster(x)))
  expect_true(checkDesignFormula(y ~ cluster(x,z,q,r)))

  expect_error(checkDesignFormula(~ cluster(x)),
               "treatment")

  expect_error(checkDesignFormula(y ~ x),
               "clustering")

  expect_error(checkDesignFormula(y ~ cluster(x) + cluster(z)),
               "only one cluster")

  expect_true(checkDesignFormula(y ~ cluster(x) + block(z)))
  expect_true(checkDesignFormula(y ~ cluster(x) + block(z, a, b, c)))

  expect_error(checkDesignFormula(y ~ cluster(x) + block(z) + block(q)),
               "only one block")

  expect_error(checkDesignFormula(y ~ cluster(x) + forcing(z)),
               "only allowed")

  expect_true(checkDesignFormula(y ~ cluster(x) + forcing(z),
                                 allowForcing = TRUE))

  expect_error(checkDesignFormula(y ~ cluster(x) + forcing(z) + forcing(q),
                                  allowForcing = TRUE),
               "only one forcing")

})
