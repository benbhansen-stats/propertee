test_that("Basic functionality", {

  # Returns the same length
  x <- rnorm(100)
  expect_identical(length(x),
                   length(areg.center(x, rep(1, length(x)))))

  # In a single block, return original
  expect_true(all.equal(x,
                        areg.center(x, rep(1, length(x)))))


})
