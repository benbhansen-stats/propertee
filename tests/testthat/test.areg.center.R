test_that("Basic functionality, grand_mean_center=TRUE", {

  mm <- matrix(rep(seq(2), 6), ncol = 1)

  # Returns the same dimension
  expect_identical(dim(mm),
                   dim(center_mm <- areg.center(mm, rep(1, nrow(mm)),
                                                grand_mean_center = TRUE)))

  # In a single block, return original
  expect_true(all.equal(mm, center_mm))

  # correct values for one column with multiple groups but no weights
  grps <- rep(seq(3), each = 4)
  expect_equal(areg.center(mm, grps, grand_mean_center = TRUE)[,1],
               rep(c(-0.5 + 1.5, 0.5 + 1.5), 6))

  # correct values for one column with multiple groups and weights
  wts <- rep(seq(2), 6)
  expect_equal(areg.center(mm, grps, wts, grand_mean_center = TRUE)[,1],
               rep(c(-2/3 + 5/3, 1/3 + 5/3), 6))
  
  # correct values for multiple columns with multiple groups and weights
  mm <- cbind(mm, mm)
  expect_identical(dim(mm), dim(center_mm <- areg.center(mm, grps, wts, grand_mean_center = TRUE)))
  expect_true(
    all.equal(center_mm,
              cbind(rep(c(-2/3 + 5/3, 1/3 + 5/3), 6),
                    rep(c(-2/3 + 5/3, 1/3 + 5/3), 6)),
              check.attributes = FALSE)
  )
  
  # correct values when a group has 0 weights for all units
  wts <- c(rep(0, 4), rep(seq(2), 4))
  expect_true(
    all.equal(areg.center(mm, grps, wts, grand_mean_center = TRUE),
              cbind(c(rep(c(1 + 5/3, 2 + 5/3), 2), rep(c(-2/3 + 5/3, 1/3 + 5/3), 4)),
                    c(rep(c(1 + 5/3, 2 + 5/3), 2), rep(c(-2/3 + 5/3, 1/3 + 5/3), 4))),
              check.attributes = FALSE)
  )
})

test_that("Basic functionality, grand_mean_center=FALSE", {

  mm <- matrix(rep(seq(2), 6), ncol = 1)
  
  # Returns the same dimension
  expect_identical(dim(mm),
                   dim(center_mm <- areg.center(mm, rep(1, nrow(mm)))))
  
  # In a single block, return original with grand mean subtracted out
  expect_true(all.equal(mm - mean(mm), center_mm))
  
  # correct values for one column with multiple groups but no weights
  grps <- rep(seq(3), each = 4)
  expect_equal(areg.center(mm, grps)[,1],
               rep(c(-0.5, 0.5), 6))
  
  # correct values for one column with multiple groups and weights
  wts <- rep(seq(2), 6)
  expect_equal(areg.center(mm, grps, wts)[,1],
               rep(c(-2/3, 1/3), 6))
  
  # correct values for multiple columns with multiple groups and weights
  mm <- cbind(mm, mm)
  expect_identical(dim(mm),
                   dim(center_mm <- areg.center(mm, grps, wts)))
  expect_true(
    all.equal(center_mm,
              cbind(rep(c(-2/3, 1/3), 6), rep(c(-2/3, 1/3), 6)),
              check.attributes = FALSE)
  )
  
  # correct values when a group has 0 weights for all units
  wts <- c(rep(0, 4), rep(seq(2), 4))
  expect_true(
    all.equal(areg.center(mm, grps, wts),
              cbind(c(rep(c(1, 2), 2), rep(c(-2/3, 1/3), 4)),
                    c(rep(c(1, 2), 2), rep(c(-2/3, 1/3), 4))),
              check.attributes = FALSE)
  )
})
