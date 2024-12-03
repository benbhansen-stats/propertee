test_that("var_table", {
  data(simdata)
  # RCT and Obs
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  tbl <- var_table(spec)

  expect_true(inherits(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(3, 2))

  expect_equal(tbl[, 1], c("Treatment", "Unit of Assignment", "Block"))

  # Only find multiple variables in cluster
  expect_equal(grepl(",", tbl[, 2]), c(FALSE, TRUE, FALSE))

  tbl2 <- var_table(spec, compress = FALSE)
  expect_true(inherits(tbl2, "matrix"))
  expect_type(tbl2, "character")
  expect_equal(colnames(tbl2)[-1], c("Variable 1", "Variable 2"))
  expect_equal(dim(tbl2), c(3, 3))
  expect_true(is.na(tbl2[1, 3]))
  expect_true(!is.na(tbl2[2, 3]))
  expect_true(is.na(tbl2[3, 3]))

  tbl3 <- var_table(spec, report_all = TRUE)
  expect_identical(tbl, tbl3)

  tbl4 <- var_table(spec, compress = FALSE, report_all = TRUE)
  expect_identical(tbl2, tbl4)

  # RD
  spec <- rd_spec(z ~ uoa(uoa1, uoa2) + forcing(force), data = simdata)

  tbl <- var_table(spec)

  expect_true(inherits(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(3, 2))

  expect_equal(tbl[, 1], c("Treatment", "Unit of Assignment", "Forcing"))

  # Only find multiple variables in cluster
  expect_equal(grepl(",", tbl[, 2]), c(FALSE, TRUE, FALSE))

  tbl2 <- var_table(spec, compress = FALSE)
  expect_true(inherits(tbl2, "matrix"))
  expect_type(tbl2, "character")
  expect_equal(dim(tbl2), c(3, 3))
  expect_true(is.na(tbl2[1, 3]))
  expect_true(!is.na(tbl2[2, 3]))
  expect_true(is.na(tbl2[3, 3]))

  tbl3 <- var_table(spec, report_all = TRUE)
  expect_true(inherits(tbl3, "matrix"))
  expect_type(tbl3, "character")
  expect_equal(dim(tbl3), c(4, 2))
  expect_true(!is.na(tbl3[1, 2]))
  expect_true(!is.na(tbl3[2, 2]))
  expect_true(is.na(tbl3[3, 2]))
  expect_true(!is.na(tbl3[4, 2]))

  tbl4 <- var_table(spec, compress = FALSE, report_all = TRUE)
  expect_true(inherits(tbl4, "matrix"))
  expect_type(tbl4, "character")
  expect_equal(dim(tbl4), c(4, 3))
  expect_true(!is.na(tbl4[1, 2]))
  expect_true(!is.na(tbl4[2, 2]))
  expect_true(is.na(tbl4[3, 2]))
  expect_true(!is.na(tbl4[4, 2]))
  expect_true(is.na(tbl4[1, 3]))
  expect_true(!is.na(tbl4[2, 3]))
  expect_true(is.na(tbl4[3, 3]))
  expect_true(is.na(tbl4[4, 3]))

  # No double-variables
  simdata$cid <- paste(simdata$uoa1, simdata$uoa2)
  spec <- rct_spec(z ~ uoa(cid), data = simdata)

  expect_identical(var_table(spec),
                   var_table(spec, compress = FALSE))

})
