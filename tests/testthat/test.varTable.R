test_that("var_table", {
  data(simdata)
  des <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)

  tbl <- var_table(des)

  expect_true(is(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(3, 2))

  expect_equal(tbl[, 1], c("Treatment", "Unit of Assignment", "Block"))

  tbl <- var_table(des, compress = FALSE)
  expect_true(is(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(3, 3))
  expect_true(is.na(tbl[1, 3]))
  expect_true(!is.na(tbl[2, 3]))
  expect_true(is.na(tbl[3, 3]))

  tbl <- var_table(des, report_all = TRUE)
  expect_true(is(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(4, 2))
  expect_true(tbl[4, 2] == "")
  expect_length(unique(tbl[, 1]), 4)

  tbl <- var_table(des, compress = FALSE, report_all = TRUE)
  expect_true(is(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(4, 3))
  expect_true(all(is.na(c(tbl[1, 3], tbl[3, 3], tbl[4, 2:3]))))

})
