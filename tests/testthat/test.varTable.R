test_that("varTable", {
  data(simdata)
  des <- RCT_Design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)

  tbl <- varTable(des)

  expect_true(is(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(3, 2))

  expect_equal(tbl[,1], c("Treatment", "Unit of Assignment", "Block"))

  tbl <- varTable(des, compress = FALSE)
  expect_true(is(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(3, 3))
  expect_true(is.na(tbl[1,3]))
  expect_true(!is.na(tbl[2,3]))
  expect_true(is.na(tbl[3,3]))

  tbl <- varTable(des, reportAll = TRUE)
  expect_true(is(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(4, 2))
  expect_true(tbl[4,2] == "")
  expect_length(unique(tbl[,1]), 4)

  tbl <- varTable(des, compress = FALSE, reportAll = TRUE)
  expect_true(is(tbl, "matrix"))
  expect_type(tbl, "character")
  expect_equal(dim(tbl), c(4, 3))
  expect_true(all(is.na(c(tbl[1, 3], tbl[3, 3], tbl[4, 2:3]))))

})
