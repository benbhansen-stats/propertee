test_that("1-way tables", {
  data(simdata)

  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  tt <- design_table(des, "treatment")

  expect_s3_class(tt, "table")
  expect_length(tt, 2)
  expect_identical(as.character(unique(des@structure$z)), names(tt))

  expect_identical(tt, design_table(des, "t"))
  expect_identical(tt, design_table(des, "Treat"))

  bt <- design_table(des, "block")

  expect_s3_class(bt, "table")
  expect_length(bt, 3)
  expect_identical(as.character(unique(des@structure$bid)), names(bt))

  expect_identical(bt, design_table(des, "bl"))
  expect_identical(bt, design_table(des, "Blo"))

  ut <- design_table(des, "units of assignment")

  expect_s3_class(ut, "table")
  expect_length(dim(ut), 1)
  expect_length(ut, 10)

  expect_identical(ut, design_table(des, "unit of as"))
  expect_identical(ut, design_table(des, "un"))
  expect_identical(ut, design_table(des, "UnitS"))
  expect_identical(ut, design_table(des, "UO"))

  # Missing element
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  expect_warning(bt <- design_table(des, "block"), "blocks not found")
  expect_s3_class(ut, "table")
  expect_equal(dim(bt), 0)

  expect_error(dtable(des, "unit"), "Design specified with")


# cluster instead of uoa
  ct <- design_table(des, "clusters")

  expect_s3_class(ct, "table")
  expect_length(dim(ct), 1)
  expect_length(ct, 10)

  expect_identical(ct, design_table(des, "cl"))
  expect_identical(ct, design_table(des, "cL"))


  # sortings
  simdata$bid <- rep(1:3, times = c(10, 30, 10))
  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  bt <- design_table(des, "block")

  expect_true(!all(bt == sort(bt, decreasing = TRUE)))
  expect_true(all(names(bt) == sort(names(bt))))

  bt2 <- design_table(des, "block", sort = TRUE)
  expect_true(!all(bt == bt2))
  expect_true(all(bt2 == sort(bt2, decreasing = TRUE)))
  expect_true(!all(names(bt2) == sort(names(bt2))))

  bt3 <- design_table(des, "block", sort = TRUE, decreasing = FALSE)
  expect_true(!all(bt == bt3))
  expect_true(all(bt3 == sort(bt3, decreasing = FALSE)))
  expect_true(!all(names(bt3) == sort(names(bt3))))


  # use_var_names

  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  expect_equal(names(dimnames(design_table(des, "t")))[1], "treatment")
  expect_equal(names(dimnames(design_table(des, "b")))[1], "blocks")
  expect_equal(names(dimnames(design_table(des, "u")))[1], "units of assignment")

  expect_equal(names(dimnames(design_table(des, "t", use_var_names = TRUE)))[1],
               "z")
  expect_equal(names(dimnames(design_table(des, "b", use_var_names = TRUE)))[1],
               "bid")
  expect_equal(names(dimnames(design_table(des, "u", use_var_names = TRUE)))[1],
               "uoa1, uoa2")

  # other arguments

  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  expect_length(design_table(des, "t", useNA = "always"), 3)

  expect_equal(names(dimnames(design_table(des, "t", dnn = list("abc")))),
               "abc")

})

test_that("2 way tables", {
  data(simdata)

  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)


  tbl <- dtable(des, "trea", "bl")
  expect_s3_class(tbl, "table")
  expect_true(all.equal(dim(tbl), c(3, 2)))

  tbl2 <- dtable(des, "bl", "treat")

  expect_identical(tbl, t(tbl2))

  # Missing entry
  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)

  expect_warning(tbl <- dtable(des, "trea", "bl"), "blocks not found")
  expect_identical(tbl, dtable(des, "tr"))

  # var names
  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  tbl <- dtable(des, "trea", "bl")
  expect_equal(names(dimnames(tbl)), c("blocks", "treatment"))
  tbl <- dtable(des, "trea", "bl", use_var_names = TRUE)
  expect_equal(names(dimnames(tbl)), c("bid", "z"))

})

test_that("Passing dnn", {
  data(simdata)
  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  tbl <- dtable(des, "trea", "bl", dnn = c("a", "b"))
  expect_identical(names(attr(tbl, "dimnames")),
                   c("a", "b"))



})