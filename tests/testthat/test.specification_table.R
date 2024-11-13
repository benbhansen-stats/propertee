test_that("1-way tables", {
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  tt <- specification_table(spec, "treatment")

  expect_s3_class(tt, "table")
  expect_length(tt, 2)
  expect_identical(as.character(unique(spec@structure$z)), names(tt))

  expect_identical(tt, specification_table(spec, "t"))
  expect_identical(tt, specification_table(spec, "Treat"))

  bt <- specification_table(spec, "block")

  expect_s3_class(bt, "table")
  expect_length(bt, 3)
  expect_identical(as.character(unique(spec@structure$bid)), names(bt))

  expect_identical(bt, specification_table(spec, "bl"))
  expect_identical(bt, specification_table(spec, "Blo"))

  ut <- specification_table(spec, "units of assignment")

  expect_s3_class(ut, "table")
  expect_length(dim(ut), 1)
  expect_length(ut, 10)

  expect_identical(ut, specification_table(spec, "unit of as"))
  expect_identical(ut, specification_table(spec, "un"))
  expect_identical(ut, specification_table(spec, "UnitS"))
  expect_identical(ut, specification_table(spec, "UO"))

  # Missing element
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  expect_warning(bt <- specification_table(spec, "block"), "blocks not found")
  expect_s3_class(ut, "table")
  expect_equal(dim(bt), 0)

  expect_error(stable(spec, "unit"), "StudySpecification specified with")


# cluster instead of uoa
  ct <- specification_table(spec, "clusters")

  expect_s3_class(ct, "table")
  expect_length(dim(ct), 1)
  expect_length(ct, 10)

  expect_identical(ct, specification_table(spec, "cl"))
  expect_identical(ct, specification_table(spec, "cL"))


  # sortings
  simdata$bid <- rep(1:3, times = c(10, 30, 10))
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  bt <- specification_table(spec, "block")

  expect_true(!all(bt == sort(bt, decreasing = TRUE)))
  expect_true(all(names(bt) == sort(names(bt))))

  bt2 <- specification_table(spec, "block", sort = TRUE)
  expect_true(!all(bt == bt2))
  expect_true(all(bt2 == sort(bt2, decreasing = TRUE)))
  expect_true(!all(names(bt2) == sort(names(bt2))))

  bt3 <- specification_table(spec, "block", sort = TRUE, decreasing = FALSE)
  expect_true(!all(bt == bt3))
  expect_true(all(bt3 == sort(bt3, decreasing = FALSE)))
  expect_true(!all(names(bt3) == sort(names(bt3))))


  # use_var_names

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  expect_equal(names(dimnames(specification_table(spec, "t")))[1], "treatment")
  expect_equal(names(dimnames(specification_table(spec, "b")))[1], "blocks")
  expect_equal(names(dimnames(specification_table(spec, "u")))[1], "units of assignment")

  expect_equal(names(dimnames(specification_table(spec, "t", use_var_names = TRUE)))[1],
               "z")
  expect_equal(names(dimnames(specification_table(spec, "b", use_var_names = TRUE)))[1],
               "bid")
  expect_equal(names(dimnames(specification_table(spec, "u", use_var_names = TRUE)))[1],
               "uoa1, uoa2")

  # other arguments

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  expect_length(specification_table(spec, "t", useNA = "always"), 3)

  expect_equal(names(dimnames(specification_table(spec, "t", dnn = list("abc")))),
               "abc")

})

test_that("2 way tables", {
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)


  tbl <- stable(spec, "trea", "bl")
  expect_s3_class(tbl, "table")
  expect_true(all.equal(dim(tbl), c(3, 2)))

  tbl2 <- stable(spec, "bl", "treat")

  expect_identical(tbl, t(tbl2))

  # Missing entry
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)

  expect_warning(tbl <- stable(spec, "trea", "bl"), "blocks not found")
  expect_identical(tbl, stable(spec, "tr"))

  # var names
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  tbl <- stable(spec, "trea", "bl")
  expect_equal(names(dimnames(tbl)), c("blocks", "treatment"))
  tbl <- stable(spec, "trea", "bl", use_var_names = TRUE)
  expect_equal(names(dimnames(tbl)), c("bid", "z"))

})

test_that("Passing dnn", {
  data(simdata)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  tbl <- stable(spec, "trea", "bl", dnn = c("a", "b"))
  expect_identical(names(attr(tbl, "dimnames")),
                   c("a", "b"))



})
