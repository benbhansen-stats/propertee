test_that("treatmentTable", {
  data(simdata)

  des <- RCT_Design(z ~ uoa(cid1, cid2), data = simdata)

  tt <- treatmentTable(des)

  expect_s3_class(tt, "table")
  expect_length(tt, 2)
  expect_identical(as.character(unique(des@structure$z)), names(tt))
  expect_identical(table(des@structure$z), tt)

  des <- RCT_Design(o ~ uoa(cid1, cid2), data = simdata)

  tt <- treatmentTable(des)

  expect_length(tt, 4)
  expect_identical(sort(as.character(unique(des@structure$o))), sort(names(tt)))
  expect_identical(sum(tt), nrow(des@structure))


  # additional arguments

  tt <- treatmentTable(des, useNA = "always")
  expect_length(tt, 5)
  expect_identical(sum(tt), nrow(des@structure))


})
