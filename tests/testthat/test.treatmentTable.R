test_that("treatment_table", {
  data(simdata)

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)

  tt <- treatment_table(des)

  expect_s3_class(tt, "table")
  expect_length(tt, 2)
  expect_identical(as.character(unique(des@structure$z)), names(tt))

  des <- rct_design(o ~ uoa(cid1, cid2), data = simdata)

  tt <- treatment_table(des)

  expect_length(tt, 4)
  expect_identical(sort(as.character(unique(des@structure$o))), sort(names(tt)))
  expect_identical(sum(tt), nrow(des@structure))


  # additional arguments

  tt <- treatment_table(des, useNA = "always")
  expect_length(tt, 5)
  expect_identical(sum(tt), nrow(des@structure))


})
