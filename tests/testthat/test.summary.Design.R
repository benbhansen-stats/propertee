test_that("summary.Design", {
  data(simdata)

  desrct <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  desrd  <-  rd_design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                       data = simdata)
  desobs <- obs_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  expect_s3_class(summary(desrct), "summary.Design")
  expect_s3_class(summary(desrd) , "summary.Design")
  expect_s3_class(summary(desobs), "summary.Design")

  expect_s4_class(summary(desrct)$Design, "Design")
  expect_s4_class(summary(desrd)$Design , "Design")
  expect_s4_class(summary(desobs)$Design, "Design")

  expect_output(print(summary(desrct)), "Randomized")
  expect_output(print(summary(desrd)) , "Discontinuity")
  expect_output(print(summary(desobs)), "Observational")

  desrct <- rct_design(o ~ cluster(cid1, cid2) + block(bid), data = simdata)
  expect_output(print(summary(desrct)), "...")
  expect_output(print(summary(desrct)), "excluded")
  expect_identical(desrct, summary(desrct)$Design)
  expect_s3_class(summary(desrct)$treatment_table, "table")

  expect_output(print(summary(desrct)),
                "group excluded")
  expect_output(print(summary(desrct), max_unit_print = 2),
                "groups excluded")

})
