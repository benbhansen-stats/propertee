test_that("summary.Design", {
  data(simdata)

  desrct <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  desrd  <-  RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force), data = simdata)
  desobs <- Obs_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  expect_s3_class(summary(desrct), "summary.Design")
  expect_s3_class(summary(desrd) , "summary.Design")
  expect_s3_class(summary(desobs), "summary.Design")

  expect_s4_class(summary(desrct)$Design, "Design")
  expect_s4_class(summary(desrd)$Design , "Design")
  expect_s4_class(summary(desobs)$Design, "Design")

  expect_output(print(summary(desrct)), "Randomized")
  expect_output(print(summary(desrd)) , "Discontinuity")
  expect_output(print(summary(desobs)), "Observational")

})
