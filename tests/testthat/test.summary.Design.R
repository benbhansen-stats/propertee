test_that("summary.Design", {
  data(simdata)

  desrct <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  desrd  <-  RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force), data = simdata)
  desobs <- Obs_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  expect_output(summary(desrct), "Randomized")
  expect_output(summary(desrd), "Discontinuity")
  expect_output(summary(desobs), "Observational")

})
