test_that("DirectAdjusted object", {

  data(simdata)
  des <- Obs_Design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- DirectAdjusted(lm(y ~ x, data = simdata, weights = ate(des)), Design = des, target = "ett")

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(is(dalm, "lm"))

  expect_identical(dalm$model$"(weights)"@Design, des)
  expect_identical(dalm$model$"(weights)"@Design, dalm@Design)
})
