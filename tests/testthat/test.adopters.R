test_that("basic adopters", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  lm(y ~ z, weights = ate(des), data = simdata)
  lm(y ~ adopters(), weights = ate(des), data = simdata)
})
