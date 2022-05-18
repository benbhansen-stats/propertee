# This isn't going to test whether anything is actually correct (see individiual
# tests for that), but instead will run each different variation of `lm` call
# and ensure that no messages, warnings, or errors are produced.

test_that("running all variations of lm", {

  # Treatment exists in data
  data(simdata)

  des1 <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)

  # Weight alone
  expect_silent(a <- lm(y ~ z, data = simdata, weights = ate(des1)))

  # Weight + adopters
  expect_silent(a <- lm(y ~ adopters(), data = simdata, weights = ett(des1)))
  expect_silent(a <- lm(y ~ adopters(des1), data = simdata,
                        weights = ett(des1)))

  camod <- lm(y ~ x, data = simdata)

  # weight + cov_adj
  expect_silent(a <- lm(y ~ z, data = simdata, weights = ett(des1),
                        offset = cov_adj(camod)))
  expect_silent(a <- lm(y ~ z, data = simdata, weights = ett(),
                        offset = cov_adj(camod, design = des1)))
  expect_silent(a <- lm(y ~ z, data = simdata, weights = ett(des1),
                        offset = cov_adj(camod, design = des1)))
  expect_silent(a <- lm(y ~ z + offset(cov_adj(camod)), data = simdata,
                        weights = ett(des1)))
#  expect_silent(a <- lm(y ~ z + offset(cov_adj(camod, design = des1)),
#                        data = simdata, weights = ett()))
  expect_silent(a <- lm(y ~ z + offset(cov_adj(camod, design = des1)),
                        data = simdata, weights = ett(des1)))


  # Treatment doesn't exist in data
  data(simdata)

  des2 <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)
  simdata$z <- NULL

  # Weight + adopters
  expect_silent(a <- lm(y ~ adopters(), data = simdata, weights = ett(des2)))
  expect_silent(a <- lm(y ~ adopters(des2), data = simdata,
                        weights = ett(des2)))



})
