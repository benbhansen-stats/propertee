test_that("basic adopters", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lm(y ~ z, weights = ate(des), data = simdata)
  mod2 <- lm(y ~ adopters(), weights = ate(des), data = simdata)
  expect_true(all(coef(mod1) == coef(mod2)))

  mod3 <- lm(y ~ z, weights = ate(), offset = cov_adj(mod1, design = des),
             data = simdata)
  mod4 <- lm(y ~ z, weights = ate(des), offset = cov_adj(mod1),
             data = simdata)
  mod5 <- lm(y ~ adopters(), weights = ate(),
             offset = cov_adj(mod1, design = des),
             data = simdata)
  mod6 <- lm(y ~ adopters(), weights = ate(des), offset = cov_adj(mod1),
             data = simdata)
})
