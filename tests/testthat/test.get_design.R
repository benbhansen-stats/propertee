test_that(".get_design", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata)

  mod1 <- lm(y ~ z, data = simdata, weights = ate(des),
             offset = cov_adj(mod))
  mod2 <- lm(y ~ z, data = simdata, weights = ate(),
             offset = cov_adj(mod, design = des))
  mod3 <- lm(y ~ z, data = simdata, weights = ate(des),
             offset = cov_adj(mod, design = des))
  mod4 <- lm(y ~ adopters(), data = simdata, weights = ate(des),
             offset = cov_adj(mod))
  mod5 <- lm(y ~ adopters(), data = simdata, weights = ate(),
             offset = cov_adj(mod, design = des))
  mod6 <- lm(y ~ adopters(), data = simdata, weights = ate(des),
             offset = cov_adj(mod, design = des))
  expect_true(all(mod1$coef == mod2$coef))
  expect_true(all(mod1$coef == mod3$coef))
  expect_true(all(mod1$coef == mod4$coef))
  expect_true(all(mod1$coef == mod5$coef))
  expect_true(all(mod1$coef == mod6$coef))

  des2 <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  expect_error(lm(y ~ adopters(), data = simdata, weights = ate(des),
                  offset = cov_adj(mod, design = des2)),
               "but differ")

  expect_error(lm( y ~ adopters(), data = simdata),
               "Unable to locate")

  # #37 offset in formula isntead of argument
  mod7 <- lm(y ~ z + offset(cov_adj(mod)), data = simdata, weights = ate(des))
  expect_true(all(mod1$coef == mod7$coef))
})

test_that(".get_design returns NULL with NULL_on_error = TRUE", {
  on.exit(cov_adj <- flexida::cov_adj)
  
  cov_adj <- function(model, newdata = NULL) {
    design <- .get_design(TRUE)
    
    return(stats::predict(model, type = "response", newdata = newdata))
  }
  
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  mod <- lm(y ~ x, data = simdata)
  expect_message(cov_adj(mod, newdata = simdata),
                 "Design in call stack.")
  expect_message(lm(y ~ z, data = simdata, offset = cov_adj(mod)),
                 "Design in call stack.")
})
