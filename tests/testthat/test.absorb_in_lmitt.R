test_that("absorb= argument", {

  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  da <- lmitt(y ~ 1, weights = ate(), data = simdata, design = des)
  expect_length(coefficients(da), 2)

  da <- lmitt(y ~ 1, weights = ate(), data = simdata, absorb = TRUE, design = des)
  expect_length(coefficients(da), 2)

  # subgroup effects
  da <- lmitt(y ~ as.factor(o), weights = ate(), data = simdata, design = des)
  expect_length(coefficients(da), 8)

  da <- lmitt(y ~ as.factor(o), weights = ate(), data = simdata, absorb = TRUE, design = des)
  expect_length(coefficients(da), 8)

})

test_that("multiple variables in blocks", {
  data(simdata)
  simdata$bid1 <- (simdata$bid > 1) + 1
  simdata$bid2 <- (simdata$bid != 2) + 2

  des <- rct_design(z ~ cluster(cid2, cid1) + block(bid1, bid2),
                   data = simdata)

  da <- lmitt(y ~ dose, data = simdata, absorb = TRUE, design = des)
  expect_true(length(da$coefficients) == 3)

})
