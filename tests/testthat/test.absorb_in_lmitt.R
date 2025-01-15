test_that("absorb= argument", {

  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  mod <- lmitt(y ~ 1, weights = ate(), data = simdata, specification = spec)
  expect_length(coefficients(mod), 3)

  mod <- lmitt(y ~ 1, weights = ate(), data = simdata, absorb = TRUE, specification = spec)
  expect_length(coefficients(mod), 3)

  # subgroup effects
  mod <- lmitt(y ~ as.factor(o), weights = ate(), data = simdata, specification = spec)
  expect_length(coefficients(mod), 12)

  mod <- lmitt(y ~ as.factor(o), weights = ate(), data = simdata, absorb = TRUE, specification = spec)
  expect_length(coefficients(mod), 12)

})

test_that("multiple variables in blocks", {
  data(simdata)
  simdata$bid1 <- (simdata$bid > 1) + 1
  simdata$bid2 <- (simdata$bid != 2) + 2

  spec <- rct_spec(z ~ cluster(uoa2, uoa1) + block(bid1, bid2),
                   data = simdata)

  mod <- lmitt(y ~ dose, data = simdata, absorb = TRUE, specification = spec)
  expect_length(coefficients(mod), 6)

})
