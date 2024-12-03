test_that("Finding cov_adj in model", {

  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  ca <- lm(y ~ x, data = simdata)

  cap <- cov_adj(ca, newdata = simdata, specification = spec)

  mod1 <- lm(y ~ z, data = simdata, weights = ate(spec), offset = cov_adj(ca))
  expect_identical(.get_cov_adj(mod1), cap)

  mod2 <- lm(y ~ z  + offset(cov_adj(ca)), data = simdata, weights = ate(spec))
  expect_identical(.get_cov_adj(mod2), cap)


  # Weird variable named offset
  names(simdata)[8] <- "offset(abc)"
  mod3 <- lm(y ~ z + `offset(abc)`, data = simdata)
  expect_null(.get_cov_adj(mod3))

  # An offset not from us
  mod4 <- lm(y ~ z + offset(force), data = simdata)
  expect_null(.get_cov_adj(mod4))
  mod5 <- lm(y ~ z, data = simdata, offset = force)
  expect_null(.get_cov_adj(mod5))

  # proper offset, but weird name cov_adj variable
  names(simdata)[8] <- "cov_adj"
  mod6 <- lm(y ~ z + offset(cov_adj), data = simdata)
  expect_null(.get_cov_adj(mod6))

  # Double offset
  names(simdata)[8] <- "x"
  mod7 <- lm(y ~ z + offset(cov_adj(ca)) + offset(x), data = simdata,
             weights = ate(spec))
  expect_identical(.get_cov_adj(mod7), cap)

  # Double offset, both with the identical cov_adj
  mod8 <- lm(y ~ z + offset(cov_adj(ca)) + offset(cap), data = simdata,
             weights = ate(spec))
  expect_identical(.get_cov_adj(mod8), cap)

  # Double offset, but with different cov_adj
  ca2 <- lm(y ~ force, data = simdata)
  mod9 <- lm(y ~ z + offset(cov_adj(ca)) + offset(cov_adj(ca2)), data = simdata,
             weights = ate(spec))
  expect_error(.get_cov_adj(mod9), "Multiple cov_adj")
})
