test_that("Same result with one or two variable for blocks", {
  data(simdata)
  spec1 <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  wts1 <- ate(spec1, data = simdata)

  mod1 <- lmitt(y ~ x, data = simdata, spec = spec1, weights = wts1)
  coeff1 <- summary(mod1)$coeff


  simdata2 <- simdata
  simdata2$bid1 <- rep(0:1, times = c(34, 16))
  simdata2$bid2 <- rep(0:1, times = c(20, 30))

  # Make sure block info is same, only diagonal entries of cross table between
  # original blocks and new blocks should have values.
  expect_true(sum(diag(with(simdata2, table(bid, paste(bid1, bid2))))) == nrow(simdata2))

  spec2 <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid1, bid2), data = simdata2)
  expect_true(ncol(blocks(spec2)) == 2)

  wts2 <- ate(spec2, data = simdata2)

  expect_true(all(wts1 == wts2))

  mod2 <- lmitt(y ~ x, data = simdata2, spec = spec2, weights = wts2)
  coeff2 <- summary(mod2)$coeff

  all(coeff1 == coeff2)

})
