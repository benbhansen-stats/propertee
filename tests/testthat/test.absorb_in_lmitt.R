test_that("absorb= argument", {

  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  da <- lmitt(y ~ dose, weights = ate(des), data = simdata)
  expect_length(coefficients(da), 2)

  da <- lmitt(y ~ dose, weights = ate(des), data = simdata, absorb = TRUE)
  expect_true(length(coefficients(da)) > 2)

})
