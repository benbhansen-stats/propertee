test_that("absorb= argument", {

  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  da <- lmitt(y ~ dose, weights = ate(), data = simdata, design = des)
  expect_length(coefficients(da), 2)

  da <- lmitt(y ~ dose, weights = ate(), data = simdata, absorb = TRUE, design = des)
  expect_true(length(coefficients(da)) > 2)

  # User passing .absorb should error
  expect_error(lmitt(y ~ .absorbed(cid1*cid2), data = simdata, design = des),
               "internal function")

})
