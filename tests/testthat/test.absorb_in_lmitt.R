test_that("absorb= argument", {

  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  da <- lmitt(y ~ 1, weights = ate(), data = simdata, design = des)
  expect_length(coefficients(da), 2)

  da <- lmitt(y ~ 1, weights = ate(), data = simdata, absorb = TRUE, design = des)
  expect_true(length(da$coefficients) > 2)

  # User passing .absorb should error
  expect_error(lmitt(y ~ .absorbed(cid1*cid2), data = simdata, design = des),
               "internal function")

})

test_that("absorbed effects aren't printed", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  da <- lmitt(y ~ dose, data = simdata, absorb = TRUE, design = des)

  expect_match(capture_output(show(da)), "assigned()")
  expect_no_match(capture_output(show(da)), ".absorbed")
  expect_match(capture_output(print(da$coefficients)), ".absorbed")

  sda <- summary(da)

  expect_match(capture_output(print(sda)), "assigned()")
  expect_no_match(capture_output(print(sda)), ".absorbed")
  expect_match(capture_output(print(sda$coefficients)), ".absorbed")

})

test_that("multiple variables in blocks", {
  data(simdata)
  simdata$bid1 <- simdata$bid
  simdata$bid2 <- c(rep(c(2, 3), each = 10), c(rep(2, 4), rep(3, 10)),
                    c(rep(2, 6), rep(3, 10)))

  des <- rct_design(z ~ cluster(cid2, cid1) + block(bid1, bid2),
                   data = simdata)

  da <- lmitt(y ~ dose, data = simdata, absorb = TRUE, design = des)
  expect_true(length(da$coefficients) > 2)

  expect_match(capture_output(show(da)), "assigned()")
  expect_no_match(capture_output(show(da)), ".absorbed")
  expect_match(capture_output(print(da$coefficients)), ".absorbed")

  sda <- summary(da)

  expect_match(capture_output(print(sda)), "assigned()")
  expect_no_match(capture_output(print(sda)), ".absorbed")
  expect_match(capture_output(print(sda$coefficients)), ".absorbed")

})
