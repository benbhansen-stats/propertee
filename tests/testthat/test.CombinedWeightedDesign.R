test_that("Combining weighted designs", {
  data(simdata)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)

  w1 <- ate(des, data = simdata[1:30,])
  w2 <- ate(des, data = simdata[31:40,])
  w3 <- ate(des, data = simdata[41:50,])

  c_w <- c(w1, w2, w3)
  expect_true(inherits(c_w, "WeightedDesign"))
  expect_length(c_w, 50)
  expect_identical(c_w, ate(des, data = simdata))

  w1e <- ett(des, data = simdata[1:30,])
  w2e <- ett(des, data = simdata[31:40,])
  w3e <- ett(des, data = simdata[41:50,])

  c_we <- c(w1e, w2e, w3e)
  expect_true(inherits(c_we, "WeightedDesign"))
  expect_length(c_we, 50)
  expect_identical(c_we, ett(des, data = simdata))

  expect_error(c(w1, 1:5), "with other")
  expect_error(c(w1, w1e), "same target")

  des2 <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  alt_w1 <- ate(des2, data = simdata)

  expect_error(c(w1, alt_w1), "which differ on elements")

  # if the first argument is compatible with WeightedDesign but isn't one (e.g.
  # numeric vector), c() will return a numeric vector
  #expect_true(inherits(c(1:5, w1), "WeightedDesign"))
})

test_that("Combining weighted designs with different dichotomys ", {
  des <- rct_design(dose ~ uoa(uoa1, uoa2), data = simdata)

  w1 <- ate(des, data = simdata[1:10, ], dichotomy = dose >= 300 ~ .)
  w2 <- ate(des, data = simdata[11:30, ], dichotomy = dose >= 200 ~ .)
  w3 <- ate(des, data = simdata[31:50, ], dichotomy = dose >= 100 ~ .)

  c_w <- c(w1, w2, w3)
  expect_true(inherits(c_w, "WeightedDesign"))
  expect_length(c_w, 50)

  expect_error(c(w1, w2, w3, force_dichotomy_equal = TRUE),
               "must be identical")
})
