test_that("ate and ett with data argument", {

  # n_clusters = n
  data(mtcars)
  mtcars <- mtcars[-c(5, 11),]
  des <- RD_Design(am ~ cluster(qsec) + forcing(am), data = mtcars)

  wdes <- ate(des, data = mtcars)

  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ate")

  expect_equal(nrow(mtcars), length(wdes))
  expect_true(all(wdes == wdes@.Data))


  # n_clusters < n
  data(simdata)
  des <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ett")

  expect_equal(nrow(simdata), length(wdes))
  expect_true(all(wdes == wdes@.Data))

})

test_that("ate and ett in lm call", {
  data(simdata)
  des <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata, weights = ate(des))

  expect_equal(mod$weights, ate(des, simdata)@.Data)

  mod <- lm(y ~ x, data = simdata, weights = ett(des))

  expect_equal(mod$weights, ett(des, simdata)@.Data)

})
