test_that("ittestimate", {

  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force), data = simdata)

  it <- ittestimate(des, simdata, "y")

  expect_s4_class(it, "DirectAdjusted")
  expect_true(is(it, "lm"))

  it_ate <- ittestimate(des, simdata, "y", "ate")
  it_ett <- ittestimate(des, simdata, "y", "ett")

  expect_identical(it$coef, it_ate$coef)
  expect_equal(it@target, "ate")
  expect_equal(it_ate@target, "ate")
  expect_equal(it_ett@target, "ett")


  expect_error(ittestimate(des, simdata, "abc"),
               "must be a column")
  expect_error(ittestimate(des, simdata, "y", "abc"),
               "invalid target")


})
