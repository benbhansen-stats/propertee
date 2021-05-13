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

test_that("clusterIds in ittestimate", {
  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force), data = simdata)

  it1 <- ittestimate(des, simdata, "y")

  simdata2 <- simdata
  colnames(simdata2)[2] <- "cccc"
  it2 <- ittestimate(des, simdata2, "y", clusterIds = list("cid2" = "cccc"))

  expect_identical(it1$coef, it2$coef)

  expect_error(ittestimate(des, simdata2, "y",
                           clusterIds = list("abc")),
               "named list")

  expect_error(ittestimate(des, simdata2, "y",
                           clusterIds = list("cid2" = "cccc", "cid2" = "bbbb")),
               "unique")

  expect_warning(it3 <- ittestimate(des, simdata2, "y",
                           clusterIds = list("cid2" = "cccc", "abc" = "cid1")),
               "not found in Design")

  expect_warning(it4 <- ittestimate(des, simdata2, "y",
                           clusterIds = list("cid2" = "cccc", "cid1" = "abc")),
               "not found in data")
  expect_identical(it3$coef, it4$coef)

  expect_warning(expect_warning(it5 <- ittestimate(des, simdata2, "y",
                           clusterIds = list("cid2" = "cccc", "abc" = "def")),
               "not found in data"), "not found in Design")
  expect_identical(it3$coef, it5$coef)

})

test_that("covariate adjustment", {
  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force), data = simdata)

  camod <- lm(y ~ x, data = simdata)

  it <- ittestimate(des, simdata, "y", covAdjModel = camod)
  expect_true(!is.null(it$model$"offset(covAdj)"))

  camod2 <- glm(y ~ x, data = simdata)

  it <- ittestimate(des, simdata, "y", covAdjModel = camod2)
  expect_true(!is.null(it$model$"offset(covAdj)"))

  expect_error(ittestimate(des, simdata, "y", covAdjModel = 1),
               "support predict")
})
