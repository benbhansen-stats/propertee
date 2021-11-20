test_that("ittestimate", {

  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

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

  expect_error(ittestimate(1, simdata, "abc"),
               "must be Design")
  expect_error(ittestimate(des, 1, "abc"),
               "must be data")
  expect_error(ittestimate(des, simdata, 1),
               "must be quoted")
})

test_that("clusterIds in ittestimate", {
  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

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
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  camod <- lm(y ~ x, data = simdata)

  it <- ittestimate(des, simdata, "y", covAdjModel = camod)
  expect_true(!is.null(it$model$"offset(covAdj)"))

  camod2 <- glm(y ~ x, data = simdata)

  it <- ittestimate(des, simdata, "y", covAdjModel = camod2)
  expect_true(!is.null(it$model$"offset(covAdj)"))

  expect_error(ittestimate(des, simdata, "y", covAdjModel = 1),
               "support predict")
})

test_that("manually passing weights", {

  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force), data = simdata)

  myweights <- rep(10.1, nrow(simdata))

  it <- ittestimate(des, simdata, "y", weights = myweights*ate(des))

  expect_true(all(it$weights == myweights*ate(des, data = simdata)))

  expect_error(ittestimate(des, simdata, "y", weights = myweights),
               "must be")


  expect_error(ittestimate(des, simdata, "y", weights = 3),
               "same length")

  expect_error(ittestimate(des, simdata, "y", weights = rep("a", nrow(simdata))),
               "numeric")

})


test_that("WeightedDesign argument instead of Design", {

  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  wdes <- ate(des, simdata)

  withdes <- ittestimate(des,  simdata, "y")
  withwdes <- ittestimate(wdes, simdata, "y")

  expect_equal(withdes$coef, withwdes$coef)
  expect_identical(withdes@Design, withwdes@Design)
  expect_equal(weights(withdes), weights(withwdes))
})

test_that("treatment types", {
  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)
  it <- ittestimate(des, simdata, "y")
  expect_length(it$coef, 2)

  des <- RD_Design(o ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)
  it <- ittestimate(des, simdata, "y")
  expect_length(it$coef, 4)

  simdata$dose <- as.ordered(simdata$dose)
  des <- RD_Design(dose ~ cluster(cid2, cid1) +
                     block(bid) + forcing(force),
                   data = simdata)
  it <- ittestimate(des, simdata, "y")
  expect_length(it$coef, 5)

})
