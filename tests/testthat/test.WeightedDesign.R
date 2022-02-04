test_that("ate and ett with data argument", {

  # n_clusters = n
  data(mtcars)
  mtcars <- mtcars[-c(5, 11),]
  des <- RD_Design(am ~ cluster(qsec) + forcing(vs), data = mtcars)

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

test_that("unitOfAssignmentIds", {

  data(simdata)
  des <- RD_Design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force), data = simdata)

  w1 <- ate(des, data = simdata)

  simdata2 <- simdata
  colnames(simdata2)[2] <- "cccc"
  w2 <- ate(des, data = simdata2, unitOfAssignmentIds = list("cid2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)

  w1 <- ett(des, data = simdata)
  w2 <- ett(des, data = simdata2, unitOfAssignmentIds = list("cid2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)

  mod1 <- lm(y ~ x, data = simdata, weights = ate(des))
  mod2 <- lm(y ~ x, data = simdata2, weights = ate(des,unitOfAssignmentIds = list("cid2" = "cccc")))

  expect_identical(mod1$coef, mod2$coef)

  mod1 <- lm(y ~ x, data = simdata, weights = ett(des))
  mod2 <- lm(y ~ x, data = simdata2, weights = ett(des,unitOfAssignmentIds = list("cid2" = "cccc")))

  expect_identical(mod1$coef, mod2$coef)

})

test_that("Ops", {
  data(simdata)
  des <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_s4_class(wdes * 2, "WeightedDesign")
  expect_s4_class(2 * wdes, "WeightedDesign")
  expect_identical(wdes * 3, 3 * wdes)
  expect_equal((wdes * 3)@.Data, wdes@.Data * 3)
  expect_error(wdes * -1, "non-negative")

  expect_s4_class(wdes / 2, "WeightedDesign")
  expect_s4_class(2 / wdes, "WeightedDesign")
  expect_equal((wdes / 3)@.Data, wdes@.Data / 3)

  expect_error(1 + wdes)
  expect_error(wdes + 1)

  expect_error(1 - wdes)
  expect_error(wdes - 1)

  expect_s4_class(wdes ^ 3, "WeightedDesign")
  expect_equal((wdes ^ 3)@.Data, wdes@.Data ^ 3)


})

test_that("show WeightedDesign", {

  data(simdata)
  des <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_silent(invisible(capture.output(expect_identical(wdes, show(wdes)))))
})

test_that("weight function", {
  data(simdata)
  des <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_s4_class(wdes, "WeightedDesign")
  expect_false(is(weights(wdes), "WeightedDesign"))

  })

test_that("Inconsistent treatment levels", {

  data(simdata)

  # Extra level in data into ett, not in design
  des <- RCT_Design(o ~ cluster(cid1, cid2), data = simdata[simdata$o != 4, ])
  expect_error(ett(des, data = simdata), "not found in Design")
  expect_error(ate(des, data = simdata), "not found in Design")

  expect_error(lm(y ~ x, data = simdata, weights = ett(des)) ,
                 "not found in Design")

  # Extra level in design,not in data to ett
  des <- RCT_Design(o ~ cluster(cid1, cid2), data = simdata)
  expect_warning(ett(des, data = simdata[simdata$o != 4, ]),
                 "not found in data")
  expect_warning(ate(des, data = simdata[simdata$o != 4, ]),
                 "not found in data")

  expect_warning(lm(y ~ x, data = simdata[simdata$o != 4, ], weights = ett(des)) ,
                 "not found in data")

  # Treatment var not found
  simdata$o <- NULL
  expect_error(ett(des, data = simdata), "not found in data")

})

test_that("ett treatment weights are correct length", {
  
  testdata <- data.frame(cid = 1:10, z = c(rep(1,4), rep(0,6)))
  
  des <- RCT_Design(z ~ cluster(cid), data = testdata)
  
  wts <- ett(des, data = testdata)
  
  expect_equal(length(wts), nrow(testdata))
  
})

test_that("ett treatment weights return numeric", {
  
  testdata <- data.frame(cid = 1:10, z = c(rep(1,4), rep(0,6)))
  
  des <- RCT_Design(z ~ cluster(cid), data = testdata)
  
  wts <- ett(des, data = testdata)
  
  expect_s4_class(wts, "numeric")
  
})

test_that("ett treatment weights = 1", {
  
  testdata <- data.frame(cid = 1:10, z = c(rep(1,4), rep(0,6)))
  
  des <- RCT_Design(z ~ cluster(cid), data = testdata)
  
  wts <- ett(des, data = testdata)
  
  expect_equal(wts@.Data[1:4], rep(1,4))
  
})

test_that("ett weights for block with P(Z = 1) = 0.5: 1", {
  
  testdata <- data.frame(cid = 1:10, z = c(rep(1,5), rep(0,5)))
  
  des <- RCT_Design(z ~ cluster(cid), data = testdata)
  
  wts <- ett(des, data = testdata)
  
  expect_equal(wts@.Data, rep(1,10))
  
})

test_that("ett weights for block with P(Z = 1) = 1/3:  1 or 0.5", {
  
  testdata <- data.frame(cid = 1:30, z = c(rep(1,10), rep(0,20)))
  
  des <- RCT_Design(z ~ cluster(cid), data = testdata)
  
  wts <- ett(des, data = testdata)
  
  expect_equal(wts@.Data, c(rep(1,10), rep(0.5, 20)))
  
})

test_that("combine two previous blocks, obtain proper weights", {
  
  testdata <- data.frame(bid = c(rep(1,10), rep(2, 30)),
                                 cid = 1:40, 
                                 z = c(rep(1,5), rep(0,5), 
                                       rep(1,10), rep(0, 20)))
  
  des <- RCT_Design(z ~ cluster(cid) + block(bid), data = testdata)
  
  wts <- ett(des, data = testdata)
  
  expect_equal(wts@.Data, c(rep(1,20), rep(0.5, 20)))
  
})

