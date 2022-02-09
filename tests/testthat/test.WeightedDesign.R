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

  simshort <- simdata[simdata$o != 4,]
  simshort$o <- droplevels(simshort$o)

  # Extra level in data into ett, not in design
  des <- RCT_Design(o ~ cluster(cid1, cid2), data = simshort)
  expect_error(ett(des, data = simdata), "not found in Design")
  expect_error(ate(des, data = simdata), "not found in Design")

  expect_error(lm(y ~ x, data = simdata, weights = ett(des)),
                 "not found in Design")

  # Extra level in design,not in data to ett
  des <- RCT_Design(o ~ cluster(cid1, cid2), data = simdata)
  # issue18 - extra expect_warning
  expect_warning(expect_warning(ett(des, data = simshort),
                                 "not found in data"))
  expect_warning(expect_warning(ate(des, data = simshort),
                                "not found in data"))

  expect_warning(expect_warning(lm(y ~ x, data = simshort,
                                   weights = ett(des)) ,
                 "not found in data"))
  # issue18 end

  # Treatment var not found
  simdata$o <- NULL
  expect_error(ett(des, data = simdata), "not found in data")

  levels(simshort$o) <- 1:6
  # issue18 extra expect_warning
  expect_warning(expect_warning(expect_warning(ate(des, data = simshort),
                                               "Empty levels"),
                                "not found in data"))
  # issue18 end

})

test_that("ett/ate treatment weights are correct length", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  des <- RCT_Design(z ~ unitOfAssignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(length(wts), nrow(testdata))

  wts <- ate(des, data = testdata)

  expect_equal(length(wts), nrow(testdata))

})

test_that("ett/ate treatment weights return numeric", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  des <- RCT_Design(z ~ unitOfAssignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_s4_class(wts, "numeric")

  wts <- ate(des, data = testdata)

  expect_s4_class(wts, "numeric")

})

test_that("ett treatment weights = 1", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  des <- RCT_Design(z ~ unitOfAssignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(wts@.Data[1:4], rep(1, 4))

})

test_that("ett weights for block with P(Z = 1) = 0.5: 1", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)))

  des <- RCT_Design(z ~ unitOfAssignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(wts@.Data, rep(1, 10))

})

test_that("ett weights for block with P(Z = 1) = 1/3:  1 or 0.5", {

  testdata <- data.frame(cid = 1:30, z = c(rep(1, 10), rep(0, 20)))

  des <- RCT_Design(z ~ unitOfAssignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(wts@.Data, c(rep(1, 10), rep(0.5, 20)))

})

test_that("ate weights for block with P(Z = 1) = 0.5: 2", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)))

  des <- RCT_Design(z ~ unitOfAssignment(cid), data = testdata)

  wts <- ate(des, data = testdata)

  expect_equal(wts@.Data, rep(2, 10))

})

test_that("ate weights for block with P(Z = 1) = 1/3:  3 or 1.5", {

  testdata <- data.frame(cid = 1:30, z = c(rep(1, 10), rep(0, 20)))

  des <- RCT_Design(z ~ unitOfAssignment(cid), data = testdata)

  wts <- ate(des, data = testdata)

  expect_equal(wts@.Data, c(rep(3, 10), rep(1.5, 20)))

})

test_that("combine two previous blocks, obtain proper weightsfor ett and ate", {

  testdata <- data.frame(bid = c(rep(1, 10), rep(2, 30)),
                         cid = 1:40,
                         z = c(rep(1, 5), rep(0, 5),
                               rep(1, 10), rep(0, 20)))

  des <- RCT_Design(z ~ unitOfAssignment(cid) + block(bid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(wts@.Data, c(rep(1, 20), rep(0.5, 20)))

  wts <- ate(des, data = testdata)

  expect_equal(wts@.Data, c(rep(2, 10), rep(3, 10), rep(1.5, 20)))

})

test_that("formula as object is able to be found in data search", {
  data(simdata)

  des1 <- RCT_Design(z ~ uoa(cid1, cid2), data = simdata)
  mod1 <- lm(y ~ x, data = simdata, weights = ate(des1))

  f <- z ~ uoa(cid1, cid2)
  des2 <- RCT_Design(f, data = simdata)
  mod2 <- lm(y ~ x, data = simdata, weights = ate(des2))

  expect_true(all(mod1$coef == mod2$coef))


})

test_that("varying treatment types", {
  data(simdata)
  s2 <- simdata

  # ways of representing binary treatments
  # numeric
  des1 <- RCT_Design(z ~ uoa(cid1, cid2), data = s2)
  mod1 <- lm(y ~ z, data = s2, weights = ate(des1))
  # factor
  s2$z <- as.factor(s2$z)
  des2 <- RCT_Design(z ~ uoa(cid1, cid2), data = s2)
  mod2 <- lm(y ~ z, data = s2, weights = ate(des2))
  expect_identical(mod1$weights, mod2$weights)
  # ordinal
  s2$z <- as.ordered(s2$z)
  des3 <- RCT_Design(z ~ uoa(cid1, cid2), data = s2)
  mod3 <- lm(y ~ z, data = s2, weights = ate(des3))
  expect_identical(mod1$weights, mod3$weights)
  # binary
  s2$z <- s2$z == 1
  des4 <- RCT_Design(z ~ uoa(cid1, cid2), data = s2)
  mod4 <- lm(y ~ z, data = s2, weights = ate(des4))
  expect_identical(mod1$weights, mod4$weights)

  # ordinal treatment
  do1 <- RCT_Design(o ~ uoa(cid1, cid2), data = s2)
  # issue18 - temporary tests and warning
  expect_warning(mo1 <- lm(y ~ o, data = s2, weights = ate(do1)))
  expect_true(all(mo1$weights == 1))
  # end issue18

  s3 <- simdata
  # repeat for ett
  # numeric
  des1 <- RCT_Design(z ~ uoa(cid1, cid2), data = s3)
  mod1 <- lm(y ~ z, data = s3, weights = ett(des1))
  # factor
  s3$z <- as.factor(s3$z)
  des2 <- RCT_Design(z ~ uoa(cid1, cid2), data = s3)
  mod2 <- lm(y ~ z, data = s3, weights = ett(des2))
  expect_identical(mod1$weights, mod2$weights)
  # ordinal
  s3$z <- as.ordered(s3$z)
  des3 <- RCT_Design(z ~ uoa(cid1, cid2), data = s3)
  mod3 <- lm(y ~ z, data = s3, weights = ett(des3))
  expect_identical(mod1$weights, mod3$weights)
  # binary
  s3$z <- s3$z == 1
  des4 <- RCT_Design(z ~ uoa(cid1, cid2), data = s3)
  mod4 <- lm(y ~ z, data = s3, weights = ett(des4))
  expect_identical(mod1$weights, mod4$weights)

  print(4)
  # ordinal treatment
  do1 <- RCT_Design(o ~ uoa(cid1, cid2), data = s3)
  # issue18 - temporary tests and warning
  expect_warning(mo1 <- lm(y ~ o, data = s3, weights = ett(do1)))
  expect_true(all(mo1$weights == 1))
  # end issue18

})
