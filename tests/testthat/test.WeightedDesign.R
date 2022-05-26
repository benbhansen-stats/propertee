test_that("internal weight function", {
  data(simdata)
  des <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)

  wdes <- .weights_calc(des, clusterdata = simdata, by = NULL, target = "ate",
                        dichotomy = NULL)
  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ate")

  expect_equal(nrow(simdata), length(wdes))
  expect_true(all(wdes == wdes@.Data))

  wdes <- .weights_calc(des, clusterdata = simdata, by = NULL, target = "ett",
                        dichotomy = NULL)
  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ett")

  expect_equal(nrow(simdata), length(wdes))
  expect_true(all(wdes == wdes@.Data))

  expect_error(.weights_calc(des, clusterdata = simdata, by = NULL,
                             target = "foo", dichotomy = NULL),
               "Invalid weight target")

  expect_error(.weights_calc(des, clusterdata = 1, by = NULL, target = "ate",
                             dichotomy = NULL),
               "`data` must be")

  expect_error(.weights_calc(des, clusterdata = simdata, by = NULL,
                             target = "ate", dichotomy = 1),
               "`dichotomy` must be")

})

test_that("dichotomy issues", {

  data(simdata)
  des <- rct_design(dose ~ uoa(cid1, cid2) + block(bid), data = simdata)

  expect_error(.weights_calc(des, clusterdata = simdata, by = NULL,
                             target = "ate", dichotomy = NULL),
               "must have a dichotomy")

  dichotomy(des) <- . ~ dose > 150

  wdes <- .weights_calc(des, clusterdata = simdata, by = NULL, target = "ate",
                        dichotomy = NULL)
  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ate")

  expect_equal(nrow(simdata), length(wdes))
  expect_true(all(wdes == wdes@.Data))

  expect_warning(wdes <- .weights_calc(des, clusterdata = simdata, by = NULL,
                                       target = "ate",
                                       dichotomy = dose > 200 ~ .),
                 "over-writing")
  expect_identical(deparse(dichotomy(wdes@Design)), "dose > 200 ~ .")

})

test_that("internal and external weight function agreement", {
  data(simdata)
  des <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)

  iwdes <- .weights_calc(des, clusterdata = simdata, by = NULL,
                         target = "ate", dichotomy = NULL)
  ewdes <- ate(des, data = simdata)
  expect_identical(iwdes, ewdes)

  iwdes <- .weights_calc(des, clusterdata = simdata, by = NULL,
                         target = "ett", dichotomy = NULL)
  ewdes <- ett(des, data = simdata)
  expect_identical(iwdes, ewdes)

})

test_that("ate and ett with data argument", {

  # n_clusters = n
  data(mtcars)
  mtcars <- mtcars[-c(5, 11), ]
  des <- rd_design(am ~ cluster(qsec) + forcing(vs), data = mtcars)

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
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

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
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata, weights = ate(des))

  expect_equal(mod$weights, ate(des, data = simdata)@.Data)

  mod <- lm(y ~ x, data = simdata, weights = ett(des))

  expect_equal(mod$weights, ett(des, data = simdata)@.Data)

})

test_that("by", {

  data(simdata)
  des <- rd_design(z ~ cluster(cid2, cid1) + block(bid) + forcing(force),
                   data = simdata)

  w1 <- ate(des, data = simdata)

  simdata2 <- simdata
  colnames(simdata2)[2] <- "cccc"
  w2 <- ate(des, data = simdata2, by = c("cid2" = "cccc"))
  w3 <- ate(des, data = simdata2, by = list("cid2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)
  expect_equal(w1@.Data, w3@.Data)

  w1 <- ett(des, data = simdata)
  w2 <- ett(des, data = simdata2, by = c("cid2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)

  mod1 <- lm(y ~ x, data = simdata, weights = ate(des))
  mod2 <- lm(y ~ x, data = simdata2,
             weights = ate(des, by = c("cid2" = "cccc")))
  mod3 <- lm(y ~ x, data = simdata2,
             weights = ate(des, by = list("cid2" = "cccc")))

  expect_identical(mod1$coef, mod2$coef)
  expect_identical(mod1$coef, mod3$coef)

  mod1 <- lm(y ~ x, data = simdata, weights = ett(des))
  mod2 <- lm(y ~ x, data = simdata2,
             weights = ett(des, by = c("cid2" = "cccc")))

  expect_identical(mod1$coef, mod2$coef)

  des <- rd_design(z ~ unit_of_assignment(cid2, cid1) +
                     block(bid) + forcing(force), data = simdata)

  w1 <- ate(des, data = simdata)

  simdata2 <- simdata
  colnames(simdata2)[2] <- "cccc"
  w2 <- ate(des, data = simdata2, by = c("cid2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)

  expect_error(ate(des, data = simdata2, by = c("cid2" = "z")),
               "variables cannot be used more than once")

})

test_that("Ops", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

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
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_silent(invisible(capture.output(expect_identical(wdes, show(wdes)))))
})
test_that("weight function", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_s4_class(wdes, "WeightedDesign")
  expect_false(is(weights(wdes), "WeightedDesign"))

  })

test_that("ett/ate treatment weights are correct length", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  des <- rct_design(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(length(wts), nrow(testdata))

  wts <- ate(des, data = testdata)

  expect_equal(length(wts), nrow(testdata))

})

test_that("ett/ate treatment weights return numeric", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  des <- rct_design(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_s4_class(wts, "numeric")

  wts <- ate(des, data = testdata)

  expect_s4_class(wts, "numeric")

})

test_that("ett treatment weights = 1", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  des <- rct_design(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(wts@.Data[1:4], rep(1, 4))

})

test_that("ett weights for block with P(Z = 1) = 0.5: 1", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)))

  des <- rct_design(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(wts@.Data, rep(1, 10))

})

test_that("ett weights for block with P(Z = 1) = 1/3:  1 or 0.5", {

  testdata <- data.frame(cid = 1:30, z = c(rep(1, 10), rep(0, 20)))

  des <- rct_design(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(wts@.Data, c(rep(1, 10), rep(0.5, 20)))

})

test_that("ate weights for block with P(Z = 1) = 0.5: 2", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)))

  des <- rct_design(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ate(des, data = testdata)

  expect_equal(wts@.Data, rep(2, 10))

})

test_that("ate weights for block with P(Z = 1) = 1/3:  3 or 1.5", {

  testdata <- data.frame(cid = 1:30, z = c(rep(1, 10), rep(0, 20)))

  des <- rct_design(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ate(des, data = testdata)

  expect_equal(wts@.Data, c(rep(3, 10), rep(1.5, 20)))

})

test_that("combine two previous blocks, obtain proper weightsfor ett and ate", {

  testdata <- data.frame(bid = c(rep(1, 10), rep(2, 30)),
                         cid = 1:40,
                         z = c(rep(1, 5), rep(0, 5),
                               rep(1, 10), rep(0, 20)))

  des <- rct_design(z ~ unit_of_assignment(cid) + block(bid), data = testdata)

  wts <- ett(des, data = testdata)

  expect_equal(wts@.Data, c(rep(1, 20), rep(0.5, 20)))

  wts <- ate(des, data = testdata)

  expect_equal(wts@.Data, c(rep(2, 10), rep(3, 10), rep(1.5, 20)))

})

test_that("formula as object is able to be found in data search", {
  data(simdata)

  des1 <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  mod1 <- lm(y ~ x, data = simdata, weights = ate(des1))

  f <- z ~ uoa(cid1, cid2)
  des2 <- rct_design(f, data = simdata)
  mod2 <- lm(y ~ x, data = simdata, weights = ate(des2))

  expect_true(all(mod1$coef == mod2$coef))


})

test_that("numeric and logical treatments work the same", {
  data(simdata)

  # ways of representing binary treatments
  # numeric
  des1 <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  mod1 <- lm(y ~ z, data = simdata, weights = ate(des1))
  # binary
  simdata$z <- simdata$z == 1
  des2 <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  mod2 <- lm(y ~ z, data = simdata, weights = ate(des2))
  expect_identical(mod1$weights, mod2$weights)
  # repeat for ett
  # numeric
  des1 <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  mod1 <- lm(y ~ z, data = simdata, weights = ett(des1))
  # binary
  simdata$z <- simdata$z == 1
  des2 <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  mod2 <- lm(y ~ z, data = simdata, weights = ett(des2))
  expect_identical(mod1$weights, mod2$weights)

})

test_that("Weighting respects ordering", {
  data(simdata)
  sd <- simdata[c(1, 11, 31), ]
  # Rows 1/2 are z=0, row 3 is z=1

  des <- rct_design(z ~ cluster(cid1), data = sd)
  w1 <- ate(des, data = sd)
  # Since first two rows are the same, they should have the same weight
  expect_true(w1[1] == w1[2])
  expect_true(w1[1] != w1[3])

  # Re-order
  sd2 <- sd[c(1, 3, 2), ]
  # Now rows 1 and 3 share treatment
  w2 <- ate(des, data = sd2)
  expect_true(w2[1] == w2[3])
  expect_true(w2[1] != w2[2])

})

test_that("Combining weighted designs", {
  data(simdata)

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)

  w1 <- ate(des, data = simdata[1:30,])
  w2 <- ate(des, data = simdata[31:40,])
  w3 <- ate(des, data = simdata[41:50,])

  c_w <- c(w1, w2, w3)
  expect_true(is(c_w, "WeightedDesign"))
  expect_length(c_w, 50)
  expect_identical(c_w, ate(des, data = simdata))

  w1e <- ett(des, data = simdata[1:30,])
  w2e <- ett(des, data = simdata[31:40,])
  w3e <- ett(des, data = simdata[41:50,])

  c_we <- c(w1e, w2e, w3e)
  expect_true(is(c_we, "WeightedDesign"))
  expect_length(c_we, 50)
  expect_identical(c_we, ett(des, data = simdata))

  expect_error(c(w1, 1:5), "with other")
  expect_error(c(w1, w1e), "same target")

  des2 <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)

  alt_w1 <- ate(des2, data = simdata)

  expect_error(c(w1, alt_w1), "which differ on elements")

  # if the first argument is compatible with WeightedDesign but isn't one (e.g.
  # numeric vector), c() will return a numeric vector
  #expect_true(is(c(1:5, w1), "WeightedDesign"))
})

test_that("Combining weighted designs with different dichotomys ", {
  des <- rct_design(dose ~ uoa(cid1, cid2), data = simdata)

  w1 <- ate(des, data = simdata, dichotomy = dose >= 300 ~ .)
  w2 <- ate(des, data = simdata, dichotomy = dose >= 200 ~ .)
  w3 <- ate(des, data = simdata, dichotomy = dose >= 100 ~ .)

  c_w <- c(w1, w2, w3)
  expect_true(is(c_w, "WeightedDesign"))
  expect_length(c_w, 150)

  expect_error(c(w1, w2, w3, force_dichotomy_equal = TRUE),
               "must be identical")
})
