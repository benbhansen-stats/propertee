test_that("internal weight function", {
  data(simdata)
  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  wdes <- propertee:::.weights_calc(des, data = simdata, by = NULL, target = "ate",
                        dichotomy = NULL)
  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ate")

  expect_equal(nrow(simdata), length(wdes))
  expect_true(all(wdes == wdes@.Data))
  
  expect_identical(deparse(stats::formula()), deparse(wdes@dichotomy))

  wdes <- propertee:::.weights_calc(des, data = simdata, by = NULL, target = "ett",
                        dichotomy = NULL)
  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ett")

  expect_equal(nrow(simdata), length(wdes))
  expect_true(all(wdes == wdes@.Data))
  
  expect_identical(deparse(stats::formula()), deparse(wdes@dichotomy))

  expect_error(propertee:::.weights_calc(des, data = simdata, by = NULL,
                             target = "foo", dichotomy = NULL),
               "Invalid weight target")

  expect_error(propertee:::.weights_calc(des, data = 1, by = NULL, target = "ate",
                             dichotomy = NULL),
               "`data` must be")

})

test_that("dichotomy issues", {

  data(simdata)
  des <- rct_design(dose ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_error(propertee:::.weights_calc(des, data = simdata, by = NULL,
                             target = "ate", dichotomy = NULL),
               "Must provide a dichotomy")

  wdes <- propertee:::.weights_calc(des, data = simdata, by = NULL, target = "ate",
                        dichotomy = . ~ dose > 150)
  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ate")

  expect_equal(nrow(simdata), length(wdes))
  expect_true(all(wdes == wdes@.Data))
  
  expect_identical(deparse(. ~ dose > 150), deparse(wdes@dichotomy))
})

test_that("internal and external weight function agreement", {
  data(simdata)
  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  iwdes <- propertee:::.weights_calc(des, data = simdata, by = NULL,
                         target = "ate", dichotomy = NULL)
  ewdes <- ate(des, data = simdata)
  expect_identical(iwdes, ewdes)

  iwdes <- propertee:::.weights_calc(des, data = simdata, by = NULL,
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

  expect_identical(deparse(stats::formula()), deparse(wdes@dichotomy))

  # n_clusters < n
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_s4_class(wdes, "WeightedDesign")
  expect_true(is.numeric(wdes@.Data))
  expect_s4_class(wdes@Design, "Design")
  expect_identical(des, wdes@Design)
  expect_type(wdes@target, "character")

  expect_equal(wdes@target, "ett")

  expect_equal(nrow(simdata), length(wdes))
  expect_true(all(wdes == wdes@.Data))
  
  expect_identical(deparse(stats::formula()), deparse(wdes@dichotomy))

})

test_that("ate and ett in lm call", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata, weights = ate(des))

  expect_equal(mod$weights, ate(des, data = simdata)@.Data)

  mod <- lm(y ~ x, data = simdata, weights = ett(des))

  expect_equal(mod$weights, ett(des, data = simdata)@.Data)

})

test_that("by", {

  data(simdata)
  des <- rd_design(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  w1 <- ate(des, data = simdata)

  simdata2 <- simdata
  colnames(simdata2)[2] <- "cccc"
  w2 <- ate(des, data = simdata2, by = c("uoa2" = "cccc"))
  w3 <- ate(des, data = simdata2, by = list("uoa2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)
  expect_equal(w1@.Data, w3@.Data)

  w1 <- ett(des, data = simdata)
  w2 <- ett(des, data = simdata2, by = c("uoa2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)

  mod1 <- lm(y ~ x, data = simdata, weights = ate(des))
  mod2 <- lm(y ~ x, data = simdata2,
             weights = ate(des, by = c("uoa2" = "cccc")))
  mod3 <- lm(y ~ x, data = simdata2,
             weights = ate(des, by = list("uoa2" = "cccc")))

  expect_identical(mod1$coef, mod2$coef)
  expect_identical(mod1$coef, mod3$coef)

  mod1 <- lm(y ~ x, data = simdata, weights = ett(des))
  mod2 <- lm(y ~ x, data = simdata2,
             weights = ett(des, by = c("uoa2" = "cccc")))

  expect_identical(mod1$coef, mod2$coef)

  des <- rd_design(z ~ unit_of_assignment(uoa2, uoa1) +
                     block(bid) + forcing(force), data = simdata)

  w1 <- ate(des, data = simdata)

  simdata2 <- simdata
  colnames(simdata2)[2] <- "cccc"
  w2 <- ate(des, data = simdata2, by = c("uoa2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)

  expect_error(ate(des, data = simdata2, by = c("uoa2" = "z")),
               "variables cannot be used more than once")

  expect_error(ate(des, data = simdata2, by = c("abc")), "named vector")

  expect_error(ate(des, data = simdata2,
                           by = c("uoa2" = "cccc", "uoa2" = "bbbb")),
               "unique")

  expect_warning(w3 <- ate(des, data = simdata2,
                           by = c("uoa2" = "cccc", "abc" = "uoa1")),
               "not found in Design")

  expect_warning(w4 <- ate(des, data = simdata2,
                           by = c("uoa2" = "cccc", "uoa1" = "abc")),
               "not found in data")
  expect_identical(w3, w4)

  expect_warning(expect_warning(w5 <- ate(des, data = simdata2,
                           by = c("uoa2" = "cccc", "abc" = "def")),
               "not found in data"), "not found in Design")
  expect_identical(w3, w5)
})

test_that("Ops", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_s4_class(wdes * 2, "WeightedDesign")
  expect_s4_class(2 * wdes, "WeightedDesign")
  expect_identical(wdes * 3, 3 * wdes)
  expect_equal((wdes * 3)@.Data, wdes@.Data * 3)

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
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_silent(invisible(capture.output(expect_identical(wdes, show(wdes)))))
})
test_that("weight function", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  wdes <- ett(des, data = simdata)

  expect_s4_class(wdes, "WeightedDesign")
  expect_false(inherits(weights(wdes), "WeightedDesign"))

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

  des1 <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ x, data = simdata, weights = ate(des1))

  f <- z ~ uoa(uoa1, uoa2)
  des2 <- rct_design(f, data = simdata)
  mod2 <- lm(y ~ x, data = simdata, weights = ate(des2))

  expect_true(all(mod1$coef == mod2$coef))


})

test_that("numeric and logical treatments work the same", {
  data(simdata)

  # ways of representing binary treatments
  # numeric
  des1 <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, data = simdata, weights = ate(des1))
  # binary
  simdata$z <- simdata$z == 1
  des2 <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  mod2 <- lm(y ~ z, data = simdata, weights = ate(des2))
  expect_identical(mod1$weights, mod2$weights)
  # repeat for ett
  # numeric
  des1 <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, data = simdata, weights = ett(des1))
  # binary
  simdata$z <- simdata$z == 1
  des2 <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  mod2 <- lm(y ~ z, data = simdata, weights = ett(des2))
  expect_identical(mod1$weights, mod2$weights)

})

test_that("Weighting respects ordering", {
  data(simdata)
  sd <- simdata[c(1, 11, 31), ]
  # Rows 1/2 are z=0, row 3 is z=1

  des <- rct_design(z ~ cluster(uoa1), data = sd)
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

test_that("Passing previously created formula", {
  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, data = simdata, weights = ett(des))
  f <- y ~ z
  mod2 <- lm(f, data = simdata, weights = ett(des))
  mod1$call <- mod2$call <- call("ls")
  expect_identical(mod1, mod2)
})

test_that("subsetting of WeightedDesign", {
  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  ww <- ate(des, data = simdata)

  expect_identical(subset(ww, rep(TRUE, 50)), ww)
  expect_identical(subset(ww, c(rep(TRUE, 10), rep(FALSE, 40)))@.Data,
                   ww@.Data[1:10])

  expect_identical(ww[], ww)
  expect_identical(ww[1:10]@.Data, ww@.Data[1:10])
})

test_that("#130 zero weights with non-varying treatment in a block", {
  data(simdata)

  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  wts <- lmitt(y ~ 1, data = simdata,
               design = des, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  wts <- lmitt(y ~ 1, data = simdata,
               design = des, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  # Make a block with no controls.
  simdata$z[simdata$bid == 3] <- 1
  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  wts <- lmitt(y ~ 1, data = simdata,
               design = des, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts[simdata$bid != 3] > 0))
  expect_true(all(wts[simdata$bid == 3] == 0))

  wts <- lmitt(y ~ 1, data = simdata,
               design = des, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts[simdata$bid != 3] > 0))
  expect_true(all(wts[simdata$bid == 3] == 0))

  # Make a block with no treatment
  simdata$z[simdata$bid == 3] <- 0
  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  wts <- lmitt(y ~ 1, data = simdata,
               design = des, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts[simdata$bid != 3] > 0))
  expect_true(all(wts[simdata$bid == 3] == 0))

  wts <- lmitt(y ~ 1, data = simdata,
               design = des, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts[simdata$bid != 3] > 0))
  expect_true(all(wts[simdata$bid == 3] == 0))


  #### Data with an NA block doesn't break
  data(STARdata)
  STARdata$starkbinary <- STARdata$stark == "small"
  STARdata$studentid <- seq_len(nrow(STARdata))
  des <- obs_design(starkbinary ~ unit_of_assignment(studentid) +
                      block(schoolidk),
                    data = STARdata, na.fail = FALSE)
  # Some students in STARdata do not have `schoolidk` and thus NA for a block

  wts <- lmitt(readk ~ 1, data = STARdata,
               design = des, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  wts <- lmitt(readk ~ 1, data = STARdata,
               design = des, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  # Make a block with no controls.
  STARdata$starkbinary[STARdata$schoolidk == 1] <- 1
  des <- obs_design(starkbinary ~ unit_of_assignment(studentid) +
                      block(schoolidk),
                    data = STARdata, na.fail = FALSE)
  wts <- lmitt(readk ~ 1, data = STARdata,
               design = des, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

  des <- obs_design(starkbinary ~ unit_of_assignment(studentid) +
                      block(schoolidk),
                    data = STARdata, na.fail = FALSE)
  wts <- lmitt(readk ~ 1, data = STARdata,
               design = des, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

  # Make a block with no treatment
  STARdata$starkbinary[STARdata$schoolidk == 1] <- 0
  des <- obs_design(starkbinary ~ unit_of_assignment(studentid) +
                      block(schoolidk),
                    data = STARdata, na.fail = FALSE)
  wts <- lmitt(readk ~ 1, data = STARdata,
               design = des, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

  des <- obs_design(starkbinary ~ unit_of_assignment(studentid) +
                      block(schoolidk),
                    data = STARdata, na.fail = FALSE)
  wts <- lmitt(readk ~ 1, data = STARdata,
               design = des, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

})

test_that("#131 numeric blocks don't cause NA weights", {
  data(simdata)
  simdata$bid[simdata$bid == 3] <- 4

  des <- obs_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  mod <-lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE, weights = "ate")
  expect_length(mod$weights, nrow(simdata))
  expect_true(all(!is.na(mod$weights)))


})

test_that(paste("weights with attention to blocks when `data` has different order",
                "from design data"), {
  set.seed(31)
  design_dat <- data.frame(id = letters[seq_len(10)],
           z = c(rep(c(1, 1, 0), 2), rep(c(0, 1), 2)),
           blk = c(rep(LETTERS[1:2], each = 3),
                   rep(LETTERS[3:4], each = 2)))
  des <- rct_design(z ~ unitid(id) + block(blk), design_dat)
  
  analysis_dat <- rbind(design_dat[design_dat$blk == "B",],
                        design_dat[design_dat$blk == "D",],
                        design_dat[design_dat$blk == "A",],
                        design_dat[design_dat$blk == "C",])
  
  wts <- .weights_calc(des, target = "ate", by = NULL, dichotomy = NULL, data = analysis_dat)
  expected_triplet_wts <- c(rep(3/2, 2), 3)
  expected_pair_wts <- rep(2, 2)
  expect_equal(wts@.Data, c(expected_triplet_wts, expected_pair_wts,
                            expected_triplet_wts, expected_pair_wts))
})
