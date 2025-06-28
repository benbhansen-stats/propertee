test_that("internal weight function", {
  data(simdata)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  wspec <- propertee:::.weights_calc(spec, data = simdata, by = NULL, target = "ate", weightAlias = "ate",
                        dichotomy = NULL)
  expect_s4_class(wspec, "WeightedStudySpecification")
  expect_true(is.numeric(wspec@.Data))
  expect_s4_class(wspec@StudySpecification, "StudySpecification")
  expect_identical(spec, wspec@StudySpecification)
  expect_type(wspec@target, "character")

  expect_equal(wspec@target, "ate")

  expect_equal(nrow(simdata), length(wspec))
  expect_true(all(wspec == wspec@.Data))

  expect_identical(deparse1(stats::formula()), deparse1(wspec@dichotomy))

  wspec <- propertee:::.weights_calc(spec, data = simdata, by = NULL, target = "ett", weightAlias = "ett",
                        dichotomy = NULL)
  expect_s4_class(wspec, "WeightedStudySpecification")
  expect_true(is.numeric(wspec@.Data))
  expect_s4_class(wspec@StudySpecification, "StudySpecification")
  expect_identical(spec, wspec@StudySpecification)
  expect_type(wspec@target, "character")

  expect_equal(wspec@target, "ett")

  expect_equal(nrow(simdata), length(wspec))
  expect_true(all(wspec == wspec@.Data))

  expect_identical(deparse1(stats::formula()), deparse1(wspec@dichotomy))

  wspec <- propertee:::.weights_calc(spec, data = simdata, by = NULL, target = "ato", weightAlias = "ato",
                        dichotomy = NULL)
  expect_s4_class(wspec, "WeightedStudySpecification")
  expect_true(is.numeric(wspec@.Data))
  expect_s4_class(wspec@StudySpecification, "StudySpecification")
  expect_identical(spec, wspec@StudySpecification)
  expect_type(wspec@target, "character")

  expect_equal(wspec@target, "ato")

  expect_equal(nrow(simdata), length(wspec))
  expect_true(all(wspec == wspec@.Data))

  expect_identical(deparse1(stats::formula()), deparse1(wspec@dichotomy))

  expect_error(propertee:::.weights_calc(spec, data = simdata, by = NULL,
                             target = "foo", dichotomy = NULL),
               "Invalid weight target")

})

test_that("dichotomy issues", {

  data(simdata)
  spec <- rct_spec(dose ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_error(propertee:::.weights_calc(spec, data = simdata, by = NULL,
                             target = "ate", weightAlias = "ate", dichotomy = NULL),
               "Must provide a dichotomy")

  wspec <- propertee:::.weights_calc(spec, data = simdata, by = NULL, target = "ate", weightAlias = "ate",
                        dichotomy = . ~ dose > 150)
  expect_s4_class(wspec, "WeightedStudySpecification")
  expect_true(is.numeric(wspec@.Data))
  expect_s4_class(wspec@StudySpecification, "StudySpecification")
  expect_identical(spec, wspec@StudySpecification)
  expect_type(wspec@target, "character")

  expect_equal(wspec@target, "ate")

  expect_equal(nrow(simdata), length(wspec))
  expect_true(all(wspec == wspec@.Data))

  expect_identical(deparse1(. ~ dose > 150), deparse1(wspec@dichotomy))
})

test_that("internal and external weight function agreement", {
  data(simdata)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  iwspec <- propertee:::.weights_calc(spec, data = simdata, by = NULL,
                         target = "ate", weightAlias = "ate", dichotomy = NULL)
  ewspec <- ate(spec, data = simdata)
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(iwspec, slot)),
                           deparse1(methods::slot(ewspec, slot)))
                } else {
                  identical(methods::slot(iwspec, slot), methods::slot(ewspec, slot))
                },
               logical(1L)))
  )

  iwspec <- propertee:::.weights_calc(spec, data = simdata, by = NULL,
                         target = "ett", weightAlias = "ett", dichotomy = NULL)
  ewspec <- ett(spec, data = simdata)
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(iwspec, slot)),
                           deparse1(methods::slot(ewspec, slot)))
               } else {
                 identical(methods::slot(iwspec, slot), methods::slot(ewspec, slot))
               },
               logical(1L)))
  )

})

test_that("ate and ett with data argument", {

  # n_clusters = n
  data(mtcars)
  mtcars <- mtcars[-c(5, 11), ]
  spec <- rd_spec(am ~ cluster(qsec) + forcing(vs), data = mtcars)

  wspec <- ate(spec, data = mtcars)

  expect_s4_class(wspec, "WeightedStudySpecification")
  expect_true(is.numeric(wspec@.Data))
  expect_s4_class(wspec@StudySpecification, "StudySpecification")
  expect_identical(spec, wspec@StudySpecification)
  expect_type(wspec@target, "character")

  expect_equal(wspec@target, "ate")

  expect_equal(nrow(mtcars), length(wspec))
  expect_true(all(wspec == wspec@.Data))

  expect_identical(deparse1(stats::formula()), deparse1(wspec@dichotomy))

  # n_clusters < n
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  wspec <- ett(spec, data = simdata)

  expect_s4_class(wspec, "WeightedStudySpecification")
  expect_true(is.numeric(wspec@.Data))
  expect_s4_class(wspec@StudySpecification, "StudySpecification")
  expect_identical(spec, wspec@StudySpecification)
  expect_type(wspec@target, "character")

  expect_equal(wspec@target, "ett")

  expect_equal(nrow(simdata), length(wspec))
  expect_true(all(wspec == wspec@.Data))

  expect_identical(deparse1(stats::formula()), deparse1(wspec@dichotomy))

})

test_that("ate, ett and ato in lm call", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata, weights = ate(spec))

  expect_equal(mod$weights@.Data, ate(spec, data = simdata)@.Data)

  mod <- lm(y ~ x, data = simdata, weights = ett(spec))

  expect_equal(mod$weights@.Data, ett(spec, data = simdata)@.Data)

  mod <- lm(y ~ x, data = simdata, weights = ato(spec))

  expect_equal(mod$weights@.Data, ato(spec, data = simdata)@.Data)

})

test_that("by", {

  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  w1 <- ate(spec, data = simdata)

  simdata2 <- simdata
  colnames(simdata2)[2] <- "cccc"
  w2 <- ate(spec, data = simdata2, by = c("uoa2" = "cccc"))
  w3 <- ate(spec, data = simdata2, by = list("uoa2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)
  expect_equal(w1@.Data, w3@.Data)

  w1 <- ett(spec, data = simdata)
  w2 <- ett(spec, data = simdata2, by = c("uoa2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)

  mod1 <- lm(y ~ x, data = simdata, weights = ate(spec))
  mod2 <- lm(y ~ x, data = simdata2,
             weights = ate(spec, by = c("uoa2" = "cccc")))
  mod3 <- lm(y ~ x, data = simdata2,
             weights = ate(spec, by = list("uoa2" = "cccc")))

  expect_identical(mod1$coef, mod2$coef)
  expect_identical(mod1$coef, mod3$coef)

  mod1 <- lm(y ~ x, data = simdata, weights = ett(spec))
  mod2 <- lm(y ~ x, data = simdata2,
             weights = ett(spec, by = c("uoa2" = "cccc")))

  expect_identical(mod1$coef, mod2$coef)

  spec <- rd_spec(z ~ unit_of_assignment(uoa2, uoa1) +
                     block(bid) + forcing(force), data = simdata)

  w1 <- ate(spec, data = simdata)

  simdata2 <- simdata
  colnames(simdata2)[2] <- "cccc"
  w2 <- ate(spec, data = simdata2, by = c("uoa2" = "cccc"))

  expect_equal(w1@.Data, w2@.Data)

  expect_error(ate(spec, data = simdata2, by = c("uoa2" = "z")),
               "variables cannot be used more than once")

  expect_error(ate(spec, data = simdata2, by = c("abc")), "named vector")

  expect_error(ate(spec, data = simdata2,
                           by = c("uoa2" = "cccc", "uoa2" = "bbbb")),
               "unique")

  expect_warning(w3 <- ate(spec, data = simdata2,
                           by = c("uoa2" = "cccc", "abc" = "uoa1")),
               "not found in StudySpecification")

  expect_warning(w4 <- ate(spec, data = simdata2,
                           by = c("uoa2" = "cccc", "uoa1" = "abc")),
               "not found in data")
  expect_identical(w3, w4)

  expect_warning(expect_warning(w5 <- ate(spec, data = simdata2,
                           by = c("uoa2" = "cccc", "abc" = "def")),
               "not found in data"), "not found in StudySpecification")
  expect_identical(w3, w5)
})

test_that("Ops", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  wspec <- ett(spec, data = simdata)

  expect_s4_class(wspec * 2, "WeightedStudySpecification")
  expect_s4_class(2 * wspec, "WeightedStudySpecification")
  expect_identical(wspec * 3, 3 * wspec)
  expect_equal((wspec * 3)@.Data, wspec@.Data * 3)

  expect_s4_class(wspec / 2, "WeightedStudySpecification")
  expect_s4_class(2 / wspec, "WeightedStudySpecification")
  expect_equal((wspec / 3)@.Data, wspec@.Data / 3)

  expect_error(1 + wspec)
  expect_error(wspec + 1)

  expect_error(1 - wspec)
  expect_error(wspec - 1)

  expect_s4_class(wspec ^ 3, "WeightedStudySpecification")
  expect_equal((wspec ^ 3)@.Data, wspec@.Data ^ 3)


})

test_that("show WeightedStudySpecification", {

  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  wspec <- ett(spec, data = simdata)

  expect_silent(invisible(capture.output(expect_identical(wspec, show(wspec)))))
})
test_that("weight function", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  wspec <- ett(spec, data = simdata)

  expect_s4_class(wspec, "WeightedStudySpecification")
  expect_false(inherits(weights(wspec), "WeightedStudySpecification"))

})

test_that("ett/ate/ato treatment weights are correct length", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(spec, data = testdata)

  expect_equal(length(wts), nrow(testdata))

  wts <- ate(spec, data = testdata)

  expect_equal(length(wts), nrow(testdata))

  wts <- ato(spec, data = testdata)

  expect_equal(length(wts), nrow(testdata))

})

test_that("ett/ate/ato treatment weights return numeric", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(spec, data = testdata)

  expect_s4_class(wts, "numeric")

  wts <- ate(spec, data = testdata)

  expect_s4_class(wts, "numeric")

  wts <- ato(spec, data = testdata)

  expect_s4_class(wts, "numeric")

})

test_that("ett treatment weights = 1", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 4), rep(0, 6)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(spec, data = testdata)

  expect_equal(wts@.Data[1:4], rep(1, 4))

})

test_that("ett weights for block with P(Z = 1) = 0.5: 1", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(spec, data = testdata)

  expect_equal(wts@.Data, rep(1, 10))

})

test_that("ett weights for block with P(Z = 1) = 1/3:  1 or 0.5", {

  testdata <- data.frame(cid = 1:30, z = c(rep(1, 10), rep(0, 20)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ett(spec, data = testdata)

  expect_equal(wts@.Data, c(rep(1, 10), rep(0.5, 20)))

})

test_that("ate weights for block with P(Z = 1) = 0.5: 2", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ate(spec, data = testdata)

  expect_equal(wts@.Data, rep(2, 10))

})

test_that("ato weights for block with P(Z = 1) = 0.5: 0.5", {

  testdata <- data.frame(cid = 1:10, z = c(rep(1, 5), rep(0, 5)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ato(spec, data = testdata)

  expect_equal(wts@.Data, rep(0.5, 10))

})


test_that("ate weights for block with P(Z = 1) = 1/3:  3 or 1.5", {

  testdata <- data.frame(cid = 1:30, z = c(rep(1, 10), rep(0, 20)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ate(spec, data = testdata)

  expect_equal(wts@.Data, c(rep(3, 10), rep(1.5, 20)))

})

test_that("ato weights for block with P(Z = 1) = 1/3:  2/3 or 1/3", {

  testdata <- data.frame(cid = 1:30, z = c(rep(1, 10), rep(0, 20)))

  spec <- rct_spec(z ~ unit_of_assignment(cid), data = testdata)

  wts <- ato(spec, data = testdata)

  expect_equal(wts@.Data, c(rep(2/3, 10), rep(1/3, 20)))

})


test_that("combine two previous blocks, obtain proper weightsfor ett and ate", {

  testdata <- data.frame(bid = c(rep(1, 10), rep(2, 30)),
                         cid = 1:40,
                         z = c(rep(1, 5), rep(0, 5),
                               rep(1, 10), rep(0, 20)))

  spec <- rct_spec(z ~ unit_of_assignment(cid) + block(bid), data = testdata)

  wts <- ett(spec, data = testdata)

  expect_equal(wts@.Data, c(rep(1, 20), rep(0.5, 20)))

  wts <- ate(spec, data = testdata)

  expect_equal(wts@.Data, c(rep(2, 10), rep(3, 10), rep(1.5, 20)))

  wts <- ato(spec, data = testdata)

  expect_equal(wts@.Data, c(rep(0.5, 10), rep(2/3, 10), rep(1/3, 20)))

})

test_that("formula as object is able to be found in data search", {
  data(simdata)

  spec1 <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ x, data = simdata, weights = ate(spec1))

  f <- z ~ uoa(uoa1, uoa2)
  spec2 <- rct_spec(f, data = simdata)
  mod2 <- lm(y ~ x, data = simdata, weights = ate(spec2))

  expect_true(all(mod1$coef == mod2$coef))


})

test_that("numeric and logical treatments work the same", {
  data(simdata)

  # ways of representing binary treatments
  # numeric
  spec1 <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, data = simdata, weights = ate(spec1))
  # binary
  simdata$z <- simdata$z == 1
  spec2 <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  mod2 <- lm(y ~ z, data = simdata, weights = ate(spec2))
  expect_identical(mod1$weights, mod2$weights)
  # repeat for ett
  # numeric
  spec1 <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, data = simdata, weights = ett(spec1))
  # binary
  simdata$z <- simdata$z == 1
  spec2 <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  mod2 <- lm(y ~ z, data = simdata, weights = ett(spec2))
  expect_identical(mod1$weights, mod2$weights)

})

test_that("Weighting respects ordering", {
  data(simdata)
  sd <- simdata[c(1, 11, 31), ]
  # Rows 1/2 are z=0, row 3 is z=1

  spec <- rct_spec(z ~ cluster(uoa1), data = sd)
  w1 <- ate(spec, data = sd)
  # Since first two rows are the same, they should have the same weight
  expect_true(w1[1] == w1[2])
  expect_true(w1[1] != w1[3])

  # Re-order
  sd2 <- sd[c(1, 3, 2), ]
  # Now rows 1 and 3 share treatment
  w2 <- ate(spec, data = sd2)
  expect_true(w2[1] == w2[3])
  expect_true(w2[1] != w2[2])

})

test_that("Passing previously created formula", {
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, data = simdata, weights = ett(spec))
  f <- y ~ z
  mod2 <- lm(f, data = simdata, weights = ett(spec))
  mod1$call <- mod2$call <- call("ls")
  expect_identical(mod1, mod2)
})

test_that("subsetting of WeightedStudySpecification", {
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  ww <- ate(spec, data = simdata)

  expect_identical(subset(ww, rep(TRUE, 50)), ww)
  expect_identical(subset(ww, c(rep(TRUE, 10), rep(FALSE, 40)))@.Data,
                   ww@.Data[1:10])

  expect_identical(ww[], ww)
  expect_identical(ww[1:10]@.Data, ww@.Data[1:10])
})

test_that("#130 zero weights with non-varying treatment in a block", {
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  wts <- lmitt(y ~ 1, data = simdata,
               specification = spec, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  wts <- lmitt(y ~ 1, data = simdata,
               specification = spec, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  # Make a block with no controls.
  simdata$z[simdata$bid == 3] <- 1
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  wts <- lmitt(y ~ 1, data = simdata,
               specification = spec, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts[simdata$bid != 3] > 0))
  expect_true(all(wts[simdata$bid == 3] == 0))

  wts <- lmitt(y ~ 1, data = simdata,
               specification = spec, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts[simdata$bid != 3] > 0))
  expect_true(all(wts[simdata$bid == 3] == 0))

  # Make a block with no treatment
  simdata$z[simdata$bid == 3] <- 0
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  wts <- lmitt(y ~ 1, data = simdata,
               specification = spec, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts[simdata$bid != 3] > 0))
  expect_true(all(wts[simdata$bid == 3] == 0))

  wts <- lmitt(y ~ 1, data = simdata,
               specification = spec, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts[simdata$bid != 3] > 0))
  expect_true(all(wts[simdata$bid == 3] == 0))


  #### Data with an NA block doesn't break
  if (requireNamespace("AER", quietly = TRUE)) {
      data(STARdata)
      STARdata <- STAR
  STARdata$starkbinary <- STARdata$stark == "small"
  STARdata$studentid <- seq_len(nrow(STARdata))
  spec <- obs_spec(starkbinary ~ unit_of_assignment(studentid) +
                      block(schoolidk),
                    data = STARdata, na.fail = FALSE)
  # Some students in STARdata do not have `schoolidk` and thus NA for a block

  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ato")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(all(wts > 0))

  # Make a block with no controls.
  STARdata$starkbinary[STARdata$schoolidk == 1] <- 1
  spec <- obs_spec(starkbinary ~ unit_of_assignment(studentid) +
                      block(schoolidk),
                    data = STARdata, na.fail = FALSE)
  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ato")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

  # Make a block with no treatment
  STARdata$starkbinary[STARdata$schoolidk == 1] <- 0
  spec <- obs_spec(starkbinary ~ unit_of_assignment(studentid) +
                      block(schoolidk),
                    data = STARdata, na.fail = FALSE)
  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ate")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ett")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)

  wts <- lmitt(readk ~ 1, data = STARdata,
               specification = spec, weights = "ato")$model$"(weights)"
  expect_true(!any(is.nan(wts)))
  expect_true(sum(wts == 0) > 0)
}
  
})

test_that("#131 numeric blocks don't cause NA weights", {
  data(simdata)
  simdata$bid[simdata$bid == 3] <- 4

  spec <- obs_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  mod <-lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE, weights = "ate")
  expect_length(mod$weights, nrow(simdata))
  expect_true(all(!is.na(mod$weights)))


})

test_that(paste("weights with attention to blocks when `data` has different order",
                "from specification data"), {
  set.seed(31)
  specification_dat <- data.frame(id = letters[seq_len(10)],
           z = c(rep(c(1, 1, 0), 2), rep(c(0, 1), 2)),
           blk = c(rep(LETTERS[1:2], each = 3),
                   rep(LETTERS[3:4], each = 2)))
  spec <- rct_spec(z ~ unitid(id) + block(blk), specification_dat)

  analysis_dat <- rbind(specification_dat[specification_dat$blk == "B",],
                        specification_dat[specification_dat$blk == "D",],
                        specification_dat[specification_dat$blk == "A",],
                        specification_dat[specification_dat$blk == "C",])

  wts <- .weights_calc(spec, target = "ate", weightAlias = "ate", by = NULL, dichotomy = NULL, data = analysis_dat)
  expected_triplet_wts <- c(rep(3/2, 2), 3)
  expected_pair_wts <- rep(2, 2)
  expect_equal(wts@.Data, c(expected_triplet_wts, expected_pair_wts,
                            expected_triplet_wts, expected_pair_wts))
})

test_that("#180 non-exhaustive dichotomies", {
  # no missing block ID's
  data(simdata)
  spec <- rct_spec(dose ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  wspec  <- .weights_calc(spec, data = simdata, by = NULL, target = "ate", weightAlias = "ate",
                         dichotomy = dose >200 ~ dose <200)

  expect_true(length(wspec) == nrow(simdata))
  expect_true(all.equal(which(wspec == 0), which(simdata$dose == 200)))
  expect_true(all.equal(which(wspec == 0),
                        which(is.na(.bin_txt(spec, simdata, dose > 200 ~ dose<200)))))
  expect_equal(
    nrow(lmitt(y~1, specification = spec, data = simdata, weights = wspec,
               dichotomy = dose >200 ~ dose <200)$model),
    nrow(simdata[simdata$dose != 200,])
  )

  # Ben's tests
  wspec2  <- propertee:::.weights_calc(spec, data = simdata, by = NULL, target = "ate", weightAlias = "ate", dichotomy = dose >200 ~ dose <200)
  expect_true(all(wspec2[simdata$dose==200]==0))
  expect_true(all(wspec2[simdata$dose!=200]!=0))

  wspec3  <- propertee:::.weights_calc(spec, data = simdata, by = NULL, target = "ett", weightAlias = "ett", dichotomy = dose >200 ~ dose <100)
  expect_true(all(wspec3[simdata$bid==3]==0))

  # missing block ID's
  simdata[simdata$uoa1 == 1 & simdata$uoa2 == 1, "bid"] <- NA_integer_
  spec <- rct_spec(dose ~ uoa(uoa1, uoa2) + block(bid), data = simdata, na.fail = FALSE)
  wspec4  <- .weights_calc(spec, data = simdata, by = NULL, target = "ate", weightAlias = "ate",
                          dichotomy = dose >200 ~ dose <200)

  expect_true(length(wspec4) == nrow(simdata))
  expect_true(all.equal(which(wspec4 == 0), which(simdata$dose == 200 | is.na(simdata$bid))))
  expect_true(all.equal(which(wspec4 == 0),
                        which(is.na(.bin_txt(spec, simdata, dose > 200 ~ dose<200)))))
  expect_equal(
    nrow(lmitt(y~1, specification = spec, data = simdata, weights = wspec4,
               dichotomy = dose >200 ~ dose <200)$model),
    nrow(simdata[simdata$dose != 200 & !is.na(simdata$bid),])
  )

  # missing unit ID's
  data(simdata)
  simdata[simdata$uoa1 == 1 & simdata$uoa2 == 1, paste0("uoa", c(1, 2))] <- NA_integer_
  spec <- rct_spec(dose ~ uoa(uoa1, uoa2) + block(bid), data = simdata, na.fail = FALSE)
  wspec5  <- .weights_calc(spec, data = simdata, by = NULL, target = "ate", weightAlias = "ate",
                          dichotomy = dose >200 ~ dose <200)

  expect_true(length(wspec5) == nrow(simdata))
  expect_true(all.equal(which(wspec5 == 0), which(simdata$dose == 200 | is.na(simdata$uoa1))))
  expect_true(all.equal(which(wspec5 == 0),
                        which(is.na(.bin_txt(spec, simdata, dose > 200 ~ dose<200)))))
  expect_equal(
    nrow(lmitt(y~1, specification = spec, data = simdata, weights = wspec4,
               dichotomy = dose >200 ~ dose <200)$model),
    nrow(simdata[simdata$dose != 200 & !is.na(simdata$uoa1),])
  )
})
