test_that("DirectAdjusted object", {

  data(simdata)
  des <- Obs_Design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- DirectAdjusted(lm(y ~ x, data = simdata, weights = ate(des)),
                         Design = des, target = "ett")

  expect_s4_class(dalm, "DirectAdjusted")
  expect_true(is(dalm, "lm"))

  expect_identical(dalm$model$"(weights)"@Design, des)
  expect_identical(dalm$model$"(weights)"@Design, dalm@Design)

  expect_error(DirectAdjusted(lm(y ~ x, data = simdata, weights = ate(des)),
                              Design = des, target = "abc"),
               "must be one of")
})

test_that("DirectAdjusted print/show", {

  data(simdata)
  des <- Obs_Design(z ~ cluster(cid2, cid1) + block(bid), data = simdata)

  dalm <- DirectAdjusted(lm(y ~ x, data = simdata, weights = ate(des)),
                         Design = des, target = "ett")

  aslm <- as(dalm, "lm")

  expect_silent(invisible(capture.output(expect_identical(print(dalm), dalm))))
  expect_silent(invisible(capture.output(expect_identical(show(dalm), dalm))))

  expect_output(print(dalm), "Coeff")
  expect_output(show(dalm), "Coeff")
})
