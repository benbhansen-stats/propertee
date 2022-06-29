# This isn't going to test whether anything is actually correct (see individiual
# tests for that), but instead will run each different variation of `lm` call
# and ensure that no messages, warnings, or errors are produced.

# ******************** IMPORTANT ********************
# If debugging things in here, REMOVE EXPECT_SILENT!!! Running it will suppress
# all output, including inside `browser()`
test_that("binary treatment, in all data", {

  # Treatment exists in data
  data(simdata)

  des <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)
  camod <- lm(y ~ x, data = simdata)

  # Weight alone
  expect_silent(as.lmitt(lm(y ~ z, data = simdata,
                                     weights = ate(des))))

  # Adopters alone
  expect_silent(a <- lm(y ~ adopters(des), data = simdata))
  # DA doesnt look in adopters for Design
  # expect_silent(as.lmitt(a))

  # cov_adj alone
  expect_silent(as.lmitt(lm(y ~ z, data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))
  expect_silent(as.lmitt(lm(y ~ z + offset(cov_adj(camod,
                                                            design = des)),
                   data = simdata), target = "ate"))


  # Weight + adopters
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des))))
  expect_silent(as.lmitt(lm(y ~ adopters(des), data = simdata,
                                     weights = ett(des))))


  # weight + cov_adj
  expect_silent(as.lmitt(lm(y ~ z, data = simdata, weights = ett(des),
                   offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ z, data = simdata, weights = ett(),
                   offset = cov_adj(camod, design = des))))
  expect_silent(as.lmitt(lm(y ~ z, data = simdata, weights = ett(des),
                   offset = cov_adj(camod, design = des))))

  # weight + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ z + offset(cov_adj(camod)), data = simdata,
                   weights = ett(des))))
  # Fails when trying to obtain Design from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ z + offset(cov_adj(camod, design = des)),
  #                      data = simdata, weights = ett())))
  expect_silent(as.lmitt(lm(y ~ z + offset(cov_adj(camod, design = des)),
                   data = simdata, weights = ett(des))))

  # adopters + cov_adj
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))
  expect_silent(as.lmitt(lm(y ~ adopters(des), data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))

  # weights + adopters + cov_adj
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des),
                                     offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(),
                                     offset = cov_adj(camod, design = des))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des),
                                     offset = cov_adj(camod, design = des))))

  # weight + adopter + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ adopters() + offset(cov_adj(camod)),
                                     data = simdata, weights = ett(des))))
  # Fails when trying to obtain Design from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ adopters() + offset(cov_adj(camod, design = des)),
  #                      data = simdata, weights = ett()))
  expect_silent(as.lmitt(lm(y ~ adopters() +
                                       offset(cov_adj(camod, design = des)),
                   data = simdata, weights = ett(des))))

})

test_that("binary treatment, not in data2", {
  # Treatment doesn't exist in data
  data(simdata)

  des <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)
  camod <- lm(y ~ x, data = simdata)
  simdata$z <- NULL

  # Adopters alone
  expect_silent(a <- lm(y ~ adopters(des), data = simdata))
  # DA doesnt look in adopters for Design
  # expect_silent(as.lmitt(a))

  # Weight + adopters
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des))))
  expect_silent(as.lmitt(lm(y ~ adopters(des), data = simdata, weights = ett(des))))

  # adopters + cov_adj
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))
  expect_silent(as.lmitt(lm(y ~ adopters(des), data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))

  # weights + adopters + cov_adj
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des),
                                     offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(),
                                     offset = cov_adj(camod, design = des))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des),
                                     offset = cov_adj(camod, design = des))))

  # weight + adopter + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ adopters() + offset(cov_adj(camod)),
                                     data = simdata, weights = ett(des))))
  # Fails when trying to obtain Design from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ adopters() + offset(cov_adj(camod, design = des)),
  #                      data = simdata, weights = ett()))
  expect_silent(as.lmitt(lm(y ~ adopters() +
                                       offset(cov_adj(camod, design = des)),
                   data = simdata, weights = ett(des))))

})

test_that("non-binary treatment, in all data, dichotomization in design", {


  # Treatment exists in data
  data(simdata)

  des <- rct_design(dose ~ uoa(cid1, cid2) + block(bid), data = simdata,
                    dichotomy = dose >= 200 ~ .)
  camod <- lm(y ~ x, data = simdata)


  # Adopters alone
  expect_silent(a <- lm(y ~ adopters(des), data = simdata))
  # DA doesnt look in adopters for Design
  # expect_silent(as.lmitt(a))

  # Weight + adopters
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des))))
  expect_silent(as.lmitt(lm(y ~ adopters(des), data = simdata,
                                     weights = ett(des))))

  # adopters + cov_adj
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))
  expect_silent(as.lmitt(lm(y ~ adopters(des), data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))

  # weights + adopters + cov_adj
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des),
                                     offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(),
                                     offset = cov_adj(camod, design = des))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des),
                                     offset = cov_adj(camod, design = des))))

  # weight + adopter + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ adopters() + offset(cov_adj(camod)),
                                     data = simdata, weights = ett(des))))
  # Fails when trying to obtain Design from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ adopters() + offset(cov_adj(camod, design = des)),
  #                      data = simdata, weights = ett()))
  expect_silent(as.lmitt(lm(y ~ adopters() +
                                       offset(cov_adj(camod, design = des)),
                   data = simdata, weights = ett(des))))

})

test_that("non-binary treatment, not in data2, dichotomization in design", {


  # Treatment exists in data
  data(simdata)

  des <- rct_design(dose ~ uoa(cid1, cid2) + block(bid), data = simdata,
                    dichotomy = dose >= 200 ~ .)
  camod <- lm(y ~ x, data = simdata)
  simdata$dose <- NULL


  # Adopters alone
  expect_silent(a <- lm(y ~ adopters(des), data = simdata))
  # DA doesnt look in adopters for Design
  # expect_silent(as.lmitt(a))

  # Weight + adopters
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des))))
  expect_silent(as.lmitt(lm(y ~ adopters(des), data = simdata,
                                     weights = ett(des))))

  # adopters + cov_adj
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))
  expect_silent(as.lmitt(lm(y ~ adopters(des), data = simdata,
                   offset = cov_adj(camod, design = des)), target = "ate"))

  # weights + adopters + cov_adj
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des),
                                     offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(),
                                     offset = cov_adj(camod, design = des))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                                     weights = ett(des),
                                     offset = cov_adj(camod, design = des))))

  # weight + adopter + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ adopters() + offset(cov_adj(camod)),
                                     data = simdata, weights = ett(des))))
  # Fails when trying to obtain Design from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ adopters() + offset(cov_adj(camod, design = des)),
  #                      data = simdata, weights = ett()))
  expect_silent(as.lmitt(lm(y ~ adopters() +
                                       offset(cov_adj(camod, design = des)),
                   data = simdata, weights = ett(des))))

})
