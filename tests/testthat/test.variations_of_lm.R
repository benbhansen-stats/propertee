# This isn't going to test whether anything is actually correct (see individiual
# tests for that), but instead will run each different variation of `lm` call
# and ensure that no messages, warnings, or errors are produced.

# ******************** IMPORTANT ********************
# If debugging things in here, REMOVE EXPECT_SILENT!!! Running it will suppress
# all output, including inside `browser()`
test_that("binary treatment, in all data", {

  # Treatment exists in data
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  camod <- lm(y ~ x, data = simdata)

  # Weight alone
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec))))
  expect_silent(as.lmitt(lm(y ~ a.(), data = simdata, weights = ate(spec))))
  expect_silent(as.lmitt(lm(y ~ z.(), data = simdata, weights = ate(spec))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata, weights = ate(spec))))

  # Adopters alone
  expect_silent(a <- lm(y ~ assigned(spec), data = simdata))
  # teeMod doesnt look in assigned for StudySpecification
  # expect_silent(as.lmitt(a))

  # cov_adj alone
  expect_silent(as.lmitt(lm(y ~ a.(), data = simdata,
                            offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ z.(spec) + offset(cov_adj(camod, specification = spec)),
                            data = simdata)))


  # Weight + assigned
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(spec), data = simdata,
                            weights = ett(spec))))


  # weight + cov_adj
  expect_silent(as.lmitt(lm(y ~ z.(), data = simdata, weights = ett(spec),
                            offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata, weights = ett(),
                            offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ a.(), data = simdata, weights = ett(spec),
                            offset = cov_adj(camod, specification = spec))))

  # weight + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ z.() + offset(cov_adj(camod)), data = simdata,
                            weights = ett(spec))))
  # Fails when trying to obtain StudySpecification from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ z.() + offset(cov_adj(camod, specification = spec)),
  #                      data = simdata, weights = ett())))
  expect_silent(as.lmitt(lm(y ~ z.() + offset(cov_adj(camod, specification = spec)),
                            data = simdata, weights = ett(spec))))

  # assigned + cov_adj
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(spec), data = simdata,
                            offset = cov_adj(camod, specification = spec))))

  # weights + assigned + cov_adj
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(spec),
                            offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(),
                            offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(spec),
                            offset = cov_adj(camod, specification = spec))))

  # weight + adopter + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ assigned() + offset(cov_adj(camod)),
                            data = simdata, weights = ett(spec))))
  # Fails when trying to obtain StudySpecification from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ assigned() + offset(cov_adj(camod, specification = spec)),
  #                      data = simdata, weights = ett()))
  expect_silent(as.lmitt(lm(y ~ assigned() +
                              offset(cov_adj(camod, specification = spec)),
                            data = simdata, weights = ett(spec))))

})

test_that("binary treatment, not in data2", {
  # Treatment doesn't exist in data
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  camod <- lm(y ~ x, data = simdata)
  simdata$z <- NULL

  # Adopters alone
  expect_silent(a <- lm(y ~ assigned(spec), data = simdata))
  # teeMod doesnt look in assigned for StudySpecification
  # expect_silent(as.lmitt(a))

  # Weight + assigned
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(spec), data = simdata,
                            weights = ett(spec))))

  # assigned + cov_adj
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(spec), data = simdata,
                            offset = cov_adj(camod, specification = spec))))

  # weights + assigned + cov_adj
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                                     weights = ett(spec),
                                     offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                                     weights = ett(),
                                     offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                                     weights = ett(spec),
                                     offset = cov_adj(camod, specification = spec))))

  # weight + adopter + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ assigned() + offset(cov_adj(camod)),
                                     data = simdata, weights = ett(spec))))
  # Fails when trying to obtain StudySpecification from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ assigned() + offset(cov_adj(camod, specification = spec)),
  #                      data = simdata, weights = ett()))
  expect_silent(as.lmitt(lm(y ~ assigned() +
                                       offset(cov_adj(camod, specification = spec)),
                   data = simdata, weights = ett(spec))))

})

test_that("non-binary treatment, in all data, dichotomization in specification", {


  # Treatment exists in data
  data(simdata)

  spec <- rct_spec(dose ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  camod <- lm(y ~ x, data = simdata)


  # Adopters alone
  expect_silent(a <- lm(y ~ assigned(spec, dichotomy = dose >= 200 ~ .), data = simdata))
  # teeMod doesnt look in assigned for StudySpecification
  # expect_silent(as.lmitt(a))

  # Weight + assigned
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                                     weights = ett(spec, dichotomy = dose >= 200 ~ .))))
  expect_silent(as.lmitt(lm(y ~ assigned(spec), data = simdata,
                                     weights = ett(spec, dichotomy = dose >= 200 ~ .))))

  # assigned + cov_adj
  expect_silent(as.lmitt(lm(y ~ assigned(dichotomy = dose >= 200 ~ .), data = simdata,
                            offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(spec, dichotomy = dose >= 200 ~ .), data = simdata,
                            offset = cov_adj(camod, specification = spec))))

  # weights + assigned + cov_adj
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(spec, dichotomy = dose >= 200 ~ .),
                            offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(dichotomy = dose >= 200 ~ .),
                            offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(spec, dichotomy = dose >= 200 ~ .),
                            offset = cov_adj(camod, specification = spec))))

  # weight + adopter + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ assigned() + offset(cov_adj(camod)),
                            data = simdata, weights = ett(spec, dichotomy = dose >= 200 ~ .))))
  # Fails when trying to obtain StudySpecification from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ assigned() + offset(cov_adj(camod, specification = spec)),
  #                      data = simdata, weights = ett()))
  expect_silent(as.lmitt(lm(y ~ assigned() +
                              offset(cov_adj(camod, specification = spec)),
                            data = simdata, weights = ett(spec, dichotomy = dose >= 200 ~ .))))

})

test_that("non-binary treatment, not in data2, dichotomization in specification", {


  # Treatment exists in data
  data(simdata)

  spec <- rct_spec(dose ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  camod <- lm(y ~ x, data = simdata)
  simdata$dose <- NULL


  # Adopters alone
  expect_silent(a <- lm(y ~ assigned(spec, dichotomy = dose >= 200 ~ .), data = simdata))
  # teeMod doesnt look in assigned for StudySpecification
  # expect_silent(as.lmitt(a))

  # Weight + assigned
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(spec, dichotomy = dose >= 200 ~ .))))
  expect_silent(as.lmitt(lm(y ~ assigned(spec), data = simdata,
                            weights = ett(spec, dichotomy = dose >= 200 ~ .))))

  # assigned + cov_adj
  expect_silent(as.lmitt(lm(y ~ assigned(dichotomy = dose >= 200 ~ .), data = simdata,
                   offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ assigned(spec, dichotomy = dose >= 200 ~ .), data = simdata,
                   offset = cov_adj(camod, specification = spec))))

  # weights + assigned + cov_adj
  expect_silent(as.lmitt(lm(y ~ assigned(), data = simdata,
                            weights = ett(spec, dichotomy = dose >= 200 ~ .),
                            offset = cov_adj(camod))))
  expect_silent(as.lmitt(lm(y ~ adopters(), data = simdata,
                            weights = ett(dichotomy = dose >= 200 ~ .),
                            offset = cov_adj(camod, specification = spec))))
  expect_silent(as.lmitt(lm(y ~ z.(), data = simdata,
                            weights = ett(spec, dichotomy = dose >= 200 ~ .),
                            offset = cov_adj(camod, specification = spec))))

  # weight + adopter + cov_adj in formula
  expect_silent(as.lmitt(lm(y ~ assigned() + offset(cov_adj(camod)),
                            data = simdata, weights = ett(spec, dichotomy = dose >= 200 ~ .))))
  # Fails when trying to obtain StudySpecification from a cov_adj inside offset in formula
  #expect_silent(as.lmitt(lm(y ~ assigned() + offset(cov_adj(camod, specification = spec)),
  #                      data = simdata, weights = ett()))
  expect_silent(as.lmitt(lm(y ~ assigned() +
                              offset(cov_adj(camod, specification = spec)),
                            data = simdata, weights = ett(spec, dichotomy = dose >= 200 ~ .))))

})
