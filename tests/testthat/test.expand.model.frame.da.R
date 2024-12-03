test_that("expand.m.f.da basics", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification = spec)
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o),
                        propertee:::.expand.model.frame_teeMod(mod, ~ o)))
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o, na.expand = TRUE),
                        propertee:::.expand.model.frame_teeMod(mod, ~ o,
                                                         na.expand = TRUE)))

  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification = spec, absorb = TRUE)
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o),
                        propertee:::.expand.model.frame_teeMod(mod, ~ o)))
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o, na.expand = TRUE),
                        propertee:::.expand.model.frame_teeMod(mod, ~ o,
                                                         na.expand = TRUE)))


})
