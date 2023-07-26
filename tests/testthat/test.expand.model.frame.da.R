test_that("expand.m.f.da basics", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des)
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o),
                        propertee:::.expand.model.frame.DA(mod, ~ o)))
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o, na.expand = TRUE),
                        propertee:::.expand.model.frame.DA(mod, ~ o,
                                                         na.expand = TRUE)))

  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE)
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o),
                        propertee:::.expand.model.frame.DA(mod, ~ o)))
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o, na.expand = TRUE),
                        propertee:::.expand.model.frame.DA(mod, ~ o,
                                                         na.expand = TRUE)))


})
