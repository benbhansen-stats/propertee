test_that("expand.m.f.da basics", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  mod <- lmitt(y ~ 1, data = simdata, design = des)
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o),
                        flexida:::.expand.model.frame.DA(mod, ~ o)))
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o, na.expand = TRUE),
                        flexida:::.expand.model.frame.DA(mod, ~ o,
                                                         na.expand = TRUE)))

  mod <- lmitt(y ~ 1, data = simdata, design = des, absorb = TRUE)
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o),
                        flexida:::.expand.model.frame.DA(mod, ~ o)))
  expect_true(all.equal(stats::expand.model.frame(mod, ~ o, na.expand = TRUE),
                        flexida:::.expand.model.frame.DA(mod, ~ o,
                                                         na.expand = TRUE)))


})
