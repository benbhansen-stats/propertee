
test_that("DirectAdjusted with SandwichLayer offset summary uses vcovDA SE's", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  damod <- lmitt(lm(y ~ z, data = simdata, weights = ate(des),
                    offset = cov_adj(cmod)))

  s <- summary(damod)
  expect_equal(s$coefficients[, 2L], sqrt(diag(vcovDA(damod))))
  expect_equal(s$coefficients[, 3L],
               damod$coefficients / sqrt(diag(vcovDA(damod))))
  expect_equal(s$coefficients[, 4L],
               2 * pt(abs(damod$coefficients / sqrt(diag(vcovDA(damod)))),
                      damod$df.residual,
                      lower.tail = FALSE))

  cmod <- lm(y ~ x, simdata, subset = z == 0)
  damod <- lmitt(lm(y ~ z, data = simdata, weights = ate(des),
                    offset = cov_adj(cmod)))
  s <- summary(damod)
  expect_equal(s$coefficients[, 2L], sqrt(diag(vcovDA(damod))))
  expect_equal(s$coefficients[, 3L],
               damod$coefficients / sqrt(diag(vcovDA(damod))))
  expect_equal(s$coefficients[, 4L],
               2 * pt(abs(damod$coefficients / sqrt(diag(vcovDA(damod)))),
                      damod$df.residual,
                      lower.tail = FALSE))
})

test_that("DirectAdjusted w/o SandwichLayer offset summary uses OLS SE's", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  damod <- lmitt(lm(y ~ z, data = simdata, weights = ate(des)))

  s <- summary(damod)
  expect_identical(s, do.call(getS3method("summary", "lm"), list(object = damod)))
})
