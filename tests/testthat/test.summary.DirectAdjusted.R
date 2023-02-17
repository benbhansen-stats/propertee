test_that("summary.DirectAdjusted basics", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  da <- lmitt(y ~ assigned(), data = simdata, design = des)
  sda <- summary(da)
  expect_s3_class(sda, "summary.DirectAdjusted")

  out <- capture_output(print(sda))

  expect_match(out, "assigned()")

  expect_no_match(out, "Residuals")
  expect_no_match(out, "R-squared")

})


test_that("DirectAdjusted with SandwichLayer offset summary uses vcovDA SE's", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des),
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
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des),
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

test_that("DirectAdjusted w/o SandwichLayer offset summary uses sandwich SE's", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))

  uoas <- apply(simdata[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_"))
  s <- summary(damod, type = "CR0", cadjust = FALSE)
  expect_equal(
    s$coefficients[, 2],
    sqrt(diag(sandwich::sandwich(damod, meat. = sandwich::meatCL, cluster = uoas,
                                 cadjust = FALSE))))
})
