test_that("summary.DirectAdjusted basics", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  da <- lmitt(y ~ 1, data = simdata, design = des, weights = ate())
  sda <- summary(da)
  expect_s3_class(sda, "summary.DirectAdjusted")

  out <- capture_output(print(sda))

  expect_match(out, "Estimate")

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
  # Removing intercept which is otherwise estimated at 0, causing warnings in
  # sqrt(diag(vcov)) calculations below
  damod <- lmitt(lm(y ~ assigned() + 0, data = simdata, weights = ate(des),
                    offset = cov_adj(cmod)))

  s <- summary(damod)
  expect_true(s$coefficients[, 2L] ==  sqrt(diag(vcovDA(damod))))
  expect_true(s$coefficients[, 3L] ==
               damod$coefficients / sqrt(diag(vcovDA(damod))))
  expect_true(s$coefficients[, 4L] ==
               2 * pt(abs(damod$coefficients / sqrt(diag(vcovDA(damod)))),
                      damod$df.residual,
                      lower.tail = FALSE))
})

test_that("DirectAdjusted w/o SandwichLayer offset summary uses sandwich SE's", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))

  uoas <- apply(simdata[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_"))
  s <- summary(damod, vcov.type = "CR0", cadjust = FALSE)
  expect_equal(
    s$coefficients[, 2],
    sqrt(diag(sandwich::sandwich(damod, meat. = sandwich::meatCL, cluster = uoas,
                                 cadjust = FALSE))))
})

test_that("vcov.type argument", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))
  expect_equal(summary(damod)$vcov.type, "CR0")
  expect_equal(summary(damod, vcov.type = "CR0")$vcov.type, "CR0")
  expect_equal(summary(damod, vcov.type = "MB0")$vcov.type, "MB0")
  expect_equal(summary(damod, vcov.type = "HC0")$vcov.type, "HC0")

  expect_equal(sum(grepl("CR0", capture.output(summary(damod)))), 1)
  expect_equal(sum(grepl("CR0",
                         capture.output(summary(damod, vcov.type = "CR0")))),
               1)

  expect_equal(sum(grepl("MB0",
                         capture.output(summary(damod, vcov.type = "MB0")))),
               1)
  expect_equal(sum(grepl("HC0",
                         capture.output(summary(damod, vcov.type = "HC0")))),
               1)

})

test_that("lmitt.form vs as.lmitt", {
  # `lmitt.formula` should report only "Treatment Effects", whereas
  # as.`lmitt`/`lmitt(lm(...` report coefficients since we don't have control
  # over what is reported.
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)

  co <- capture.output(summary(lmitt(y ~ 1, data = simdata, design = des)))
  expect_true(any(grepl("Treatment Effects", co)))
  expect_false(any(grepl("Coefficients", co)))

  co <- capture.output(summary(as.lmitt(lm(y ~ z.(des), data = simdata),
                                        design = des)))
  expect_false(any(grepl("Treatment Effects", co)))
  expect_true(any(grepl("Coefficients", co)))

  co <- capture.output(summary(lmitt(lm(y ~ z.(des), data = simdata),
                                        design = des)))
  expect_false(any(grepl("Treatment Effects", co)))
  expect_true(any(grepl("Coefficients", co)))
})

test_that("issues with coefficients", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)
  mod <- lmitt(y ~ 1, data = simdata, design =des, subset = simdata$z == 0)
  s <- summary(mod)
  co <- capture.output(s)
  expect_true(any(grepl("Treatment Effects", co)))

  expect_false(any(grepl("calculated via type", co)))

  dalm <- new("DirectAdjusted",
              lm(y ~ 0, data = simdata, weights = ate(des)),
              Design = des,
              lmitt_fitted = TRUE)

  expect_true(any(grepl("No Coefficients", capture.output(summary(dalm)))))

})

test_that("catching bug with summary(as.lmitt", {
  data(simdata)
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force), simdata)

  mod <- lm(y ~ assigned(des), data = simdata)
  ss <- summary(as.lmitt(mod, des))
  expect_true(is(ss, "summary.DirectAdjusted"))
})

test_that("#119 flagging as NaN", {
  ### factor moderator variable
  data(simdata)
  copy_simdata <- simdata
  copy_simdata$o_fac <- as.factor(copy_simdata$o)
  des <- rct_design(z ~ cluster(cid1, cid2), copy_simdata)

  ### lmitt.formula
  damod <- lmitt(y ~ o_fac, data = copy_simdata, design = des)
  expect_warning(cf <- coefficients(summary(damod)),
                 "will be returned as NaN: o_fac1, o_fac3")

  # Issue is in subgroup o=1, so *find that* entry in the vcov matrix
  nan_cf <- which(grepl("z._o_fac1", rownames(cf)))
  expect_true(all(is.nan(cf[nan_cf, 2:4])))
  expect_true(all(!is.nan(cf[, 1])))
  expect_true(all(!is.nan(cf[-nan_cf, ])))

  #### lmitt.lm
  damod <- lmitt(lm(y ~ o_fac + o_fac:assigned(des), data = copy_simdata), design = des)
  cf <- coefficients(summary(damod))

  #****************************************
  ### Setting these to NaN manually only for testing purposes!
  cf[1, 2:4] <- NaN
  ### Remove these once #119 is addressed!!!!!
  #****************************************

  # Issue is in subgroup o=1, so the first entry in the vcov matrix
  expect_true(all(is.nan(cf[1, 2:4])))
  expect_true(all(!is.nan(cf[, 1])))
  expect_true(all(!is.nan(cf[-1, ])))

  ### valid continuous moderator variable
  damod <- lmitt(y ~ o, data = copy_simdata, design = des)
  cf <- coefficients(summary(damod))
  expect_true(all(!is.nan(cf)))
  
  ### invalid continuous moderator variable
  copy_simdata$invalid_o <- 0
  copy_simdata$invalid_o[(copy_simdata$cid1 == 2 & copy_simdata$cid2 == 2) |
                      (copy_simdata$cid1 == 2 & copy_simdata$cid2 == 1)] <- 1
  damod <- lmitt(y ~ invalid_o, data = copy_simdata, design = des)
  expect_warning(cf <- coefficients(summary(damod)), "will be returned as NaN: invalid_o")
  nan_cf <- which(grepl("z._invalid_o", cf))
  expect_true(all(is.nan(cf[nan_cf, 2:4])))
  expect_true(all(!is.nan(cf[, nan_cf])))
  expect_true(all(!is.nan(cf[-nan_cf, ])))
})


test_that("print.summary isn't confused by bad naming", {
  data(simdata)
  simdata$abz.c <- as.factor(simdata$o)
  des <- rct_design(z ~ unitid(cid1, cid2), simdata)

  mod <- lmitt(y ~ abz.c, data = simdata, design = des)
  expect_warning(co <- capture.output(summary(mod)), "enough degrees of freedom")
  cos <- strsplit(trimws(co), "[[:space:]]+")

  expect_true(!any(grepl("^abz\\.", co)))
  expect_equal(sum(grepl("^z\\.", co)), 4)

  # to force ` in variable names via as.factor
  mod <- lmitt(y ~ as.factor(abz.c), data = simdata, design = des)
  expect_warning(co <- capture.output(summary(mod)), "enough degrees of freedom")
  expect_true(!any(grepl("^`abz\\.", co)))
  expect_equal(sum(grepl("^`z\\.", co)), 4)

})
