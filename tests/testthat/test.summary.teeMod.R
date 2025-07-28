test_that("summary.teeMod basics", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification = spec, weights = ate())
  summod <- summary(mod)
  expect_s3_class(summod, "summary.teeMod")

  out <- capture_output(print(summod))

  expect_match(out, "Estimate")

  expect_no_match(out, "Residuals")
  expect_no_match(out, "R-squared")

})


test_that("teeMod with SandwichLayer offset summary uses vcov_tee SE's", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)
  cmod <- lm(y ~ x, simdata)
  ssmod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec),
                    offset = cov_adj(cmod)))

  s <- summary(ssmod, vcov.type = "CR0")
  vc <- vcov_tee(ssmod, type = "CR0")
  dof <- length(unique(paste(simdata$uoa1, simdata$uoa2, sep = "_")))-1
  ix <- row.names(vc) != "offset:(Intercept)"
  expect_equal(s$coefficients[, 2L], sqrt(diag(vc)[ix]))
  expect_equal(s$coefficients[, 3L], ssmod$coefficients[ix] / sqrt(diag(vc)[ix]))
  expect_equal(s$coefficients[, 4L],
               2 * pt(abs(ssmod$coefficients[ix] / sqrt(diag(vc)[ix])),
                      dof,
                      lower.tail = FALSE))

  cmod <- lm(y ~ x, simdata, subset = z == 0)
  # Removing intercept which is otherwise estimated at 0, causing warnings in
  # sqrt(diag(vcov)) calculations below
  ssmod <- lmitt(lm(y ~ assigned() + 0, data = simdata, weights = ate(spec),
                    offset = cov_adj(cmod)))

  s <- summary(ssmod, vcov.type = "CR0")
  vc <- vcov_tee(ssmod, type = "CR0")
  ix <- row.names(vc) != "offset:(Intercept)"
  expect_equal(s$coefficients[, 2L], sqrt(diag(vc)[ix]))
  expect_equal(s$coefficients[, 3L], ssmod$coefficients[ix] / sqrt(diag(vc)[ix]))
  expect_equal(s$coefficients[, 4L],
               2 * pt(abs(ssmod$coefficients[ix] / sqrt(diag(vc)[ix])),
                      dof,
                      lower.tail = FALSE))
})

test_that("teeMod w/o SandwichLayer offset summary uses sandwich SE's", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)
  ssmod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec)))

  uoas <- apply(simdata[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))
  s <- summary(ssmod, vcov.type = "CR0", cadjust = FALSE)
  expect_equal(
    s$coefficients[, 2],
    sqrt(diag(sandwich::sandwich(ssmod, meat. = sandwich::meatCL, cluster = uoas,
                                 cadjust = FALSE))))
})

test_that("vcov.type argument", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)
  ssmod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec)))
  expect_equal(summary(ssmod)$vcov.type, "HC0")
  expect_equal(summary(ssmod, vcov.type = "CR0")$vcov.type, "CR0")
  expect_equal(summary(ssmod, vcov.type = "MB0")$vcov.type, "MB0")
  expect_equal(summary(ssmod, vcov.type = "HC0")$vcov.type, "HC0")

  expect_equal(sum(grepl("HC0", capture.output(summary(ssmod)))), 1)
  expect_equal(sum(grepl("CR0",
                         capture.output(summary(ssmod, vcov.type = "CR0")))),
               1)

  expect_equal(sum(grepl("MB0",
                         capture.output(summary(ssmod, vcov.type = "MB0")))),
               1)
  expect_equal(sum(grepl("HC0",
                         capture.output(summary(ssmod, vcov.type = "HC0")))),
               1)

})

test_that("lmitt.form vs as.lmitt", {
  # `lmitt.formula` should report only "Treatment Effects", whereas
  # as.`lmitt`/`lmitt(lm(...` report coefficients since we don't have control
  # over what is reported.
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)

  co <- capture.output(summary(lmitt(y ~ 1, data = simdata, specification = spec)))
  expect_true(any(grepl("Treatment Effects", co)))
  expect_false(any(grepl("Coefficients", co)))

  co <- capture.output(summary(as.lmitt(lm(y ~ z.(spec), data = simdata),
                                        specification = spec)))
  expect_false(any(grepl("Treatment Effects", co)))
  expect_true(any(grepl("Coefficients", co)))

  co <- capture.output(summary(lmitt(lm(y ~ z.(spec), data = simdata),
                                        specification = spec)))
  expect_false(any(grepl("Treatment Effects", co)))
  expect_true(any(grepl("Coefficients", co)))


  ### With interactions, `as.lmitt` should report all coefficients
  mod <- as.lmitt(lm(
    y ~ as.factor(o) + as.factor(o):z.(spec), data = simdata),
    specification = spec)
  suppressWarnings(co <- capture.output(summary(mod), "NaNs"))
  expect_equal(sum(grepl("as.factor(o)", co, fixed = TRUE)),
               8)

})

test_that("issues with coefficients", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)
  mod <- lmitt(y ~ 1, data = simdata, specification =spec, subset = simdata$z == 0)
  s <- summary(mod, vcov.type = "CR0")
  co <- capture.output(s)
  expect_true(any(grepl("Treatment Effects", co)))

  expect_true(any(grepl("calculated via type", co)))

  sslm <- new("teeMod",
              lm(y ~ 0, data = simdata, weights = ate(spec)),
              StudySpecification = spec,
              ctrl_means_model = lm(1~0),
              lmitt_fitted = TRUE)

  expect_true(any(grepl("No Coefficients", capture.output(summary(sslm)))))

})

test_that("catching bug with summary(as.lmitt", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), simdata)

  mod <- lm(y ~ assigned(spec), data = simdata)
  ss <- summary(as.lmitt(mod, spec))
  expect_true(is(ss, "summary.teeMod"))
})

test_that("print.summary isn't confused by bad naming", {
  data(simdata)
  simdata$abz.c <- as.factor(simdata$o)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)

  mod <- lmitt(y ~ abz.c, data = simdata, specification = spec)
  expect_warning(co <- capture.output(summary(mod)), "sufficient degrees of freedom")
  cos <- strsplit(trimws(co), "[[:space:]]+")

  expect_true(!any(grepl("^abz\\.", co)))
  expect_equal(sum(grepl("^z\\.", co)), 4)

  # to force ` in variable names via as.factor
  mod <- lmitt(y ~ as.factor(abz.c), data = simdata, specification = spec)
  expect_warning(co <- capture.output(summary(mod)), "sufficient degrees of freedom")
  expect_true(!any(grepl("^`abz\\.", co)))
  expect_equal(sum(grepl("^`z\\.", co)), 4)

})
