test_that("basic assigned", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, weights = ate(spec), data = simdata)
  mod2 <- lm(y ~ assigned(), weights = ate(spec), data = simdata)
  expect_true(all(coef(mod1) == coef(mod2)))
  expect_true(all(mod2$model$`assigned()` %in% c(0, 1)))
  expect_true(all.equal(mod2$model$`assigned()`, simdata$z))
  
  mod5 <- lm(y ~ assigned(), weights = ate(),
             offset = cov_adj(mod1, specification = spec),
             data = simdata)
  mod6 <- lm(y ~ assigned(), weights = ate(spec), offset = cov_adj(mod1),
             data = simdata)
  expect_true(all(coef(mod5) == coef(mod6)))
  expect_true(all(mod5$model$`assigned()` %in% c(0, 1)))
  expect_true(all.equal(mod5$model$`assigned()`, simdata$z))
  expect_true(all(mod6$model$`assigned()` %in% c(0, 1)))
  expect_true(all.equal(mod6$model$`assigned()`, simdata$z))
  
  nonbin_spec <- rct_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod7 <- lm(y ~ assigned(nonbin_spec), simdata)
  expect_equal(length(mod7$coefficients), 2)
  expect_true(length(unique(mod7$model$`assigned(nonbin_spec)`)) > 2)
  expect_equal(mod7$model$`assigned(nonbin_spec)`, simdata$dose)
})

test_that("with dichotomy passed directly", {
  
  data(simdata)
  spec <- obs_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  wts <- ate(spec, data = simdata, dichotomy = dose < 250 ~ .)
  m1 <- lm(y ~ assigned(spec, dichotomy = dose < 250 ~ .), data = simdata,
           weights = wts)
  m2 <- lm(y ~ assigned(spec, dichotomy = dose < 250 ~ .), data = simdata)
  
  expect_true(all(m1$model$`assigned()` %in% 0:1))
  expect_true(all(m2$model$`assigned()` %in% 0:1))
})

test_that("with dichotomy found in weights", {
  
  data(simdata)
  spec <- obs_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  m1 <- lm(y ~ assigned(), data = simdata,
           weights = ate(spec, dichotomy = dose < 250 ~ .))
  
  expect_true(all(m1$model$`assigned()` %in% 0:1))
})

test_that("Missing data", {
  
  data(simdata)
  simdata$z[1:4] <- NA
  
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  expect_equal(length(assigned(spec, data = simdata)),
               nrow(simdata))
  
  simdata$dose[1:4] <- NA
  spec <- rct_spec(dose ~ uoa(uoa1, uoa2), data = simdata)
  expect_equal(length(assigned(spec, data = simdata, dichotomy = dose > 250 ~ .)),
               nrow(simdata))
  
  
  if (requireNamespace("AER", quietly = TRUE)) {
    data(STARdata)
    STARdata <- STAR
    STARdata$id <- seq_len(nrow(STARdata))
    STARdata$stark <- as.character(STARdata$stark)
    spec <- rct_spec(stark ~ unitid(id), data = STARdata)
    mod <- lm(readk ~ birth + lunchk + assigned(), data = STARdata,
              weights = ate(spec, dichotomy = stark == "small" ~ .))
  }
  
})

test_that("aliases", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ assigned(), weights = ate(spec), data = simdata)
  mod2 <- lm(y ~ adopters(), weights = ate(spec), data = simdata)
  mod3 <- lm(y ~ a.(), weights = ate(spec), data = simdata)
  mod4 <- lm(y ~ z.(), weights = ate(spec), data = simdata)
  expect_true(all(coef(mod1) == coef(mod2)))
  expect_true(all(coef(mod1) == coef(mod3)))
  expect_true(all(coef(mod1) == coef(mod4)))
})

test_that("no data argument is equivalent to treatment() #88", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  
  expect_warning(a1 <- assigned(spec), "Extracting treatment")
  a2 <- assigned(spec, simdata)
  expect_identical(a1, treatment(spec)[, 1])
  expect_true(length(a2) > length(a1))
})
