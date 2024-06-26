test_that("basic assigned", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ z, weights = ate(des), data = simdata)
  mod2 <- lm(y ~ assigned(), weights = ate(des), data = simdata)
  expect_true(all(coef(mod1) == coef(mod2)))
  expect_true(all(mod2$model$`assigned()` %in% c(0, 1)))
  expect_true(all.equal(mod2$model$`assigned()`, simdata$z))

  mod5 <- lm(y ~ assigned(), weights = ate(),
             offset = cov_adj(mod1, design = des),
             data = simdata)
  mod6 <- lm(y ~ assigned(), weights = ate(des), offset = cov_adj(mod1),
             data = simdata)
  expect_true(all(coef(mod5) == coef(mod6)))
  expect_true(all(mod5$model$`assigned()` %in% c(0, 1)))
  expect_true(all.equal(mod5$model$`assigned()`, simdata$z))
  expect_true(all(mod6$model$`assigned()` %in% c(0, 1)))
  expect_true(all.equal(mod6$model$`assigned()`, simdata$z))
  
  nonbin_des <- rct_design(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod7 <- lm(y ~ assigned(nonbin_des), simdata)
  expect_equal(length(mod7$coefficients), 2)
  expect_true(length(unique(mod7$model$`assigned(nonbin_des)`)) > 2)
  expect_equal(mod7$model$`assigned(nonbin_des)`, simdata$dose)
})

test_that("with dichotomy passed directly", {
  
  data(simdata)
  des <- obs_design(dose ~ cluster(uoa1, uoa2), data = simdata)
  wts <- ate(des, data = simdata, dichotomy = dose < 250 ~ .)
  m1 <- lm(y ~ assigned(des, dichotomy = dose < 250 ~ .), data = simdata,
           weights = wts)
  m2 <- lm(y ~ assigned(des, dichotomy = dose < 250 ~ .), data = simdata)
  
  expect_true(all(m1$model$`assigned()` %in% 0:1))
  expect_true(all(m2$model$`assigned()` %in% 0:1))
})

test_that("with dichotomy found in weights", {

  data(simdata)
  des <- obs_design(dose ~ cluster(uoa1, uoa2), data = simdata)
  m1 <- lm(y ~ assigned(), data = simdata,
           weights = ate(des, dichotomy = dose < 250 ~ .))

  expect_true(all(m1$model$`assigned()` %in% 0:1))
})

test_that("Missing data", {

  data(simdata)
  simdata$z[1:4] <- NA

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  expect_equal(length(assigned(des, data = simdata)),
               nrow(simdata))

  simdata$dose[1:4] <- NA
  des <- rct_design(dose ~ uoa(uoa1, uoa2), data = simdata)
  expect_equal(length(assigned(des, data = simdata, dichotomy = dose > 250 ~ .)),
               nrow(simdata))


  data(STARdata)
  STARdata$id <- seq_len(nrow(STARdata))
  STARdata$stark <- as.character(STARdata$stark)
  des <- rct_design(stark ~ unitid(id), data = STARdata)
  mod <- lm(readk ~ birth + lunchk + assigned(), data = STARdata,
          weights = ate(des, dichotomy = stark == "small" ~ .))

})

test_that("aliases", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lm(y ~ assigned(), weights = ate(des), data = simdata)
  mod2 <- lm(y ~ adopters(), weights = ate(des), data = simdata)
  mod3 <- lm(y ~ a.(), weights = ate(des), data = simdata)
  mod4 <- lm(y ~ z.(), weights = ate(des), data = simdata)
  expect_true(all(coef(mod1) == coef(mod2)))
  expect_true(all(coef(mod1) == coef(mod3)))
  expect_true(all(coef(mod1) == coef(mod4)))
})

test_that("no data argument is equivalent to treatment() #88", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  expect_warning(a1 <- assigned(des), "Extracting treatment")
  a2 <- assigned(des, simdata)
  expect_identical(a1, treatment(des)[, 1])
  expect_true(length(a2) > length(a1))
})
