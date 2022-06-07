test_that("basic adopters", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  mod1 <- lm(y ~ z, weights = ate(des), data = simdata)
  mod2 <- lm(y ~ adopters(), weights = ate(des), data = simdata)
  expect_true(all(coef(mod1) == coef(mod2)))

  mod3 <- lm(y ~ z, weights = ate(), offset = cov_adj(mod1, design = des),
             data = simdata)
  mod4 <- lm(y ~ z, weights = ate(des), offset = cov_adj(mod1),
             data = simdata)
  mod5 <- lm(y ~ adopters(), weights = ate(),
             offset = cov_adj(mod1, design = des),
             data = simdata)
  mod6 <- lm(y ~ adopters(), weights = ate(des), offset = cov_adj(mod1),
             data = simdata)
})

test_that("with dichotomy", {

  data(simdata)
  des <- obs_design(dose ~ cluster(cid1, cid2), data = simdata,
                    dichotomy = dose < 250 ~ .)
  m1 <- lm(y ~ adopters(), data = simdata, weights = ate(des))

  expect_true(all(m1$model$`adopters()` %in% 0:1))
})

test_that("Missing data", {

  data(simdata)
  simdata$z[1:4] <- NA

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  expect_equal(length(adopters(des, data = simdata)),
               nrow(simdata))

  simdata$dose[1:4] <- NA
  des <- rct_design(dose ~ uoa(cid1, cid2), data = simdata,
                    dichotomy = dose > 250 ~ .)
  expect_equal(length(adopters(des, data = simdata)),
               nrow(simdata))


  data(STARdata)
  STARdata$id <- seq_len(nrow(STARdata))
  STARdata$stark <- as.character(STARdata$stark)
  des <- rct_design(stark ~ unitid(id), data = STARdata,
                    dichotomy = stark == "small" ~ .)
  mod <- lm(readk ~ birth + lunchk + adopters(), data = STARdata,
          weights = ate(des))

})
