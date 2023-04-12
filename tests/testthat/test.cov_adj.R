old_opt <- options()
on.exit(options(old_opt))
options(flexida_warn_on_conditional_treatment = FALSE)

data(simdata)
data(STARdata)
STARdata$id <- seq_len(nrow(STARdata))
Q_wo_nulls <- simdata
Q_w_nulls <- STARdata
Q_partial_overlap <- simdata[simdata$bid == 3,]

# test whether the values are the same as calling predict on the quasiexperimental
# data if called outside of `lm`
test_ca <- function(ca, cov_mod,  Q_) {
  expect_true(inherits(ca, "numeric"))
  expect_true(inherits(ca, "vector"))
  expect_equal(ca@.Data, as.numeric(stats::predict(cov_mod, Q_)))
  expect_equal(ca@fitted_covariance_model, cov_mod)
  return(NULL)
}

test_that("cov_adj outside of lm call specifying newdata and design, data has NULLs", {
  des <- rct_design(stark == "small" ~ unitid(id), data = Q_w_nulls)
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)

  expect_warning(cov_adj(cmod, design = des, newdata = Q_w_nulls),
                 "adjustments are NA")
  ca <- suppressWarnings(cov_adj(cmod, design = des, newdata = Q_w_nulls))
  test_ca(ca, cmod, Q_w_nulls)
  expect_true(inherits(ca, "SandwichLayer"))
})

test_that("cov_adj outside of lm call specifying newdata and design, data has no NULLs", {
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_wo_nulls)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  ca <- cov_adj(cmod, design = des, newdata = Q_wo_nulls)
  test_ca(ca, cmod, Q_wo_nulls)
  expect_true(inherits(ca, "SandwichLayer"))
})

test_that("cov_adj outside of lm call specifying newdata and design, data has partial overlap", {
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_partial_overlap)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  ca <- cov_adj(cmod, design = des, newdata = Q_partial_overlap)
  test_ca(ca, cmod, Q_partial_overlap)
  expect_true(inherits(ca, "SandwichLayer"))
})

test_that("cov_adj outside of lm call specifying newdata but no design, data has NULLs", {
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  ca <- cov_adj(cmod, newdata = Q_w_nulls)
  test_ca(ca, cmod, Q_w_nulls)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

test_that("cov_adj outside of lm call specifying newdata but no design, data has no NULLs", {
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  ca <- cov_adj(cmod, newdata = Q_wo_nulls)
  test_ca(ca, cmod, Q_wo_nulls)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

test_that("cov_adj outside of lm call specifying newdata but no design, data has partial overlap", {
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  ca <- cov_adj(cmod, newdata = Q_partial_overlap)
  test_ca(ca, cmod, Q_partial_overlap)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

test_that("cov_adj outside of lm call not specifying newdata or design, data has NULLs", {
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  w <- capture_warnings(cov_adj(cmod))
  expect_true(any(vapply(w, grepl, logical(1),
                         pattern = "quasiexperimental data in the call stack")))
  ca <- suppressWarnings(cov_adj(cmod))
  test_ca(ca, cmod, Q_w_nulls)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

test_that("cov_adj outside of lm call not specifying newdata or design, data has no NULLs", {
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  w <- capture_warnings(cov_adj(cmod))
  expect_true(any(vapply(w, grepl, logical(1),
                         pattern = "quasiexperimental data in the call stack")))
  ca <- suppressWarnings(cov_adj(cmod))
  test_ca(ca, cmod, Q_wo_nulls)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

test_that("cov_adj outside of lm call not specifying newdata or design, data has partial overlap", {
  cmod <- lm(y ~ x, data = Q_partial_overlap)
  w <- capture_warnings(cov_adj(cmod))
  expect_true(any(vapply(w, grepl, logical(1),
                         pattern = "quasiexperimental data in the call stack")))
  ca <- suppressWarnings(cov_adj(cmod))
  test_ca(ca, cmod, Q_partial_overlap)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

test_that("cov_adj as offset with weights, data has NULLs", {
  des <- rct_design(stark == "small" ~ unitid(id), data = Q_w_nulls)
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  expect_warning(lm(readk ~ stark == "small", data = Q_w_nulls,
                    offset = cov_adj(cmod), weights = ate(des)),
                 "adjustments are NA")
  m <- suppressWarnings(lm(readk ~ stark == "small", data = Q_w_nulls,
                           offset = cov_adj(cmod), weights = ate(des)))

  form <- paste("~", paste(unique(c(all.vars(formula(cmod)), all.vars(formula(m)))),
                           collapse = "+"))
  test_ca(m$model$`(offset)`,
          cmod,
          stats::model.frame(as.formula(form), data = Q_w_nulls))
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

test_that("cov_adj as offset with weights, data has no NULLs", {
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_wo_nulls)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  m <- lm(y ~ z, data = Q_wo_nulls, offset = cov_adj(cmod), weights = ate(des))
  test_ca(m$model$`(offset)`, cmod, Q_wo_nulls)
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

test_that("cov_adj as offset with weights, data has partial overlap", {
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_partial_overlap)
  cmod <- lm(y ~ x, data = Q_partial_overlap)
  m <- lm(y ~ z, data = Q_partial_overlap, offset = cov_adj(cmod), weights = ate(des))
  test_ca(m$model$`(offset)`, cmod, Q_partial_overlap)
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

test_that("cov_adj as offset specified w/ newdata and design, no weights, data has NULLs", {
  des <- rct_design(stark == "small" ~ unitid(id), data = Q_w_nulls)
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  expect_warning(lm(readk ~ stark == "small", data = Q_w_nulls,
                    offset = cov_adj(cmod, newdata = Q_w_nulls, design = des)),
                 "adjustments are NA")
  m <- suppressWarnings(
    lm(readk ~ stark == "small", data = Q_w_nulls,
       offset = cov_adj(cmod, newdata = Q_w_nulls, design = des)))
  all_vars <- c(all.vars(cmod$call$formula[-2]), all.vars(m$call$formula))
  keep_idx <- apply(is.na(Q_w_nulls[, all_vars]), 1, sum) == 0
  test_ca(m$model$`(offset)`, cmod, Q_w_nulls[keep_idx, ])
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

test_that("cov_adj as offset specified w/ newdata and design, no weights, data has no NULLs", {
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_wo_nulls)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  m <- lm(y ~ z, data = Q_wo_nulls,
          offset = cov_adj(cmod, newdata = Q_wo_nulls, design = des))
  test_ca(m$model$`(offset)`, cmod, Q_wo_nulls)
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

test_that("cov_adj as offset specified w/ newdata and design, no weights, data has partial overlap", {
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_partial_overlap)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  m <- lm(y ~ z, data = Q_partial_overlap,
          offset = cov_adj(cmod, newdata = Q_partial_overlap, design = des))
  test_ca(m$model$`(offset)`, cmod, Q_partial_overlap)
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

test_that("cov_adj as offset specified w/ no newdata nor design, no weights, data has NULLs", {
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  m <- lm(readk ~ stark == "small", data = Q_w_nulls, offset = cov_adj(cmod))
  all_vars <- c(all.vars(cmod$call$formula[-2]), all.vars(m$call$formula))
  keep_idx <- apply(is.na(Q_w_nulls[, all_vars]), 1, sum) == 0
  test_ca(m$model$`(offset)`, cmod, Q_w_nulls[keep_idx, ])
  expect_true(inherits(m$model$`(offset)`, "PreSandwichLayer"))
})

test_that("cov_adj as offset specified w/ no newdata nor design, no weights, data has no NULLs", {
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  m <- lm(y ~ z, data = Q_wo_nulls, offset = cov_adj(cmod))
  test_ca(m$model$`(offset)`, cmod, Q_wo_nulls)
  expect_true(inherits(m$model$`(offset)`, "PreSandwichLayer"))
})

test_that("cov_adj as offset specified w/ no newdata nor design, no weights, data has partial overlap", {
  cmod <- lm(y ~ x, data = Q_partial_overlap)
  m <- lm(y ~ z, data = Q_partial_overlap, offset = cov_adj(cmod))
  test_ca(m$model$`(offset)`, cmod, Q_partial_overlap)
  expect_true(inherits(m$model$`(offset)`, "PreSandwichLayer"))
})

test_that("cov_adj variance estimates for orthogonal predictors", {
  library(sandwich)
  set.seed(230274373)
  k <- 25
  n <- 4 * k
  x1 <- c(rep(1, 2 * k), rep(-1, 2 * k))
  x2 <- c(rep(1, k), rep(-1, 2 * k), rep(1, k))
  z <- c(rep(1, k), rep(-1, k), rep(1, k), rep(-1, k))

  # all orthogonal by design
  expect_equal((t(x1) %*% x2)[1], 0)
  expect_equal((t(x1) %*% z)[1], 0)
  expect_equal((t(x2) %*% z)[1], 0)

  xy <- function(b0 =  1, b1 = 2, b2 = 3, b3 = -1, sigma2 = 4) {
    y <- b0 + b1 * x1 + b2 * x2 + b3 * z + rnorm(n, sd = sqrt(sigma2))
    data.frame(id = 1:n, x1 = x1, x2 = x2, z = z, y = y)
  }

  df <- xy()
  mboth <- glm(y ~ x1 + x2 + z, data = df)
  m1 <- glm(y ~ x1, data = df)

  # we'll start by manually doing the offest Y - \hat Y
  m2 <- glm(y ~ x2 + z, data = df, offset = predict(m1))

  # Because of orthogonality, we should get $\hat \beta_0= \hat \alpha_0$,
  # $\hat \beta_1 = \hat \alpha_1$, and $\hat beta_2 = \hat gamma_1$.

  expect_equal(coef(mboth), c(coef(m1), coef(m2)[2:3]))

  # now repeat with cov_adj
  design <- rct_design(z == 1 ~ unitid(id), data = df)
  m2ca <- glm(y ~ x2 + z, data = df, offset = cov_adj(m1, design = design))

  expect_equal(coef(mboth), c(coef(m1), coef(m2ca)[2:3]))

  ## naive case
  m2naive <- glm(y ~ x2 + z, data = df)

  expect_false(all(vcov(m2naive) == vcov(m2ca)))

  hc_both  <- vcovHC(mboth, type = "HC0")
  hc_naive <- vcovHC(m2naive, type = "HC0")
  hc_m2ca  <- vcovHC(m2ca, type = "HC0")

  # trim down the bigger model to just the variables in the second stage
  in_two <- c("(Intercept)", "x2", "z")
  hc_both_trimmed <- hc_both[in_two, in_two]

  expect_equal(hc_m2ca, hc_both_trimmed)
  expect_false(all(hc_naive == hc_m2ca))


})


test_that("cov_adj variance estimates for correlated predictors", {
  library(sandwich)
  set.seed(230274373)
  k <- 25
  n <- 4 * k
  x1 <- c(rep(1, 2 * k), rep(-1, 2 * k))
  z <- sample(c(rep(1, 2 * k), rep(0, 2 * k)), prob = x1 + 2)

  xy <- function(b0 =  1, b1 = 2, b2 = 3, sigma2 = 4) {
    y <- b0 + b1 * x1 + b2 * z + rnorm(n, sd = sqrt(sigma2))
    data.frame(id = 1:n, x1 = x1,  z = z, y = y)
  }

  df <- xy()
  mboth <- lm(y ~ x1 + z, data = df)
  m1 <- lm(y ~ x1, data = df) # includes an intercept

  # we'll start by manually doing the offset Y - \hat Y
  # the intercept is in the "x" variables, not the z variables
  m2 <- lm(y ~ z - 1, data = df, offset = predict(m1))

  ## verify some of the results of the VarianceEstimation vignette
  X <- model.matrix(m1)
  Z <- model.matrix(m2)
  Id <- diag(n)

  # some quantities we will need
  XtXinv <- solve(t(X) %*% X)
  ZtZinv <- solve(t(Z) %*% Z)
  H <- Id - X %*% XtXinv %*% t(X)
  G <- Id - Z %*% ZtZinv %*% t(Z)

  a <- solve(t(X) %*% G %*% X)
  b <- solve(t(Z) %*% H %*% Z)
  s2_1 <- (1/(n - 3)) * t(df$y) %*% (Id - X %*% a %*% t(X) %*% G - Z %*% b %*% t(Z) %*% H) %*% df$y
  s2_1 <- s2_1[1,1]
  expect_equal(summary(mboth)$sigma, sqrt(s2_1))

  # variance we will get from the combined model (just subsetting for the 'z' variable)
  var_beta_1 <- s2_1 * b
  expect_equal(vcov(mboth)[3,3], var_beta_1[1,1])

  alpha_2 <- (XtXinv %*% t(X) %*% df$y)[,1]
  beta_2 <- (ZtZinv %*% t(Z) %*% H %*% df$y)[1,1]
  expect_equal(coef(m1), alpha_2)
  expect_equal(coef(m2), beta_2)

  s2_2 <- (1 / (n - 1)) * t(df$y) %*% H %*% G %*% H %*% df$y
  s2_2 <- s2_2[1,1]
  expect_equal(summary(m2)$sigma, sqrt(s2_2))

  r_2 <- resid(m2)
  expect_equal(summary(m2)$sigma, sqrt((1/(n - 1)) * sum(r_2^2)))
  expect_equal(r_2, df$y - X %*% alpha_2 - Z %*% beta_2, ignore_attr = TRUE)
  expect_equal(r_2, G %*% H %*% df$y, ignore_attr = TRUE)

  # now repeat with cov_adj
  design <- rct_design(z ~ unitid(id), data = df)
  m2ca <- lm(y ~ assigned() - 1, data = df, offset = cov_adj(m1, design = design), weights = ate(design))

  m2ca_da <- as.lmitt(m2ca)

  ## TODO: What are we guaranteeing about vcov and sandwich about m2ca_da?
})

options(old_opt)

# This currently errors with
# Error in model.frame.default(formula = readk ~ birth + lunchk, data =
#STARdata, : variable lengths differ (found for '(weights)')
# mod <- lm(readk ~ birth + lunchk, data = STARdata,
#          offset = cov_adj(cmod, newdata = STARdata), weights = ate(des))



test_that("Basics of replacing treatment variable with reference level", {
  data(simdata)

  # Binary treatment
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  camod <- lm(y ~ x + z, data = simdata)
  ca <- cov_adj(camod, newdata = simdata, design = des)

  simdata2 <- simdata
  simdata2$z <- 0

  manual <- predict(camod, newdata = simdata2)
  expect_true(all(manual == ca))

  ### Let's just make sure we're not getting spurious positive results
  simdata2 <- simdata
  manual <- predict(camod, newdata = simdata2)
  expect_false(all(manual == ca))

  # Numeric treatment
  des <- rct_design(dose ~ cluster(cid1, cid2), data = simdata)
  camod <- lm(y ~ x + dose, data = simdata)
  ca <- cov_adj(camod, newdata = simdata, design = des)

  simdata2 <- simdata
  simdata2$dose <- 50

  manual <- predict(camod, newdata = simdata2)
  expect_true(all(manual == ca))

  # Factor treatment
  ## simdata$dose <- as.factor(simdata$dose)
  ## des <- rct_design(dose ~ cluster(cid1, cid2), data = simdata)
  ## camod <- lm(y ~ x + dose, data = simdata)
  ## ca <- cov_adj(camod, newdata = simdata, design = des)

  ## simdata2 <- simdata
  ## simdata2$dose <- levels(simdata2$dose)[1]

  ## manual <- predict(camod, newdata = simdata2)
  ## expect_true(all(manual == ca))
  # Current build does NOT allow factor treatment so this will fail

})
