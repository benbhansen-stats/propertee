old_opt <- options()
on.exit(options(old_opt))
options(propertee_warn_on_conditional_treatment = FALSE)

data(simdata)
  if (requireNamespace("AER", quietly = TRUE)) {
      data("STAR", package="AER")
      STARdata <- STAR
      STARdata$id <- seq_len(nrow(STARdata))
      Q_w_nulls <- STARdata
      }
Q_wo_nulls <- simdata
Q_partial_overlap <- simdata[simdata$bid == 3,]

# test whether the values are the same as calling predict on the direct adjustment
# data if called outside of `lm`
test_ca <- function(ca, cov_mod,  Q_) {
  expect_true(inherits(ca, "numeric"))
  expect_true(inherits(ca, "vector"))
  expect_equal(ca@.Data, as.numeric(stats::predict(cov_mod, Q_)))
  expect_equal(ca@fitted_covariance_model, cov_mod)
  return(NULL)
}

  if (requireNamespace("AER", quietly = TRUE)) {
test_that("cov_adj outside of lm call specifying newdata and specification, data has NULLs", {
  spec <- rct_spec(stark == "small" ~ unitid(id), data = Q_w_nulls)
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)

  expect_warning(cov_adj(cmod, specification = spec, newdata = Q_w_nulls),
                 "adjustments are NA")
  ca <- suppressWarnings(cov_adj(cmod, specification = spec, newdata = Q_w_nulls))
  test_ca(ca, cmod, Q_w_nulls)
  expect_true(inherits(ca, "SandwichLayer"))
})
}

test_that("cov_adj outside of lm call specifying newdata and specification, data has no NULLs", {
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_wo_nulls)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  ca <- cov_adj(cmod, specification = spec, newdata = Q_wo_nulls)
  test_ca(ca, cmod, Q_wo_nulls)
  expect_true(inherits(ca, "SandwichLayer"))
})

test_that("cov_adj outside of lm call specifying newdata and specification, data has partial overlap", {
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_partial_overlap)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  ca <- cov_adj(cmod, specification = spec, newdata = Q_partial_overlap)
  test_ca(ca, cmod, Q_partial_overlap)
  expect_true(inherits(ca, "SandwichLayer"))
})

  if (requireNamespace("AER", quietly = TRUE)) {
test_that("cov_adj outside of lm call specifying newdata but no specification, data has NULLs", {
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  ca <- cov_adj(cmod, newdata = Q_w_nulls)
  test_ca(ca, cmod, Q_w_nulls)
  expect_true(inherits(ca, "PreSandwichLayer"))
})
  }

test_that("cov_adj outside of lm call specifying newdata but no specification, data has no NULLs", {
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  ca <- cov_adj(cmod, newdata = Q_wo_nulls)
  test_ca(ca, cmod, Q_wo_nulls)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

test_that("cov_adj outside of lm call specifying newdata but no specification, data has partial overlap", {
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  ca <- cov_adj(cmod, newdata = Q_partial_overlap)
  test_ca(ca, cmod, Q_partial_overlap)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

  if (requireNamespace("AER", quietly = TRUE)) {
test_that("cov_adj outside of lm call not specifying newdata or specification, data has NULLs", {
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  w <- capture_warnings(cov_adj(cmod))
  expect_true(any(vapply(w, grepl, logical(1),
                         pattern = "direct adjustment data in the call stack")))
  ca <- suppressWarnings(cov_adj(cmod))
  test_ca(ca, cmod, Q_w_nulls)
  expect_true(inherits(ca, "PreSandwichLayer"))
})
}
test_that("cov_adj outside of lm call not specifying newdata or specification, data has no NULLs", {
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  w <- capture_warnings(cov_adj(cmod))
  expect_true(any(vapply(w, grepl, logical(1),
                         pattern = "direct adjustment data in the call stack")))
  ca <- suppressWarnings(cov_adj(cmod))
  test_ca(ca, cmod, Q_wo_nulls)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

test_that("cov_adj outside of lm call not specifying newdata or specification, data has partial overlap", {
  cmod <- lm(y ~ x, data = Q_partial_overlap)
  w <- capture_warnings(cov_adj(cmod))
  expect_true(any(vapply(w, grepl, logical(1),
                         pattern = "direct adjustment data in the call stack")))
  ca <- suppressWarnings(cov_adj(cmod))
  test_ca(ca, cmod, Q_partial_overlap)
  expect_true(inherits(ca, "PreSandwichLayer"))
})

  if (requireNamespace("AER", quietly = TRUE)) {
test_that("cov_adj as offset with weights, data has NULLs", {
  spec <- rct_spec(stark == "small" ~ unitid(id), data = Q_w_nulls)
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  expect_warning(lm(readk ~ stark == "small", data = Q_w_nulls,
                    offset = cov_adj(cmod), weights = ate(spec)),
                 "adjustments are NA")
  m <- suppressWarnings(lm(readk ~ stark == "small", data = Q_w_nulls,
                           offset = cov_adj(cmod), weights = ate(spec)))

  form <- paste("~", paste(unique(c(all.vars(formula(cmod)), all.vars(formula(m)))),
                           collapse = "+"))
  test_ca(m$model$`(offset)`,
          cmod,
          stats::model.frame(as.formula(form), data = Q_w_nulls))
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})
  }

test_that("cov_adj as offset with weights, data has no NULLs", {
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_wo_nulls)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  m <- lm(y ~ z, data = Q_wo_nulls, offset = cov_adj(cmod), weights = ate(spec))
  test_ca(m$model$`(offset)`, cmod, Q_wo_nulls)
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

test_that("cov_adj as offset with weights, data has partial overlap", {
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_partial_overlap)
  cmod <- lm(y ~ x, data = Q_partial_overlap)
  m <- lm(y ~ z, data = Q_partial_overlap, offset = cov_adj(cmod), weights = ate(spec))
  test_ca(m$model$`(offset)`, cmod, Q_partial_overlap)
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

  if (requireNamespace("AER", quietly = TRUE)) {
test_that("cov_adj as offset specified w/ newdata and specification, no weights, data has NULLs", {
  spec <- rct_spec(stark == "small" ~ unitid(id), data = Q_w_nulls)
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  expect_warning(lm(readk ~ stark == "small", data = Q_w_nulls,
                    offset = cov_adj(cmod, newdata = Q_w_nulls, specification = spec)),
                 "adjustments are NA")
  m <- suppressWarnings(
    lm(readk ~ stark == "small", data = Q_w_nulls,
       offset = cov_adj(cmod, newdata = Q_w_nulls, specification = spec)))
  all_vars <- c(all.vars(cmod$call$formula[-2]), all.vars(m$call$formula))
  keep_idx <- apply(is.na(Q_w_nulls[, all_vars]), 1, sum) == 0
  test_ca(m$model$`(offset)`, cmod, Q_w_nulls[keep_idx, ])
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})
  }

test_that("cov_adj as offset specified w/ newdata and specification, no weights, data has no NULLs", {
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_wo_nulls)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  m <- lm(y ~ z, data = Q_wo_nulls,
          offset = cov_adj(cmod, newdata = Q_wo_nulls, specification = spec))
  test_ca(m$model$`(offset)`, cmod, Q_wo_nulls)
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

test_that("cov_adj as offset specified w/ newdata and specification, no weights, data has partial overlap", {
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_partial_overlap)
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  m <- lm(y ~ z, data = Q_partial_overlap,
          offset = cov_adj(cmod, newdata = Q_partial_overlap, specification = spec))
  test_ca(m$model$`(offset)`, cmod, Q_partial_overlap)
  expect_true(inherits(m$model$`(offset)`, "SandwichLayer"))
})

  if (requireNamespace("AER", quietly = TRUE)) {
test_that("cov_adj as offset specified w/ no newdata nor specification, no weights, data has NULLs", {
  cmod <- lm(readk ~ gender + ethnicity, data = Q_w_nulls)
  m <- lm(readk ~ stark == "small", data = Q_w_nulls, offset = cov_adj(cmod))
  all_vars <- c(all.vars(cmod$call$formula[-2]), all.vars(m$call$formula))
  keep_idx <- apply(is.na(Q_w_nulls[, all_vars]), 1, sum) == 0
  test_ca(m$model$`(offset)`, cmod, Q_w_nulls[keep_idx, ])
  expect_true(inherits(m$model$`(offset)`, "PreSandwichLayer"))
})
  }

test_that("cov_adj as offset specified w/ no newdata nor specification, no weights, data has no NULLs", {
  cmod <- lm(y ~ x, data = Q_wo_nulls)
  m <- lm(y ~ z, data = Q_wo_nulls, offset = cov_adj(cmod))
  test_ca(m$model$`(offset)`, cmod, Q_wo_nulls)
  expect_true(inherits(m$model$`(offset)`, "PreSandwichLayer"))
})

test_that("cov_adj as offset specified w/ no newdata nor specification, no weights, data has partial overlap", {
  cmod <- lm(y ~ x, data = Q_partial_overlap)
  m <- lm(y ~ z, data = Q_partial_overlap, offset = cov_adj(cmod))
  test_ca(m$model$`(offset)`, cmod, Q_partial_overlap)
  expect_true(inherits(m$model$`(offset)`, "PreSandwichLayer"))
})

test_that("cov_adj with named `by` argument called within `lmitt`", {
  on.exit(Q_wo_nulls <- Q_wo_nulls[, !(colnames(Q_wo_nulls) == "obs_id")])
  Q_wo_nulls$obs_id <- seq_len(nrow(Q_wo_nulls))

  cmod <- lm(y ~ x, Q_wo_nulls)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), Q_wo_nulls)
  m <- lmitt(y ~ 1, specification = spec, data = Q_wo_nulls,
             offset = cov_adj(cmod, by = "obs_id"))
  test_ca(m$model$`(offset)`, cmod, Q_wo_nulls)
})

test_that("cov_adj with `by` argument that's the same as the default called within `lmitt`", {
  cmod <- lm(y ~ x, Q_wo_nulls)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), Q_wo_nulls)
  m <- lmitt(y ~ 1, specification = spec, data = Q_wo_nulls,
             offset = cov_adj(cmod, by = c("uoa1", "uoa2")))
  test_ca(m$model$`(offset)`, cmod, Q_wo_nulls)
})

test_that(paste("cov_adj without `lmitt` call with named `by` argument and cmod",
                "fit without `data` arg"), {
  on.exit(Q_wo_nulls <- Q_wo_nulls[, !(colnames(Q_wo_nulls) == "obs_id")])
  C_data <- Q_wo_nulls
  Q_wo_nulls$obs_id <- seq_len(nrow(Q_wo_nulls))
  C_data$uid <- seq_len(nrow(C_data))

  cmod <- lm(C_data$y ~ C_data$x)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), Q_wo_nulls)

  expect_error(
    expect_warning(
      expect_warning(
        expect_warning(
          cov_adj(cmod, specification = spec, by = c("obs_id" = "uid")),
          "No call to"
        ),
        "Unable to detect"
      ),
      "Could not find direct adjustment data"
    ),
    "`model` must be fit"
  )
})

test_that(paste("cov_adj with unnamed `by` argument without `lmitt` call and cmod",
                "fit with `data` arg"), {
  C_data <- Q_wo_nulls[, !(colnames(Q_wo_nulls) %in% c("uoa1", "uoa2"))]
  C_data$c1 <- Q_wo_nulls$uoa1
  C_data$c2 <- Q_wo_nulls$uoa2

  cmod <- lm(y ~ x, Q_wo_nulls)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), Q_wo_nulls)

  expect_warning(
    expect_warning(
      expect_warning(
        ca <- cov_adj(cmod, specification = spec, by = c("uoa1", "uoa2")),
        "No call to"
      ),
      "Unable to detect"
    ),
    "Could not find direct adjustment data"
  )

  test_ca(ca, cmod, Q_wo_nulls)
})

test_that("cov_adj coltypes for treatment", {
  set.seed(52)
  n <- 10
  cmod_data <- data.frame(x = rnorm(n), y = rnorm(n), uid = 10 + seq_len(10))
  cmod <- lm(y ~ x, cmod_data)
  specdata <- data.frame(x = rnorm(n), a = factor(rep(c(0, 1), each = 5)),
                        uid = seq_len(10))
  spec <- rct_spec(a ~ unitid(uid), specdata)
  ca <- cov_adj(cmod, specdata, spec)
  expect_true(inherits(ca, "SandwichLayer"))

  specdata$a <- as.character(specdata$a)
  spec <- rct_spec(a ~ unitid(uid), specdata)
  expect_warning(cov_adj(cmod, specdata, spec))
})

test_that(paste(".update_ca_model_formula with NULL `by` and NULL `specification` returns",
                "returns original model formula"), {
  set.seed(819)
  df <- data.frame(x = rnorm(10), y = rnorm(10), uid = seq_len(10))

  model <- lm(y ~ x, df)
  expect_equal(deparse1(.update_ca_model_formula(model)), "y ~ x")
})

test_that(paste(".update_ca_model_formula with named `by`"), {
  set.seed(819)
  df <- data.frame(x = rnorm(10), y = rnorm(10), uid = seq_len(10))

  model <- lm(y ~ x, df)
  expect_equal(deparse1(.update_ca_model_formula(model, by = c("uoa1" = "uid"))),
               "y ~ x + uoa1")
})

test_that(paste(".update_ca_model_formula with unnamed `by`"), {
  set.seed(819)
  df <- data.frame(x = rnorm(10), y = rnorm(10), uid = seq_len(10))

  model <- lm(y ~ x, df)
  expect_equal(deparse1(.update_ca_model_formula(model, by = "uid")),
               "y ~ x + uid")
})

test_that(paste(".update_ca_model_formula with NULL `by` and non-NULL `specification`"), {
  set.seed(819)
  df <- data.frame(x = rnorm(10), y = rnorm(10), uid = seq_len(10))
  specification_df <- data.frame(z = rep(c(0, 1), each = 5), uid = seq_len(10))
  spec <- rct_spec(z ~ unitid(uid), specification_df)

  model <- lm(y ~ x, df)
  expect_equal(deparse1(.update_ca_model_formula(model, specification = spec)),
               "y ~ x + uid")
})

test_that(paste(".update_ca_model_formula with non-NULL `by` and non-NULL `specification`"), {
  set.seed(819)
  df <- data.frame(x = rnorm(10), y = rnorm(10), uid = seq_len(10))
  specification_df <- data.frame(z = rep(c(0, 1), each = 5), clust = rep(c(1, 2), each = 5), uid = seq_len(10))
  spec <- rct_spec(z ~ cluster(clust), specification_df)

  model <- lm(y ~ x, df)
  expect_equal(deparse1(.update_ca_model_formula(model, by = "uid", specification = spec)),
               "y ~ x + clust + uid")
})

test_that("cov_adj variance estimates for orthogonal predictors", {
  library(sandwich)
  set.seed(230274373)
  k <- 25
  n <- 4 * k
  x1 <- c(rep(1, 2 * k), rep(-1, 2 * k))
  x2 <- c(rep(1, k), rep(-1, 2 * k), rep(1, k))
  z <- c(rep(1, k), rep(-1, k), rep(1, k), rep(-1, k))

  # all orthogonal by specification
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
  specification <- rct_spec(z == 1 ~ unitid(id), data = df)
  m2ca <- glm(y ~ x2 + z, data = df, offset = cov_adj(m1, specification = specification))

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

## We'll go on to check agreement in cases with non-orthogonal
## predictors below, after verifying that when a treatment variable
## appears in the covariance model its values are replaced with a
## reference level before generating the predictions/offsets for the
## computation of diff-in-diffs.

options(old_opt)




test_that("Basics of replacing treatment variable with reference level", {
  data(simdata)

  # Binary treatment
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  camod <- lm(y ~ x + z, data = simdata)
  ca <- cov_adj(camod, newdata = simdata, specification = spec)

  simdata2 <- simdata
  simdata2$z <- 0

  manual <- predict(camod, newdata = simdata2)
  expect_true(all(manual == ca))

  ### Let's just make sure we're not getting spurious positive results
  simdata2 <- simdata
  manual <- predict(camod, newdata = simdata2)
  expect_false(all(manual == ca))

  # Numeric treatment
  spec <- rct_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  camod <- lm(y ~ x + dose, data = simdata)
  ca <- cov_adj(camod, newdata = simdata, specification = spec)

  simdata2 <- simdata
  simdata2$dose <- 50

  manual <- predict(camod, newdata = simdata2)
  expect_true(all(manual == ca))

  # Factor treatment
  ## simdata$dose <- as.factor(simdata$dose)
  ## spec <- rct_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  ## camod <- lm(y ~ x + dose, data = simdata)
  ## ca <- cov_adj(camod, newdata = simdata, specification = spec)

  ## simdata2 <- simdata
  ## simdata2$dose <- levels(simdata2$dose)[1]

  ## manual <- predict(camod, newdata = simdata2)
  ## expect_true(all(manual == ca))
  # Current build does NOT allow factor treatment so this will fail

})

test_that("cov_adj sets treatment =0 when there is an interaction in cmod", {
  n <- 1000
  Rslope <- 0.5
  Xslope <- 0.2
  R <- sample(seq(-1,1,length.out=30),n,replace=TRUE)
  Z <- R>0
  x <- rnorm(n)
  P <- plogis(Rslope*R+Xslope*x)

  Y <- rbinom(1000,1,P)
  dat <- data.frame(R=R,Z=Z,x=x,Y=Y)

  dat$id=1:nrow(dat)
  spec=rd_spec(Z~forcing(R)+unitid(id),data=dat)

  modHet=glm(Y~R+Z*x,family=binomial,data=dat)

  ## propertee cov_adj
  ##(offset is in the @.Data slot)
  ca=cov_adj(modHet, newdata=dat,specification=spec )

  ##what we want cov_adj to be
  mm0 = model.matrix(modHet)
  mm0[,'ZTRUE'] <- 0
  mm0[,'ZTRUE:x'] <- 0

  myca=plogis(crossprod(t(mm0),coef(modHet))[,1])

  expect_true(all.equal(myca,ca@.Data, check.attributes=FALSE, check.names=FALSE))
})

test_that("two stage lm estimates, SEs reproduce 1 stage as appropriate",
{ ## Rationale for expecting these equivalences is given in 
  ## developer_docs/LinearModelVarianceEstimation.Rmd 
  library(sandwich)
  data(simdata)

  # Binary treatment
  suppressWarnings(spec <- rct_spec(z ~ 1, data = simdata))
  camod <- lm(y ~ x + z, data = simdata)
  ca <- cov_adj(camod, newdata = simdata, specification = spec)
  ddmod  <- lm(y~ z.(), offset=ca, data=simdata)
  expect_equal(coef(ddmod)["z.()"], coef(camod)["z"], ignore_attr=TRUE)
  all.equal(vcovHC(ddmod, type="HC0")["z.()","z.()"],
            vcovHC(camod, type="HC0")["z","z"], 
            check.attributes=FALSE) |>  isTRUE() |> 
    expect_false()# misses camod sampling variability
  ddmod_tee1 <- as.lmitt(ddmod)
  expect_equal(coef(ddmod_tee1)["z.()"], 
               coef(camod)["z"], ignore_attr=TRUE)
### Currently failing:
###  expect_equal(vcov_tee(ddmod_tee1, type="HC0")["z.()","z.()"],
###               vcovHC(camod, type="HC0")["z","z"],
###              ignore_attr=TRUE)

  ddmod_tee2 <- lmitt(y~1, offset=ca, specification=spec, data=simdata)
  expect_equal(coef(ddmod_tee2)["z."], coef(camod)["z"], ignore_attr=TRUE)
### Currently failing:
###  expect_equal(vcov_tee(ddmod_tee2, type="HC0")["z.","z."],
###                 vcovHC(camod, type="HC0")["z","z"], 
###                 ignore_attr=TRUE)
})
