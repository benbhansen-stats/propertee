test_that("PreSandwichLayer covariance model incompatible with model.matrix", {
  N <- 100
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  invalid_cmod <- list("a" = c(1, 2, 3),
                       "terms" = list("b" = c(4, 5, 6)))

  expect_error(new("PreSandwichLayer",
                   offset,
                   fitted_covariance_model = invalid_cmod,
                   prediction_gradient = pred_gradient),
               "must have a valid 'terms' attribute")
})

test_that("PreSandwichLayer covariance model incompatible with sandwich package", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  cmod <- lm(y ~ x, df)
  class(cmod) <- "new_lm"

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  expect_error(new("PreSandwichLayer",
                   offset,
                   fitted_covariance_model = cmod,
                   prediction_gradient = pred_gradient),
               "extracting vcov elements not applicable")
})

test_that("PreSandwichLayer prediction gradient is not a numeric matrix", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  cmod <- lm(y ~ x, df)

  offset <- rep(1, N)
  pred_gradient <- as.matrix(cbind(1, as.character(df$x)))

  expect_error(new("PreSandwichLayer",
                   offset,
                   fitted_covariance_model = cmod,
                   prediction_gradient = pred_gradient),
               "must be a numeric matrix")
})

test_that("PreSandwichLayer prediction gradient has invalid number of rows", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  cmod <- lm(y ~ x, df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N - 1, ncol = 2)

  expect_error(new("PreSandwichLayer",
                   offset,
                   fitted_covariance_model = cmod,
                   prediction_gradient = pred_gradient),
               "same dimension along axis 1")
})

test_that("PreSandwichLayer prediction gradient has invalid number of columns", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  cmod <- lm(y ~ x, df)
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 1)

  expect_error(new("PreSandwichLayer",
                   offset,
                   fitted_covariance_model = cmod,
                   prediction_gradient = pred_gradient),
               "same number of columns")
})

test_that("SandwichLayer has NA's", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  des <- rct_design(z ~ unitid(uid), df)
  keys <- df[, "uid", drop = FALSE]
  
  offset <- rep(1, N)
  offset[N] <- NA_real_
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  pred_gradient[N,] <- NA_real_
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  
  expect_warning(new("SandwichLayer",
                     psl,
                     keys = keys,
                     Design = des),
                 "adjustments are NA")
})

test_that("SandwichLayer keys doesn't have the same row count as covariance model data", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  des <- rct_design(z ~ unitid(uid), df)
  keys <- df[seq_len(N-1), "uid", drop = FALSE]
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  expect_error(new("SandwichLayer",
                   psl,
                   keys = keys,
                   Design = des),
               "to fit the covariance model")
})

test_that("SandwichLayer created correctly", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  des <- rct_design(z ~ unitid(uid), df)
  keys <- df[, "uid", drop = FALSE]
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  expect_true(inherits(new("SandwichLayer",
                     psl,
                     keys = keys,
                     Design = des),
                 "SandwichLayer"))

  keys$uid <- NA_integer_

  expect_true(inherits(new("SandwichLayer",
                     psl,
                     keys = keys,
                     Design = des),
                 "SandwichLayer"))
})

test_that("as.SandwichLayer not called with a PreSandwichLayer", {
  set.seed(20)
  N <- 100
  df <- data.frame("z" = rbinom(N, 1, 0.5), "uid" = seq_len(N))
  des <- rct_design(z ~ unitid(uid), df)
  
  expect_error(as.SandwichLayer(seq_len(N), des),
               "must be a `PreSandwichLayer`")
})

test_that("as.SandwichLayer not fit with a data argument", {
  set.seed(20)
  N <- 100
  x <- rnorm(N); y <- rnorm(N)
  design_df <- data.frame("z" = rbinom(N, 1, 0.5), "uid" = seq_len(N))
  cmod <- lm(y ~ x)
  des <- rct_design(z ~ unitid(uid), design_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  expect_error(as.SandwichLayer(psl, des),
               "must be fit using a `data` argument")
})

test_that("as.SandwichLayer `by` is not a named vector", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  des <- rct_design(z ~ unitid(uid), df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  by <- c("unitid")

  expect_error(as.SandwichLayer(psl, des, by),
               "must be named vector")
})

test_that("as.SandwichLayer missing desvar columns from covariance model data", {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  design_df <- data.frame("z" = rbinom(N, 1, 0.5), "uid" = seq_len(N))
  cmod <- lm(y ~ x, cmod_df)
  des <- rct_design(z ~ unitid(uid), design_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  expect_error(as.SandwichLayer(psl, des),
               "columns uid are missing")
})

test_that("as.SandwichLayer used correctly", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  des <- rct_design(z ~ unit_of_assignment(uid), df)
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  expect_true(inherits(as.SandwichLayer(psl, des), "SandwichLayer"))

  des <- rct_design(z ~ cluster(uid), df)
  expect_true(inherits(as.SandwichLayer(psl, des), "SandwichLayer"))

  des <- rct_design(z ~ unitid(uid), df)
  expect_true(inherits(as.SandwichLayer(psl, des), "SandwichLayer"))
})

test_that("as.SandwichLayer used correctly with `by`", {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "uid" = seq_len(N))
  cmod <- lm(y ~ x, cmod_df)
  design_df <- data.frame("uoa1" = seq_len(N), "z" = rbinom(N, 1, 0.5))
  des <- rct_design(z ~ unit_of_assignment(uoa1), design_df)
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  by <- c("uoa1" = "uid")
  sl <- as.SandwichLayer(psl, des, by)
  
  expect_true(inherits(sl, "SandwichLayer"))
  expect_equal(colnames(sl@keys), as.character(by))
})

test_that(paste("as.SandwichLayer produces NA rows in `keys` for non-NA uoa",
                "values not found in the design"), {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "uid" = seq_len(N) - N/4)
  cmod <- lm(y ~ x, cmod_df)
  design_df <- data.frame("uid" = seq_len(N), "z" = rbinom(N, 1, 0.5))
  des <- rct_design(z ~ unit_of_assignment(uid), design_df)
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, des)
  
  expect_equal(nrow(sl@keys),
               nrow(sl@fitted_covariance_model$model))
  expect_true(nrow(sl@keys[is.na(sl@keys$uid), , drop = FALSE]) > 0)
  expect_equal(nrow(sl@keys[is.na(sl@keys$uid), , drop = FALSE]),
               nrow(cmod_df[cmod_df$uid <= 0,]))
  expect_equal(which(is.na(sl@keys$uid)), seq_len(25))
})


test_that("as.SandwichLayer produces NA rows in `keys` for NA uoa values", {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N),
                        "uid" = c(rep(NA_integer_, N/4), seq_len(N - N/4)))
  cmod <- lm(y ~ x, cmod_df)
  design_df <- data.frame("uid" = seq_len(N), "z" = rbinom(N, 1, 0.5))
  des <- rct_design(z ~ unit_of_assignment(uid), design_df)
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, des)

  expect_equal(nrow(sl@keys),
               nrow(sl@fitted_covariance_model$model))
  expect_true(nrow(sl@keys[is.na(sl@keys$uid), , drop = FALSE]) > 0)
  expect_equal(nrow(sl@keys[is.na(sl@keys$uid), , drop = FALSE]),
               nrow(cmod_df[is.na(cmod_df$uid),]))
  expect_equal(which(is.na(sl@keys$uid)), seq_len(25))
})

test_that("show_sandwich_layer works", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  des <- rct_design(z ~ unitid(uid), df)
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  out <- capture.output(show(as.vector(psl)))
  caout <- capture.output(show(psl))
  expect_identical(out, caout)

  sl <- as.SandwichLayer(psl, des)
  out <- capture.output(show(as.vector(sl)))
  caout <- capture.output(show(sl))
  expect_identical(out, caout)
})

test_that("subsetting PreSandwich and SandwichLayer works", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  des <- rct_design(z ~ unitid(uid), df)
  
  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  keys <- df[, "uid", drop = FALSE]
  sl <- new("SandwichLayer",
            psl,
            keys = keys,
            Design = des)

  no_subset_psl <- subset(psl, rep(TRUE, length(offset)))
  no_subset_sl <- subset(sl, rep(TRUE, length(offset)))
  expect_identical(no_subset_psl@.Data, psl@.Data)
  expect_identical(no_subset_psl@prediction_gradient,
                   psl@prediction_gradient[1:100,])
  expect_identical(no_subset_sl@.Data, sl@.Data)
  expect_identical(no_subset_sl@prediction_gradient,
                   sl@prediction_gradient[1:100,])

  expect_true(inherits(no_subset_psl, "PreSandwichLayer"))
  expect_true(inherits(no_subset_sl, "SandwichLayer"))

  expect_identical(no_subset_psl@fitted_covariance_model,
                   psl@fitted_covariance_model)
  expect_identical(no_subset_sl@fitted_covariance_model,
                   sl@fitted_covariance_model)

  subset_length <- floor(length(offset) / 2)
  subset_psl <- subset(psl, c(rep(TRUE, subset_length),
                              rep(FALSE, length(offset) - subset_length)))
  subset_sl <- subset(sl, c(rep(TRUE, subset_length),
                            rep(FALSE, length(offset) - subset_length)))
  expect_identical(subset_psl@.Data, psl@.Data[1:subset_length])
  expect_identical(subset_sl@.Data, sl@.Data[1:subset_length])

  expect_true(inherits(subset_psl, "PreSandwichLayer"))
  expect_true(inherits(subset_sl, "SandwichLayer"))

  expect_identical(subset_psl@fitted_covariance_model,
                   psl@fitted_covariance_model)
  expect_identical(subset_psl@prediction_gradient,
                   psl@prediction_gradient[1:subset_length,])
  expect_identical(subset_sl@keys, sl@keys)
  expect_identical(subset_sl@Design, sl@Design)

  no_subset_psl <- psl[]
  no_subset_sl <- sl[]
  expect_identical(no_subset_psl@.Data, psl@.Data)
  expect_identical(no_subset_sl@.Data, sl@.Data)

  expect_true(inherits(no_subset_psl, "PreSandwichLayer"))
  expect_true(inherits(no_subset_sl, "SandwichLayer"))

  expect_identical(psl[1:10]@.Data, psl@.Data[1:10])
  expect_identical(sl[1:10]@.Data, sl@.Data[1:10])
  expect_identical(psl[1:10]@prediction_gradient,
                   psl@prediction_gradient[1:10,])
  expect_identical(sl[1:10]@prediction_gradient,
                   sl@prediction_gradient[1:10,])
})

test_that(".get_ca_and_prediction_gradient newdata is a matrix", {
  expect_error(.get_ca_and_prediction_gradient("a", matrix(1)),
               "must be a dataframe")
})

test_that(".get_ca_and_prediction_gradient model doesn't have a terms method", {
  expect_error(.get_ca_and_prediction_gradient(list("coefficients" = c(1., 1.)),
                                               data.frame("x" = 1)),
               "must have `terms`")
})

test_that(".get_ca_and_prediction_gradient model frame missing cmod columns", {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "x2" = rnorm(N), "y" = rnorm(N))
  pred_df <- data.frame("x1" = rnorm(N))
  cmod <- lm(y ~ x1 + x2, df)
  
  expect_error(.get_ca_and_prediction_gradient(cmod, pred_df),
               "'x2' not found")
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "for `lm` object"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  cmod <- lm(y ~ x, df)
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod)

  expect_equal(ca_and_grad$ca, cmod$fitted.values)
  expect_equal(ca_and_grad$prediction_gradient, stats::model.matrix(cmod))
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "for `lm` object with new data"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  pred_df <- data.frame("x" = rnorm(N))
  cmod <- lm(y ~ x, df)
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod, pred_df)

  expect_equal(ca_and_grad$ca,
               drop(stats::model.matrix(
                      formula(stats::delete.response(terms(cmod))), pred_df) %*%
                    cmod$coefficients))
  expect_equal(ca_and_grad$prediction_gradient,
               stats::model.matrix(
                 formula(stats::delete.response(terms(cmod))), pred_df))
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "when formula is a symbol"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  cmod_form <- y ~ x
  cmod <- lm(cmod_form, df)
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod)
  
  expect_equal(ca_and_grad$ca, cmod$fitted.values)
  expect_equal(ca_and_grad$prediction_gradient, stats::model.matrix(cmod))
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "when variables are modified in the formula"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "x2" = rpois(N, 1) + 1, "y" = rnorm(N))
  cmod_form <- y ~ stats::poly(x1, 3) + log(x2)
  cmod <- lm(cmod_form, data = df)
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod)
  
  expect_equal(ca_and_grad$ca, cmod$fitted.values)
  expect_equal(ca_and_grad$prediction_gradient, stats::model.matrix(cmod))
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "for `glm` object"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "x2" = rpois(N, 1) + 1, "y" = rnorm(N))
  cmod <- glm(x2 ~ x1, data = df, family = stats::poisson())
  mm <- stats::model.matrix(cmod)
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod)
  
  expect_equal(ca_and_grad$ca,
              drop(exp(mm %*% cmod$coefficients)))
  expect_equal(ca_and_grad$prediction_gradient,
              cmod$family$mu.eta(drop(exp(mm %*% cmod$coefficients))) * mm)
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "for `glm` object with new data"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "x2" = rpois(N, 1) + 1, "y" = rnorm(N))
  cmod <- glm(x2 ~ x1, data = df, family = stats::poisson())
  pred_df <- data.frame("x1" = rnorm(N))
  mm <- stats::model.matrix(formula(stats::delete.response(terms(cmod))),
                            data = pred_df)
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod, pred_df)
  
  expect_equal(ca_and_grad$ca,
               drop(exp(mm %*% cmod$coefficients)))
  expect_equal(ca_and_grad$prediction_gradient,
               cmod$family$mu.eta(drop(exp(mm %*% cmod$coefficients))) * mm)
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output when",
                "NA's are present"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  pred_df <- data.frame("x" = rnorm(N))
  cmod <- lm(y ~ x, df)

  pred_df[N, c("x")] <- NA_real_
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod, pred_df)
  pred_gradient <- stats::model.matrix(
    formula(stats::delete.response(terms(cmod))),
    stats::model.frame(pred_df, na.action = na.pass))
  
  expect_equal(length(ca_and_grad$ca), N)
  expect_equal(sum(is.na(ca_and_grad$ca)), 1)
  expect_equal(dim(ca_and_grad$prediction_gradient), dim(pred_gradient))
  expect_equal(ca_and_grad$prediction_gradient[seq_len(N-1),],
               pred_gradient[seq_len(N-1),])
})
