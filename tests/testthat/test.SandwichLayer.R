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
  spec <- rct_spec(z ~ unitid(uid), df)
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
                     StudySpecification = spec),
                 "adjustments are NA")
})

test_that("SandwichLayer keys doesn't have the same row count as covariance model data", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  spec <- rct_spec(z ~ unitid(uid), df)
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
                   StudySpecification = spec),
               "to fit the covariance adjustment model")
})

test_that("SandwichLayer created correctly", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  spec <- rct_spec(z ~ unitid(uid), df)
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
                     StudySpecification = spec),
                 "SandwichLayer"))

  keys$uid <- NA_integer_

  expect_true(inherits(new("SandwichLayer",
                     psl,
                     keys = keys,
                     StudySpecification = spec),
                 "SandwichLayer"))
})

test_that("as.SandwichLayer not called with a PreSandwichLayer", {
  set.seed(20)
  N <- 100
  df <- data.frame("z" = rbinom(N, 1, 0.5), "uid" = seq_len(N))
  spec <- rct_spec(z ~ unitid(uid), df)

  expect_error(as.SandwichLayer(seq_len(N), spec),
               "must be a `PreSandwichLayer`")
})

test_that("as.SandwichLayer not fit with a data argument", {
  set.seed(20)
  N <- 100
  x <- rnorm(N); y <- rnorm(N)
  specification_df <- data.frame("z" = rbinom(N, 1, 0.5), "uid" = seq_len(N))
  cmod <- lm(y ~ x)
  spec <- rct_spec(z ~ unitid(uid), specification_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  expect_error(as.SandwichLayer(psl, spec),
               "must be fit using a `data` argument")
})

test_that("as.SandwichLayer missing specvar columns from covariance model data", {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  specification_df <- data.frame("z" = rbinom(N, 1, 0.5), "uid" = seq_len(N))
  cmod <- lm(y ~ x, cmod_df)
  spec <- rct_spec(z ~ unitid(uid), specification_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  expect_error(as.SandwichLayer(psl, spec),
               "Columns uid are missing")
})

test_that("as.SandwichLayer used correctly with NULL `by` and `Q_data`", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  spec <- rct_spec(z ~ unit_of_assignment(uid), df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  expect_true(inherits(sl1 <- as.SandwichLayer(psl, spec), "SandwichLayer"))
  expect_equal(length(setdiff(colnames(sl1@keys), c("uid", "in_Q"))), 0)
  expect_true(all(sl1@keys$in_Q == 1))
  expect_equal(sl1@keys$uid, df$uid)

  spec <- rct_spec(z ~ cluster(uid), df)
  expect_true(inherits(sl2 <- as.SandwichLayer(psl, spec), "SandwichLayer"))
  expect_equal(length(setdiff(colnames(sl2@keys), c("uid", "in_Q"))), 0)
  expect_true(all(sl2@keys$in_Q == 1))
  expect_equal(sl2@keys$uid, df$uid)

  spec <- rct_spec(z ~ unitid(uid), df)
  expect_true(inherits(sl3 <- as.SandwichLayer(psl, spec), "SandwichLayer"))
  expect_equal(length(setdiff(colnames(sl3@keys), c("uid", "in_Q"))), 0)
  expect_true(all(sl3@keys$in_Q == 1))
  expect_equal(sl3@keys$uid, df$uid)
})

test_that("as.SandwichLayer matches on `by` column when uoa columns don't (`by` has no names)", {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame(x = rnorm(N), y = rnorm(N), uid = seq_len(N), clust = rep(NA_integer_, N))
  cmod <- lm(y ~ x, cmod_df)
  specification_df <- data.frame(uid = seq_len(N), clust = rep(c(1, 2), each = N/2),
                          z = rep(c(0, 1), each = N/2))
  spec <- rct_spec(z ~ unit_of_assignment(clust), specification_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, spec, "uid", specification_df)

  expect_true(inherits(sl, "SandwichLayer"))
  expect_equal(length(setdiff(colnames(sl@keys), c("clust", "uid", "in_Q"))), 0)
  expect_true(all(sl@keys$in_Q))
})

test_that("as.SandwichLayer matches on `by` column when uoa columns don't (`by` has names)", {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame(x = rnorm(N), y = rnorm(N), uid = seq_len(N), clust = rep(NA_integer_, N))
  cmod <- lm(y ~ x, cmod_df)
  specification_df <- data.frame(uoa1 = seq_len(N), clust = rep(c(1, 2), each = N/2),
                          z = rep(c(0, 1), each = N/2))
  spec <- rct_spec(z ~ cluster(clust), specification_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  by <- c("uoa1" = "uid")
  sl <- as.SandwichLayer(psl, spec, by, specification_df)

  expect_true(inherits(sl, "SandwichLayer"))
  expect_equal(length(setdiff(colnames(sl@keys), c("clust", "uid", "in_Q"))), 0)
  expect_true(all(sl@keys$in_Q))
})

test_that("as.SandwichLayer used correctly with unnamed `by` and NULL `Q_data`", {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "uoa1" = seq_len(N))
  cmod <- lm(y ~ x, cmod_df)
  specification_df <- data.frame("uoa1" = seq_len(N), "z" = rbinom(N, 1, 0.5))
  spec <- rct_spec(z ~ unit_of_assignment(uoa1), specification_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  expect_warning(
    expect_warning(
      expect_warning(
        sl <- as.SandwichLayer(psl, spec, "uoa1"),
        "No call to"
      ),
      "Unable to detect"
    ),
    "Could not find direct adjustment data"
  )

  expect_true(inherits(sl, "SandwichLayer"))
  expect_equal(length(setdiff(colnames(sl@keys), c("uoa1", "in_Q"))), 0)
  expect_true(all(sl@keys$in_Q == 1))
  expect_equal(sl@keys$uoa1, cmod_df$uoa1)
})

test_that("as.SandwichLayer failed by", {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "uoa1" = seq_len(N))
  cmod <- lm(y ~ x, cmod_df)
  specification_df <- data.frame("uoa1" = seq_len(N), "z" = rbinom(N, 1, 0.5))
  spec <- rct_spec(z ~ unit_of_assignment(uoa1), specification_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  expect_error(
    expect_warning(
      expect_warning(
        expect_warning(
          as.SandwichLayer(psl, specification = spec, by = "not_uoa"),
          "No call to"
        ),
        "Unable to detect"
      ),
      "Could not find direct adjustment data"
    ),
    "Could not find columns"
  )
})

test_that(paste("as.SandwichLayer produces correct ID's for univariate uoa ID's",
                "not starting at 1"), {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N),
                        "uid" = c(paste0("0400", seq_len(floor(N/2))),
                                  rep(NA_character_, N - floor(N/2))))
  cmod <- lm(y ~ x, cmod_df)
  specification_df <- data.frame("uid" = paste0("0400", seq_len(N)),
                          "z" = rbinom(N, 1, 0.5))
  spec <- rct_spec(z ~ unit_of_assignment(uid), specification_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, spec)

  expect_equal(nrow(sl@keys),
               nrow(sl@fitted_covariance_model$model))
  expect_equal(sl@keys[, 1],
               c(paste0("0400", seq_len(floor(N/2))),
                 rep(NA_character_, N - floor(N/2))))
})

test_that(paste("as.SandwichLayer produces correct ID's for univariate uoa ID's",
                "not starting at 1"), {
  set.seed(20)
  N <- 100
  cmod_df <- data.frame("x" = rnorm(N), "y" = rnorm(N),
                        "classid" = c(rep(c(1, 2), each = floor(N/4)),
                                      rep(NA_character_, N - floor(N/2))),
                        "schoolid" = c(rep("04001", 2 * floor(N/4)),
                                       rep(NA_character_, N - floor(N/2))))
  cmod <- lm(y ~ x, cmod_df)
  specification_df <- data.frame("classid" = rep(c(1, 2), each = floor(N/4)),
                          "schoolid" = rep("04001", 2 * floor(N/4)),
                          "z" = rep(c(0, 1), each = floor(N/4)))
  spec <- rct_spec(z ~ cluster(classid, schoolid), specification_df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, spec)

  expected_keys <- c(rep(paste0("04001_", c(1, 2)), each = floor(N/4)),
                     rep(NA_character_, N - 2 * floor(N/4)))

  expect_equal(nrow(sl@keys),
               nrow(sl@fitted_covariance_model$model))
  expect_equal(sl@keys[, "classid"],
               c(rep(c(1, 2), each = floor(N/4)),
                 rep(NA_character_, N - floor(N/2))))
  expect_equal(sl@keys[, "schoolid"],
               c(rep("04001", 2 * floor(N/4)),
                 rep(NA_character_, N - floor(N/2))))
})

test_that("show_sandwich_layer works", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N), "z" = rbinom(N, 1, 0.5),
                   "uid" = seq_len(N))
  cmod <- lm(y ~ x, df)
  spec <- rct_spec(z ~ unitid(uid), df)

  offset <- rep(1, N)
  pred_gradient <- matrix(1, nrow = N, ncol = 2)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)

  out <- capture.output(show(as.vector(psl)))
  caout <- capture.output(show(psl))
  expect_identical(out, caout)

  sl <- as.SandwichLayer(psl, spec)
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
  spec <- rct_spec(z ~ unitid(uid), df)

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
            StudySpecification = spec)

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
  expect_identical(subset_sl@StudySpecification, sl@StudySpecification)

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

test_that(".make_PreSandwichLayer newdata is a matrix", {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "y" = rnorm(N))
  cmod <- lm(y ~ x1, df)
  expect_error(.make_PreSandwichLayer(cmod, matrix(1)),
               "must be a data.frame")
})

test_that(".make_PreSandwichLayer model doesn't have a terms method", {
  expect_error(.make_PreSandwichLayer(list("coefficients" = c(1., 1.)),
                                               data.frame("x" = 1)),
               "must have `terms`")
})

test_that(".make_PreSandwichLayer model frame missing cmod columns", {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "x2" = rnorm(N), "y" = rnorm(N))
  pred_df <- data.frame("x1" = rnorm(N))
  cmod <- lm(y ~ x1 + x2, df)

  expect_error(.make_PreSandwichLayer(cmod, pred_df),
               "'x2' not found")
})

test_that(".make_PreSandwichLayer returns expected output for `lm` object", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  cmod <- lm(y ~ x, df)
  psl <- .make_PreSandwichLayer(cmod, df)

  expect_true(all.equal(psl@.Data, cmod$fitted.values, check.attributes = FALSE))
  expect_true(all.equal(psl@prediction_gradient, stats::model.matrix(cmod),
                        check.attributes = FALSE))
})

test_that(".make_PreSandwichLayer returns expected output for `lm` object with new data", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  pred_df <- data.frame("x" = rnorm(N))
  cmod <- lm(y ~ x, df)
  psl <- .make_PreSandwichLayer(cmod, pred_df)

  expect_true(all.equal(
    psl@.Data,
    drop(stats::model.matrix(formula(stats::delete.response(terms(cmod))), pred_df) %*%
           cmod$coefficients),
    check.attributes = FALSE
  ))
  expect_true(all.equal(
    psl@prediction_gradient,
    stats::model.matrix(formula(stats::delete.response(terms(cmod))), pred_df)
  ))
})

test_that(".make_PreSandwichLayer returns expected output when formula is a symbol", {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  cmod_form <- y ~ x
  cmod <- lm(cmod_form, df)
  psl <- .make_PreSandwichLayer(cmod, df)

  expect_true(all.equal(psl@.Data, cmod$fitted.values, check.attributes = FALSE))
  expect_true(all.equal(psl@prediction_gradient, stats::model.matrix(cmod),
                        check.attributes = FALSE))
})

test_that(paste(".make_PreSandwichLayer returns expected output",
                "when variables are modified in the formula"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "x2" = rpois(N, 1) + 1, "y" = rnorm(N))
  cmod_form <- y ~ stats::poly(x1, 3) + log(x2)
  cmod <- lm(cmod_form, data = df)
  psl <- .make_PreSandwichLayer(cmod, df)

  expect_true(all.equal(psl@.Data, cmod$fitted.values, check.attributes = FALSE))
  expect_true(all.equal(psl@prediction_gradient, stats::model.matrix(cmod),
                        check.attributes = FALSE))
})

test_that(".make_PreSandwichLayer returns expected output for `glm` object", {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "x2" = rpois(N, 1) + 1, "y" = rnorm(N))
  cmod <- glm(x2 ~ x1, data = df, family = stats::poisson())
  mm <- stats::model.matrix(cmod)
  psl <- .make_PreSandwichLayer(cmod, df)

  expect_true(all.equal(psl@.Data, drop(exp(mm %*% cmod$coefficients)),
                        check.attributes = FALSE))
  expect_true(all.equal(psl@prediction_gradient,
                        cmod$family$mu.eta(drop(mm %*% cmod$coefficients)) * mm,
                        check.attributes = FALSE))
})

test_that(paste(".make_PreSandwichLayer returns expected output",
                "for `glm` object with new data"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x1" = rnorm(N), "x2" = rpois(N, 1) + 1, "y" = rnorm(N))
  cmod <- glm(x2 ~ x1, data = df, family = stats::poisson())
  pred_df <- data.frame("x1" = rnorm(N))
  mm <- stats::model.matrix(formula(stats::delete.response(terms(cmod))),
                            data = pred_df)
  psl <- .make_PreSandwichLayer(cmod, pred_df)

  expect_true(all.equal(psl@.Data, drop(exp(mm %*% cmod$coefficients)),
                        check.attributes = FALSE))
  expect_true(all.equal(psl@prediction_gradient,
                        cmod$family$mu.eta(drop(mm %*% cmod$coefficients)) * mm,
                        check.attributes = FALSE))
})

if (requireNamespace("robustbase", quietly = TRUE)) {
  test_that(".make_PreSandwichLayer glmrob", {
    set.seed(30)
    x = sample(rep(c(-1, 1), each = 15))
    y = c(0.5 * x[1:29] - 2 + 1e-2 * rnorm(29), 12 * x[30] + 1e-2 * rnorm(1))
    moddata <- data.frame(x = x, y = y)

    suppressWarnings(
      mod <- robustbase::glmrob(y ~ x, moddata, family = stats::gaussian(),
                                control = robustbase::glmrobMqle.control(tcc = 1.5))
    )

    expect_true(!all(mod$w.r == 1.))
    expect_true(all.equal(.make_PreSandwichLayer(mod, moddata)@.Data,
                          mod$fitted.values,
                          check.attributes = FALSE))
  })
}

test_that(paste(".make_PreSandwichLayer returns expected output when",
                "NA's are present"), {
  set.seed(20)
  N <- 100
  df <- data.frame("x" = rnorm(N), "y" = rnorm(N))
  pred_df <- data.frame("x" = rnorm(N))
  cmod <- lm(y ~ x, df)

  pred_df[N, c("x")] <- NA_real_
  psl <- .make_PreSandwichLayer(cmod, pred_df)
  pred_gradient <- stats::model.matrix(
    formula(stats::delete.response(terms(cmod))),
    stats::model.frame(pred_df, na.action = na.pass))

  expect_equal(length(psl@.Data), N)
  expect_equal(sum(is.na(psl@.Data)), 1)
  expect_equal(dim(psl@prediction_gradient), dim(pred_gradient))
  expect_true(all.equal(psl@prediction_gradient[seq_len(N-1),],
                        pred_gradient[seq_len(N-1),],
                        check.attributes = FALSE))
})

test_that(".sanitize_C_ids fails with invalid `cluster` argument", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  ssmod <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))

  expect_error(.sanitize_C_ids(ssmod$model$`(offset)`, id_col = "uid"),
               "uid could not be found")
})

test_that(".sanitize_C_ids with full UOA info", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  ssmod <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))

  ids <- .sanitize_C_ids(ssmod$model$`(offset)`)
  expected_ids <- apply(simdata[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))
  expect_equal(ids, expected_ids)
})

test_that(".sanitize_C_ids with partial UOA info", {
  data(simdata)
  cmod_data <- data.frame("x" = rnorm(10), "y" = rnorm(10),
                          "uoa1" = rep(c(1, 2), each = 5),  "uoa2" = NA)

  cmod <- lm(y ~ x, cmod_data)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  ssmod <- lmitt(y ~ 1, data = simdata, specification = spec,
                offset = cov_adj(cmod))

  ids <- .sanitize_C_ids(ssmod$model$`(offset)`)
  expect_equal(length(ids), nrow(cmod_data))
  expect_equal(length(unique(ids)), 2)
})

test_that(".sanitize_C_ids with no UOA info", {
  data(simdata)
  cmod_data <- data.frame("x" = rnorm(10), "y" = rnorm(10), "uoa1" = NA,  "uoa2" = NA)

  cmod <- lm(y ~ x, cmod_data)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), simdata)
  ssmod <- lmitt(y ~ 1, data = simdata, specification = spec,
                offset = cov_adj(cmod))

  ids <- .sanitize_C_ids(ssmod$model$`(offset)`)
  expect_equal(length(ids), nrow(cmod_data))
  expect_equal(length(unique(ids)), 1)
})

test_that(".sanitize_C_ids miscellaneous errors", {
  expect_error(.sanitize_C_ids(2), "x must be a `SandwichLayer`")

  n <- 10
  df <- data.frame("x" = rnorm(n), "a" = rep(c(0, 1), each = 5), "y" = rnorm(n),
                   "cid" = sample(seq_len(n)))
  cmod <- lm(y ~ x, df)
  spec <- rct_spec(a ~ cluster(cid), df)
  sl <- cov_adj(cmod, newdata = df, specification = spec)
  num_C_ids <- .sanitize_C_ids(sl, sorted = TRUE)
  expect_true(all.equal(num_C_ids$x, seq_len(n), check.attributes = FALSE))

  df$cid <- sample(letters[1:n])
  cmod <- lm(y ~ x, df)
  spec <- rct_spec(a ~ cluster(cid), df)
  sl <- cov_adj(cmod, newdata = df, specification = spec)
  char_C_ids <- .sanitize_C_ids(sl, sorted = TRUE)
  expect_true(all.equal(char_C_ids$x, letters[1:n], check.attributes = FALSE))
})
