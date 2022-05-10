set.seed(20)
cont_x <- rnorm(100)
cat_x <- rbinom(100, 2, 0.2)
y <- cont_x - cat_x + rnorm(100)
x <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                "cont_x" = cont_x,
                "cat_x" = cat_x,
                "y" = y)
xstar <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                    "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                    "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                    "cont_x" = rnorm(100),
                    "cat_x" = rbinom(100, 2, 0.2))

des <- rct_design(t ~ uoa(uoa1, uoa2), data = x)
by <- c("teacher" = "uoa1", "student" = "uoa2")
cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x)
offset <- stats::predict(cmod, xstar)
pred_gradient <- model.matrix(cmod)
keys <- xstar[, c("uoa1", "uoa2", "t")]
invalid_cmod <- list("a" = c(1, 2, 3))

expectSandwichLayerError <- function(err_msg,
                                     offset,
                                     fitted_covariance_model,
                                     prediction_gradient,
                                     keys) {
  expect_error(new("SandwichLayer",
                   offset,
                   fitted_covariance_model = fitted_covariance_model,
                   prediction_gradient = prediction_gradient,
                   keys = keys),
               err_msg)
}

test_that("SandwichLayer covariance model incompatible with model.matrix", {
  expectSandwichLayerError("must have a 'terms' attribute",
                           offset, invalid_cmod, pred_gradient, keys)
})

test_that("SandwichLayer covariance model incompatible with sandwich package", {
  on.exit(invalid_cmod <- list("a", "b", "c"))

  invalid_cmod$terms <- cmod$terms
  expectSandwichLayerError("extracting vcov elements not applicable",
                           offset, invalid_cmod, pred_gradient, keys)
})
  
test_that("SandwichLayer prediction gradient is not a numeric matrix", {
  on.exit(pred_gradient <- model.matrix(cmod))

  pred_gradient <- matrix(as.character(cbind(cont_x, cat_x)), ncol = 2)
  expectSandwichLayerError("must be a numeric matrix",
                           offset, cmod, pred_gradient, keys)
})

test_that("SandwichLayer prediction gradient has invalid number of rows", {
  on.exit(pred_gradient <- model.matrix(cmod))

  pred_gradient <- as.matrix(cbind(cont_x, cat_x)[2:100,])
  expectSandwichLayerError("same number of rows",
                           offset, cmod, pred_gradient, keys)
})

test_that("SandwichLayer prediction gradient has invalid number of columns", {
  on.exit(pred_gradient <- model.matrix(cmod))

  pred_gradient <- as.matrix(cont_x)
  expectSandwichLayerError("same number of columns",
                           offset, cmod, pred_gradient, keys)
})

test_that("SandwichLayer prediction gradient has NA's", {
  on.exit(offset <- stats::predict(cmod, xstar))

  offset <- c(offset[1:99], NA_real_)
  expect_warning(new("SandwichLayer",
                     offset,
                     fitted_covariance_model = cmod,
                     prediction_gradient = pred_gradient,
                     keys = keys),
                 "has NA values")
})

test_that("SandwichLayer created correctly", {
  expect_true(is(new("SandwichLayer",
                     offset,
                     fitted_covariance_model = cmod,
                     prediction_gradient = pred_gradient,
                     keys = keys),
                 "SandwichLayer"))
})

test_that("as.SandwichLayer not called with a PreSandwichLayer", {
  expect_error(as.SandwichLayer(offset, des),
               "must be a `PreSandwichLayer`")
})

test_that("as.SandwichLayer not fit with a data argument", {
  on.exit(rm(psl))
  on.exit(cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x))
  
  cmod <- lm(y ~ cont_x + as.factor(cat_x))
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  expect_error(as.SandwichLayer(psl, des),
               "must be fit using a `data` argument")
})

test_that("as.SandwichLayer `by` is not a named vector", {
  on.exit(by <- c("uoa1" = "teacher", "uoa2" = "student"))
  on.exit(rm(psl))

  by <- c("teacher", "student")
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  expect_error(as.SandwichLayer(psl, des, by),
               "must be a named vector")
})

test_that("as.SandwichLayer missing desvars from cov mod data", {
  on.exit(rm(psl))
  on.exit(offset <- stats::predict(cmod, xstar))
  
  x_missing_desvar <- x[, colnames(x)[!(colnames(x) %in% c("uoa1", "uoa2"))]]
  new_cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x_missing_desvar)
  offset <- stats::predict(new_cmod, xstar)
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = new_cmod,
             prediction_gradient = pred_gradient)
  expect_error(as.SandwichLayer(psl, des),
               "columns \"uoa1\", \"uoa2\" are missing")
})

test_that("as.SandwichLayer used correctly", {
  on.exit(des <- rct_design(t ~ uoa(uoa1, uoa2), data = x))

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  expect_true(is(as.SandwichLayer(psl, des), "SandwichLayer"))

  des <- rct_design(t ~ cluster(uoa1, uoa2), data = x)
  expect_true(is(as.SandwichLayer(psl, des), "SandwichLayer"))
  
  des <- rct_design(t ~ unitid(uoa1, uoa2), data = x)
  expect_true(is(as.SandwichLayer(psl, des), "SandwichLayer"))
})

test_that("as.SandwichLayer used correctly with `by`", {
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, des, by)
  expect_true(is(sl, "SandwichLayer"))
  expect_equal(colnames(sl@keys), names(by))
})