set.seed(20)
cont_x <- rnorm(100)
cat_x <- rbinom(100, 2, 0.2)
y <- cont_x - cat_x + rnorm(100)
xstar <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                    "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                    "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                    "cont_x" = rnorm(100),
                    "cat_x" = rbinom(100, 2, 0.2))

cmod <- lm(y ~ cont_x + as.factor(cat_x))
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

