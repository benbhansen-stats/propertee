test_that("Sandwich Layer validity", {
  set.seed(20)
  cont_x <- rnorm(100)
  cat_x <- rbinom(100, 2, 0.2)
  y <- cont_x - cat_x + rnorm(100)
  xstar <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                      "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                      "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                      "cont_x" = rnorm(100),
                      "cat_x" = rbinom(100, 2, 0.2))
  
  valid_cmod <- lm(y ~ cont_x + as.factor(cat_x))
  invalid_cmod <- list("a" = c(1, 2, 3))
  offset <- stats::predict(valid_cmod, xstar)
  pred_gradient <- as.matrix(cbind(cont_x, cat_x))
  keys <- cbind(xstar[, c("uoa1", "uoa2", "t")], "offset" = offset)
  
  # covariance model incompatible with model.matrix
  expect_error(new("SandwichLayer",
                   offset,
                   fitted_covariance_model = invalid_cmod,
                   prediction_gradient = pred_gradient,
                   keys = keys),
               "must have a 'terms' attribute")
  
  # covariance model incompatible with sandwich package functions
  invalid_cmod$terms <- valid_cmod$terms
  expect_error(new("SandwichLayer",
                   offset,
                   fitted_covariance_model = invalid_cmod,
                   prediction_gradient = pred_gradient,
                   keys = keys),
               "extracting vcov elements not applicable")
  
  # prediction gradient is not a numeric matrix
  expect_error(new("SandwichLayer",
                   offset,
                   fitted_covariance_model = valid_cmod,
                   prediction_gradient = matrix(as.character(cbind(cont_x, cat_x)), ncol = 2),
                   keys = keys),
               "must be a numeric matrix")
  
  # prediction gradient does not have valid number of rows
  expect_error(new("SandwichLayer",
                   offset,
                   fitted_covariance_model = valid_cmod,
                   prediction_gradient = as.matrix(cbind(cont_x, cat_x)[2:100,]),
                   keys = keys),
               "same number of rows")
  
})