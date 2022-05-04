test_that("Sandwich Layer validity", {
  set.seed(20)
  x <- rnorm(100)
  y <- x + rnorm(100)
  xstar <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                      "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                      "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                      "x" = rnorm(100))
  
  valid_lm_cmod <- lm(y ~ x)
  invalid_cmod <- list("a" = c(1, 2, 3))
  offset <- stats::predict(valid_lm_cmod, xstar)
  pred_gradient <- x
  keys <- cbind(xstar[, c("uoa1", "uoa2", "t")], "offset" = offset)
  
  expect_error(new("SandwichLayer",
                   offset,
                   fitted_covariance_model = invalid_cmod,
                   prediction_gradient = pred_gradient,
                   keys = keys),
               "must have a 'terms' attribute")
})