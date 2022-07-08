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
cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x)
offset <- stats::predict(cmod, xstar)
pred_gradient <- model.matrix(~ cont_x + as.factor(cat_x), xstar)
keys <- xstar[, c("uoa1", "uoa2")]
invalid_cmod <- list("a" = c(1, 2, 3))

expect_layer_error <- function(layer_class = c("PreSandwichLayer", "SandwichLayer"),
                               err_msg,
                               offset,
                               ...) {
  expect_error(new(layer_class,
                   offset,
                   ...),
               err_msg)
}

test_that("PreSandwichLayer covariance model incompatible with model.matrix", {
  expect_layer_error("PreSandwichLayer",
                     "must have a 'terms' attribute",
                     offset,
                     fitted_covariance_model = invalid_cmod,
                     prediction_gradient = pred_gradient)
})

test_that("PreSandwichLayer covariance model incompatible with sandwich package", {
  on.exit(invalid_cmod <- list("a", "b", "c"))

  invalid_cmod$terms <- cmod$terms
  expect_layer_error("PreSandwichLayer",
                     "extracting vcov elements not applicable",
                     offset,
                     fitted_covariance_model = invalid_cmod,
                     prediction_gradient = pred_gradient)
})

test_that("PreSandwichLayer prediction gradient is not a numeric matrix", {
  on.exit(pred_gradient <- model.matrix(cmod))

  pred_gradient <- matrix(as.character(cbind(cont_x, cat_x)), ncol = 2)
  expect_layer_error("PreSandwichLayer",
                     "must be a numeric matrix",
                     offset,
                     fitted_covariance_model = cmod,
                     prediction_gradient = pred_gradient)
})

test_that("PreSandwichLayer prediction gradient has invalid number of rows", {
  on.exit(pred_gradient <- model.matrix(cmod))

  pred_gradient <- as.matrix(cbind(cont_x, cat_x)[2:100,])
  expect_layer_error("PreSandwichLayer",
                     "same dimension along axis 1",
                     offset,
                     fitted_covariance_model = cmod,
                     prediction_gradient = pred_gradient)
})

test_that("PreSandwichLayer prediction gradient has invalid number of columns", {
  on.exit(pred_gradient <- model.matrix(cmod))

  pred_gradient <- as.matrix(cont_x)
  expect_layer_error("PreSandwichLayer",
                     "same number of columns",
                     offset,
                     fitted_covariance_model = cmod,
                     prediction_gradient = pred_gradient)
})

test_that("SandwichLayer has NA's", {
  on.exit(offset <- stats::predict(cmod, xstar))

  offset <- c(offset[1:99], NA_real_)
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
  on.exit(keys <- xstar[, c("uoa1", "uoa2", "t")])

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  keys <- keys[1:99,]
  expect_layer_error("SandwichLayer",
                     "to fit the covariance model",
                     psl,
                     keys = keys,
                     Design = des)
})

test_that("SandwichLayer created correctly", {
  on.exit(keys <- xstar[, c("uoa1", "uoa2", "t")])

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  expect_true(inherits(new("SandwichLayer",
                     psl,
                     keys = keys,
                     Design = des),
                 "SandwichLayer"))

  keys[100,] <- rep(NA_integer_, ncol(keys))
  expect_true(inherits(new("SandwichLayer",
                     psl,
                     keys = keys,
                     Design = des),
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
  on.exit(rm(psl))

  by <- c("uoa1", "uoa2")
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  expect_error(as.SandwichLayer(psl, des, by),
               "must be named vector")
})

test_that("as.SandwichLayer missing desvar columns from covariance model data", {
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
               "columns uoa1, uoa2 are missing")
})

test_that("as.SandwichLayer used correctly", {
  on.exit(des <- rct_design(t ~ uoa(uoa1, uoa2), data = x))

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  expect_true(inherits(as.SandwichLayer(psl, des), "SandwichLayer"))

  des <- rct_design(t ~ cluster(uoa1, uoa2), data = x)
  expect_true(inherits(as.SandwichLayer(psl, des), "SandwichLayer"))

  des <- rct_design(t ~ unitid(uoa1, uoa2), data = x)
  expect_true(inherits(as.SandwichLayer(psl, des), "SandwichLayer"))
})

test_that("as.SandwichLayer used correctly with `by`", {
  on.exit(x <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                          "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                          "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                          "cont_x" = cont_x,
                          "cat_x" = cat_x,
                          "y" = y))
  on.exit(cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x), add = TRUE)

  colnames(x) <- c("school", "teacher", colnames(x)[-(1:2)])
  cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x)
  by <- c("uoa1" = "school", "uoa2" = "teacher")

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, des, by)
  expect_true(inherits(sl, "SandwichLayer"))
  expect_equal(colnames(sl@keys), as.character(by))
})

test_that("as.SandwichLayer produces NA rows in `keys` for non-NA uoa values", {
  on.exit(x <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                          "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                          "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                          "cont_x" = cont_x,
                          "cat_x" = cat_x,
                          "y" = y))
  on.exit(cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x), add = TRUE)

  x[x$uoa1 == 1 & x$uoa2 == 1, c("uoa1", "uoa2")] <- c(0, 0)
  cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, des)
  expect_equal(nrow(sl@keys),
               nrow(sl@fitted_covariance_model$model))
  expect_true(nrow(sl@keys[is.na(sl@keys$uoa1),]) > 0)
  expect_equal(nrow(x[x$uoa1 == 0 & x$uoa2 == 0,]),
               nrow(sl@keys[is.na(sl@keys$uoa1),]))
  expect_equal(1:25, which(is.na(sl@keys$uoa1)))
})


test_that("as.SandwichLayer produces NA rows in `keys` for NA uoa values", {
  on.exit(x <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                          "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                          "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                          "cont_x" = cont_x,
                          "cat_x" = cat_x,
                          "y" = y))
  on.exit(cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x), add = TRUE)

  x[x$uoa1 == 1 & x$uoa2 == 1, c("uoa1", "uoa2")] <- c(NA, NA)
  cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x)

  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- as.SandwichLayer(psl, des)
  expect_equal(nrow(sl@keys),
               nrow(sl@fitted_covariance_model$model))
  expect_true(nrow(sl@keys[is.na(sl@keys$uoa1),]) > 0)
  expect_equal(nrow(x[is.na(x$uoa1) & is.na(x$uoa2),]),
               nrow(sl@keys[is.na(sl@keys$uoa1),]))
  expect_equal(1:25, which(is.na(sl@keys$uoa1)))
})

test_that("show_sandwich_layer works", {
  on.exit(des <- rct_design(t ~ uoa(uoa1, uoa2), data = x))

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
  psl <- new("PreSandwichLayer",
             offset,
             fitted_covariance_model = cmod,
             prediction_gradient = pred_gradient)
  sl <- new("SandwichLayer",
            psl,
            keys = keys,
            Design = des)

  no_subset_psl <- subset(psl, rep(TRUE, length(offset)))
  no_subset_sl <- subset(sl, rep(TRUE, length(offset)))
  expect_identical(no_subset_psl, psl)
  expect_identical(no_subset_sl, sl)

  expect_true(inherits(no_subset_psl, "PreSandwichLayer"))
  expect_true(inherits(no_subset_sl, "SandwichLayer"))

  expect_identical(capture_output(str(no_subset_psl)),
                   capture_output(str(psl)))
  expect_identical(capture_output(str(no_subset_sl)),
                   capture_output(str(sl)))

  subset_length <- floor(length(offset) / 2)
  subset_psl <- subset(psl, c(rep(TRUE, subset_length),
                              rep(FALSE, length(offset) - subset_length)))
  subset_sl <- subset(sl, c(rep(TRUE, subset_length),
                            rep(FALSE, length(offset) - subset_length)))
  expect_identical(subset_psl@.Data, psl@.Data[1:subset_length])
  expect_identical(subset_sl@.Data, sl@.Data[1:subset_length])

  expect_true(inherits(subset_psl, "PreSandwichLayer"))
  expect_true(inherits(subset_sl, "SandwichLayer"))

  expect_identical(subset_psl@fitted_covariance_model, psl@fitted_covariance_model)
  expect_identical(subset_psl@prediction_gradient, psl@prediction_gradient)
  expect_identical(subset_sl@keys, sl@keys)
  expect_identical(subset_sl@Design, sl@Design)

  no_subset_psl <- psl[]
  no_subset_sl <- sl[]
  expect_identical(no_subset_psl, psl)
  expect_identical(no_subset_sl, sl)

  expect_true(inherits(no_subset_psl, "PreSandwichLayer"))
  expect_true(inherits(no_subset_sl, "SandwichLayer"))

  expect_identical(psl[1:10]@.Data, psl@.Data[1:10])
  expect_identical(sl[1:10]@.Data, sl@.Data[1:10])
})

test_that(".get_ca_and_prediction_gradient newdata is a matrix", {
  expect_error(.get_ca_and_prediction_gradient(cmod, as.matrix(xstar)),
               "must be a dataframe")
})

test_that(".get_ca_and_prediction_gradient model doesn't have a model.matrix method", {
  expect_error(.get_ca_and_prediction_gradient(list("coefficients" = c(1., 1.)),
                                               xstar),
               "must have a `call`")
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "for `lm` object"), {
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod)
  expect_equal(ca_and_grad$ca, cmod$fitted.values)
  expect_equal(ca_and_grad$prediction_gradient, stats::model.matrix(cmod))
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "for `lm` object with new data"), {
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod, xstar)
  expect_equal(ca_and_grad$ca, offset)
  expect_equal(ca_and_grad$prediction_gradient, pred_gradient)
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "when formula is a symbol"), {
  on.exit(cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x))

  cmod_form <- y ~ cont_x + as.factor(cat_x)
  cmod <- lm(cmod_form, data = x)
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod)
  expect_equal(ca_and_grad$ca, cmod$fitted.values)
  expect_equal(ca_and_grad$prediction_gradient, stats::model.matrix(cmod))
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "for `glm` object"), {
  on.exit(cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x))

  cmod <- glm(as.factor(cat_x) ~ cont_x, data = x, family = stats::binomial())
  mm <- stats::model.matrix(cmod)

  ca_and_grad <- .get_ca_and_prediction_gradient(cmod)
  expect_equal(ca_and_grad$ca,
              drop(exp(mm %*% coef(cmod)) / (1 + exp(mm %*% coef(cmod)))))
  expect_equal(ca_and_grad$prediction_gradient,
              cmod$family$mu.eta(drop(exp(mm %*% coef(cmod)) /
                                        (1 + exp(mm %*% coef(cmod))))) * mm)
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output",
                "for `glm` object with new data"), {
  on.exit(cmod <- lm(y ~ cont_x + as.factor(cat_x), data = x))

  cmod <- glm(as.factor(cat_x) ~ cont_x, data = x, family = stats::binomial())
  mm <- stats::model.matrix(cmod, data = xstar)

  ca_and_grad <- .get_ca_and_prediction_gradient(cmod, xstar)
  expect_equal(ca_and_grad$ca,
               drop(exp(mm %*% coef(cmod)) / (1 + exp(mm %*% coef(cmod)))))
  expect_equal(ca_and_grad$prediction_gradient,
               cmod$family$mu.eta(drop(exp(mm %*% coef(cmod)) /
                                         (1 + exp(mm %*% coef(cmod))))) * mm)
})

test_that(paste(".get_ca_and_prediction_gradient returns expected output when",
                "NA's are present"), {
  on.exit(xstar <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                              "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                              "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                              "cont_x" = rnorm(100),
                              "cat_x" = rbinom(100, 2, 0.2)))

  xstar[100, c("cont_x", "cat_x")] <- NA_real_
  ca_and_grad <- .get_ca_and_prediction_gradient(cmod, xstar)
  expect_equal(length(ca_and_grad$ca), 100)
  expect_equal(sum(is.na(ca_and_grad$ca)), 1)
  expect_equal(dim(ca_and_grad$prediction_gradient), dim(pred_gradient[1:99,]))
  expect_equal(ca_and_grad$prediction_gradient[1:99,], pred_gradient[1:99,])
})
