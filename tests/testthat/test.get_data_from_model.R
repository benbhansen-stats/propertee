test_that("Obtaining data for weights", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata, weights = ate(des))

  #Handling multiple models in the stack
  mod1 <- lm(predict(lm(y~x, data = simdata, weights = ate(des))) ~ 1)
  mod2 <- lm(predict(lm(y~x, data = simdata)) ~ 1, weights = ate(des),
             data = simdata)
  expect_true(mod1$coefficients[1] != mod2$coefficients[1])

  # linear glm and lm are equivalent
  mod3 <- glm(y ~ x, data = simdata, weights = ate(des))
  expect_equal(mod$coefficients,
               mod3$coefficients)

  # Works using `with`
  mod4 <- with(simdata,
               lm(y ~ x, weights = ate(des)))
  expect_equal(mod$coefficients,
               mod4$coefficients)

  # Works using `attach`
  attach(simdata)
  mod5 <- lm(y ~ x, weights = ate(des))
  expect_equal(mod$coefficients,
               mod5$coefficients)
  detach(simdata)

  # Allow simple manipulation of weights
  mod6 <- lm(y ~ x, data = simdata, weights = sqrt(ate(des)))


  # Test for fallback if no model.frame call

  w <- capture_warnings(with(simdata, ate(des)))
  expect_equal(length(w), 4)
  expect_true(any(sapply(w, grepl, pattern = "No call")))
  expect_true(any(sapply(w, grepl, pattern = "trying fallback")))

  expect_error(.get_data_from_model("weights", 1),
               "must be a formula")

  expect_error(.get_data_from_model("abc", 1))

  # This will error in a nonstandard way; just want to ensure it fails
  f <- function() {
    frm <- y ~ x
    .get_data_from_model("weights", str2lang("frm"))
  }
  capture_warnings(expect_error(f()))

  f2 <- function() {
    frm <- 1
    .get_data_from_model("weights", str2lang("frm"))
  }
  expect_error(f2(), "unable to convert")

})

test_that("get_data_from_model finds data in lmitt object's eval env", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  damod <- lmitt(y ~ 1, design = des, data = simdata)

  lmitt_func <- function(x, ...) {
    rhs <- formula(x)[[3L]]
    eval(rhs, environment(formula(x)))
  }

  #expect_warning(lmitt_func(damod), "No call to `model.frame`")
  out <- model.matrix(formula(damod), environment(formula(damod))$data)
  expect_true(is.numeric(out))
  expect_equal(dim(out), c(nrow(simdata), 2))
  expect_true(all(out[, 2] == damod$model[, 2]))
})

test_that(paste("get_data_from_model frame recursion doesn't return",
                "data from lmitt found in global env"), {
  rlang_trace_top_env <- getOption("rlang_trace_top_env")
  on.exit(options(rlang_trace_top_env = rlang_trace_top_env))

  options(rlang_trace_top_env = environment())

  data(simdata)
  data <- simdata
  des <- rct_design(z ~ cluster(cid1, cid2), data)
  damod <- lmitt(y ~ 1, design = des, data = data)

  non_lmitt_func <- function(data = NULL) {
    if (is.null(data)) data <- .get_data_from_model("weights",
                                                    y ~ assigned(design = des))
    data$newcol <- 5
    return(data)
  }

  expect_warning(expect_warning(non_lmitt_func(),
                                "No call to `model.frame`"),
                 "Unable to detect")
  out <- suppressWarnings(non_lmitt_func())
  expect_identical(out, cbind(simdata, newcol = 5))
})
