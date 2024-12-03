test_that(".get_spec", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata)

  mod1 <- lm(y ~ z, data = simdata, weights = ate(spec),
             offset = cov_adj(mod))
  mod2 <- lm(y ~ z, data = simdata, weights = ate(),
             offset = cov_adj(mod, specification = spec))
  mod3 <- lm(y ~ z, data = simdata, weights = ate(spec),
             offset = cov_adj(mod, specification = spec))
  mod4 <- lm(y ~ assigned(), data = simdata, weights = ate(spec),
             offset = cov_adj(mod))
  mod5 <- lm(y ~ assigned(), data = simdata, weights = ate(),
             offset = cov_adj(mod, specification = spec))
  mod6 <- lm(y ~ assigned(), data = simdata, weights = ate(spec),
             offset = cov_adj(mod, specification = spec))
  mod7 <- lm(y ~ assigned(), data = simdata, weights = ate(spec),
             offset = propertee::cov_adj(mod))
  expect_true(all(mod1$coef == mod2$coef))
  expect_true(all(mod1$coef == mod3$coef))
  expect_true(all(mod1$coef == mod4$coef))
  expect_true(all(mod1$coef == mod5$coef))
  expect_true(all(mod1$coef == mod6$coef))
  expect_true(all(mod1$coef == mod7$coef))

  spec2 <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  expect_error(lm(y ~ assigned(), data = simdata, weights = ate(spec),
                  offset = cov_adj(mod, specification = spec2)),
               "differing `StudySpecification`")

  expect_error(lm( y ~ assigned(), data = simdata),
               "Unable to locate")

  expect_error(lm(y ~ assigned(), data = simdata, weights = ate(),
                  offset = cov_adj(mod)),
               "Inf Recusion")

  # #37 offset in formula isntead of argument
  mod7 <- lm(y ~ z + offset(cov_adj(mod)), data = simdata, weights = ate(spec))
  expect_true(all(mod1$coef == mod7$coef))
})

test_that(".get_spec returns NULL with NULL_on_error = TRUE", {
  on.exit(cov_adj <- propertee::cov_adj)

  cov_adj <- function(model, newdata = NULL) {
    specification <- .get_spec(TRUE)

    return(specification)
  }

  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  mod <- lm(y ~ x, data = simdata)
  expect_true(is.null(cov_adj(mod, newdata = simdata)))
  mod1 <- lm(y ~ z, data = simdata, offset = cov_adj(mod))
  expect_true(is.null(mod1$offset))
})

test_that(".get_spec finds specification in expand.model.frame call", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  ssmod <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)

  # ensure we're taking the StudySpecification object from the model not the existing object
  # in the environment
  uoanames <- var_names(spec, "u")
  spec <- rct_spec(z ~ cluster(uoa1), data = simdata,
                    subset = simdata$uoa1 %in% c(1, 3, 5))
  dat <- .expand.model.frame_teeMod(ssmod, uoanames)
  expect_true(is.data.frame(dat))
  expect_true(all(uoanames %in% colnames(dat)))
})

test_that(".get_spec finds specification in model.frame call", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  ssmod <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)

  mf1 <- stats::model.frame(ssmod)
  mf2 <- stats::model.frame(formula(ssmod), ssmod$call$data)
  mf3 <- stats::model.frame(terms(ssmod), ssmod$call$data)

  expect_equal(mf2, mf3)
  expect_equal(lapply(mf1@.Data, as.numeric),
               lapply(as.data.frame(cbind(mf2, weights = ssmod$weights))@.Data,
                      as.numeric))
})
