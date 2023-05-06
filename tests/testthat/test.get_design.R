test_that(".get_design", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  mod <- lm(y ~ x, data = simdata)

  mod1 <- lm(y ~ z, data = simdata, weights = ate(des),
             offset = cov_adj(mod))
  mod2 <- lm(y ~ z, data = simdata, weights = ate(),
             offset = cov_adj(mod, design = des))
  mod3 <- lm(y ~ z, data = simdata, weights = ate(des),
             offset = cov_adj(mod, design = des))
  mod4 <- lm(y ~ assigned(), data = simdata, weights = ate(des),
             offset = cov_adj(mod))
  mod5 <- lm(y ~ assigned(), data = simdata, weights = ate(),
             offset = cov_adj(mod, design = des))
  mod6 <- lm(y ~ assigned(), data = simdata, weights = ate(des),
             offset = cov_adj(mod, design = des))
  expect_true(all(mod1$coef == mod2$coef))
  expect_true(all(mod1$coef == mod3$coef))
  expect_true(all(mod1$coef == mod4$coef))
  expect_true(all(mod1$coef == mod5$coef))
  expect_true(all(mod1$coef == mod6$coef))

  des2 <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  expect_error(lm(y ~ assigned(), data = simdata, weights = ate(des),
                  offset = cov_adj(mod, design = des2)),
               "differing `Design`")

  expect_error(lm( y ~ assigned(), data = simdata),
               "Unable to locate")

  expect_error(lm(y ~ assigned(), data = simdata, weights = ate(),
                  offset = cov_adj(mod)),
               "Inf Recusion")

  # #37 offset in formula isntead of argument
  mod7 <- lm(y ~ z + offset(cov_adj(mod)), data = simdata, weights = ate(des))
  expect_true(all(mod1$coef == mod7$coef))
})

test_that(".get_design returns NULL with NULL_on_error = TRUE", {
  on.exit(cov_adj <- flexida::cov_adj)

  cov_adj <- function(model, newdata = NULL) {
    design <- .get_design(TRUE)

    return(design)
  }

  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  mod <- lm(y ~ x, data = simdata)
  expect_true(is.null(cov_adj(mod, newdata = simdata)))
  mod1 <- lm(y ~ z, data = simdata, offset = cov_adj(mod))
  expect_true(is.null(mod1$offset))
})

test_that(".get_design finds design in expand.model.frame call", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  damod <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)

  # ensure we're taking the Design object from the model not the existing object
  # in the environment
  uoanames <- var_names(des, "u")
  des <- rct_design(z ~ cluster(cid1), data = simdata,
                    subset = simdata$cid1 %in% c(1, 3, 5))
  dat <- .expand.model.frame.DA(damod, uoanames)
  expect_true(is.data.frame(dat))
  expect_true(all(uoanames %in% colnames(dat)))
})

test_that(".get_design finds design in model.frame call", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  damod <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)

  mf1 <- stats::model.frame(damod)
  mf2 <- stats::model.frame(formula(damod), eval(damod$call$data, environment(formula(damod))))
  mf3 <- stats::model.frame(terms(damod), eval(damod$call$data, environment(formula(damod))))

  expect_equal(mf2, mf3)
  expect_equal(lapply(mf1@.Data, as.numeric),
               lapply(as.data.frame(cbind(mf2, weights = damod$weights))@.Data,
                      as.numeric))
})
