test_correct_b12 <- function(m) {
  uoanames <- var_names(m@Design, "u")
  zname <- var_names(m@Design, "t")

  # est eqns for lm are wts * resids * dmatrix
  cmod <- m$model$`(offset)`@fitted_covariance_model
  cmod_eqns <- Reduce(
    rbind,
    by(cmod$residuals * model.matrix(cmod),
       lapply(uoanames, function(col) m$model$`(offset)`@keys[,col]),
       colSums))

  Q <- stats::expand.model.frame(m, uoanames)
  msk <- !is.na(.merge_preserve_order(
    Q,
    merge(unique(m$model$`(offset)`@keys), m@Design@structure),
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)[paste0(zname, ".y")])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))

  expect_equal(.get_b12(m),
               t(cmod_eqns) %*% m_eqns)
}

test_that("variance helper functions fail without a DirectAdjusted model", {
  data(simdata)
  cmod <- lm(y ~ z, data = simdata)

  expect_error(vcovDA(cmod), "must be a DirectAdjusted")
  expect_error(.get_b12(cmod), "must be a DirectAdjusted")
  expect_error(.get_a22_inverse(cmod), "must be a DirectAdjusted")
  expect_error(.get_b22(cmod), "must be a DirectAdjusted")
  expect_error(.get_a22_inverse(cmod), "must be a DirectAdjusted")
  expect_error(.get_a11_inverse(cmod), "must be a DirectAdjusted")
  expect_error(.get_b11(cmod), "must be a DirectAdjusted")
  expect_error(.get_a21(cmod), "must be a DirectAdjusted")
})

test_that(paste(".get_b12, .get_a11_inverse, .get_b11, .get_a21 used with DirectAdjusted model",
                "without SandwichLayer offset"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  offset <- stats::predict(cmod, simdata)
  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = offset)
  )

  expect_error(vcovDA(m), "must have an offset of class")
  expect_error(.get_b12(m), "must have an offset of class")
  expect_error(.get_a11_inverse(m), "must have an offset of class")
  expect_error(.get_b11(m), "must have an offset of class")
  expect_error(.get_a21(m), "must have an offset of class")
})

test_that(paste(".get_b12 returns expected B_12 for individual-level",
                "experimental data identical to cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
   lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  test_correct_b12(m)
})

test_that(paste(".get_b12 returns expected B_12 for cluster-level",
                "experimental data whose rows fully overlap with cov model data"), {
  data(simdata)
  cluster_ids <- unique(simdata[, c("cid1", "cid2")])
  Q_cluster <- data.frame(Reduce(
    rbind,
    by(simdata,
       list(simdata$cid1, simdata$cid2),
       function(x) {colMeans(x[, c("cid1", "cid2", "x", "y", "z")])}),
  ), row.names = NULL)

  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_cluster)

  m <- as.lmitt(
   lm(y ~ z, data = Q_cluster, weights = ate(des), offset = cov_adj(cmod))
  )

  test_correct_b12(m)
})

test_that(paste(".get_b12 returns expected B_12 for individual-level",
                "experimental data that is a subset of cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata, subset = simdata$cid2 == 1)
  weighted_design <- ate(des, data = simdata[simdata$cid2 == 1,])
  m <- as.lmitt(
   lm(y ~ z, data = simdata[simdata$cid2 == 1,],
      weights = weighted_design, offset = cov_adj(cmod))
  )

  test_correct_b12(m)
})

test_that(paste(".get_b12 returns expected B_12 for cluster-level",
                "experimental data that is a subset of cov model data"), {
  data(simdata)
  subset_cluster_ids <- unique(simdata[simdata$cid2 == 1, c("cid1", "cid2")])
  Q_cluster_subset <- data.frame(Reduce(
    rbind,
    by(simdata,
       list(simdata$cid1, simdata$cid2),
       function(x) {colMeans(x[, c("cid1", "cid2", "x", "y", "z")])}),
  ), row.names = NULL)

  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = Q_cluster_subset)
  weighted_design <- ate(des, data = Q_cluster_subset)

  m <- as.lmitt(
    lm(y ~ z, data = Q_cluster_subset,
       weights = weighted_design, offset = cov_adj(cmod))
  )

  test_correct_b12(m)
})

test_that(paste(".get_b12 returns expected B_12 for experimental",
                "data that has no overlap with cov model data"), {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  C_no_cluster_ids <- simdata
  C_no_cluster_ids[, var_names(des, "u")] <- NA
  cmod <- lm(y ~ x, data = C_no_cluster_ids)

  m <- as.lmitt(
    lm(y ~ z + force, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  expect_equal(dim(.get_b12(m)), c(2, 3))
  expect_true(all(.get_b12(m) == 0))
})

test_that(paste(".get_b12 returns 0 when there's only one cluster overlapping",
                "between the covariance and direct adjustment model samples"), {
  data(simdata)
  cmod <- lm(y ~ x, simdata, subset = cid1 == 1)
  des <- rct_design(z ~ cluster(cid1), simdata, subset = simdata$cid1 %in% c(1, 5))

  msk <- simdata$cid1 %in% c(1, 5)
  m <- as.lmitt(
    lm(y ~ z, simdata[msk,],
       weights = ate(des, data = simdata[msk,]),
       offset = cov_adj(cmod, newdata = simdata[msk,]))
  )

  b12 <- .get_b12(m)
  expect_equal(dim(b12), c(2, 2))
  expect_equal(b12, matrix(0, nrow = 2, ncol = 2))
})

test_that(paste(".get_b12 returns expected value for B12 when no intercept is",
                "included in the direct adjustment model"), {
  data(simdata)
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)

  m <- as.lmitt(lm(y ~ z - 1, simdata, weights = ate(des),
                   offset = cov_adj(cmod)))

  test_correct_b12(m)
})

test_that(".get_a22_inverse returns correct value for lm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(lm(y ~ z, data = simdata, weights = ate(des)))

  fim <- crossprod(stats::model.matrix(m) * m$weights, stats::model.matrix(m))
  expect_equal(.get_a22_inverse(m), solve(fim))
})

test_that(".get_a22_inverse returns correct value for glm fit with Gaussian family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(glm(y ~ z, data = simdata, weights = ate(des)))

  dispersion <- sum((m$weights * m$residuals)^2) / sum(m$weights)
  fim <- crossprod(stats::model.matrix(m) * m$weights / dispersion, stats::model.matrix(m))
  expect_equal(.get_a22_inverse(m), solve(fim))
})

test_that(".get_a22_inverse returns correct value for glm fit with poisson family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(
    glm(round(exp(y)) ~ z, data = simdata, weights = ate(des),
        family = stats::poisson())
  )

  fim <- crossprod(stats::model.matrix(m) * exp(m$linear.predictors) * m$weights,
                   stats::model.matrix(m))
  expect_equal(.get_a22_inverse(m), solve(fim))
})

test_that(".get_a22_inverse returns correct value for glm fit with quasipoisson family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(
    glm(round(exp(y)) ~ z, data = simdata, weights = ate(des),
        family = stats::quasipoisson())
  )

  dispersion <- sum((m$weights * m$residuals)^2) / sum(m$weights)
  fim <- crossprod(stats::model.matrix(m) * exp(m$linear.predictors) * m$weights / dispersion,
                   stats::model.matrix(m))
  expect_equal(.get_a22_inverse(m), solve(fim))
})

test_that(".get_b22 returns correct value for lm object w/o offset", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(lm(y ~ z, data = simdata, weights = ate(des)))
  nq <- sum(summary(m)$df[1L:2L])
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)

  expect_equal(.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
  expect_equal(.get_b22(m, type = "HC1"),
               vmat * nuoas / (nuoas - 1L) * (nq - 1L) / (nq - 2L))
})

test_that(".get_b22 returns correct value for lm object w offset", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  nq <- sum(summary(m)$df[1L:2L])
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)

  expect_equal(.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
  expect_equal(.get_b22(m, type = "HC1"),
               vmat * nuoas / (nuoas - 1L) * (nq - 1L) / (nq - 2L))
})


test_that(".get_b22 allows custom cluster argument to meatCL", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  nq <- sum(summary(m)$df[1L:2L])
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  uoas <- Reduce(function(...) paste(..., sep = "_"),
                 stats::expand.model.frame(m, uoanames)[, uoanames])
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)

  expect_equal(.get_b22(m, cluster = uoas, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
})

test_that(".get_b22 with one clustering column", {
  data(simdata)

  simdata[simdata$cid1 == 4, "z"] <- 0
  simdata[simdata$cid1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(cid1), data = simdata)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(simdata$cid1)
  nuoas <- length(levels(uoas))

  nq <- sum(summary(m)$df[1L:2L])
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoa_matrix <- stats::model.matrix(as.formula("~ -1 + as.factor(cid1)"),
                                    stats::expand.model.frame(m, "cid1")[, "cid1", drop = FALSE])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)

  expect_equal(.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
})

test_that(".get_b22 returns corrrect value for glm fit with Gaussian family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(glm(y ~ z, data = simdata, weights = ate(des)))
  nq <- sum(summary(m)$df[1L:2L])
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)

  expect_equal(.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for poisson glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    glm(round(exp(y)) ~ z, data = simdata, weights = ate(des),
        family = stats::poisson())
  )
  nq <- sum(summary(m)$df[1L:2L])
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)

  expect_equal(.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for quasipoisson glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    glm(round(exp(y)) ~ z, data = simdata, weights = ate(des),
        family = stats::quasipoisson())
  )
  nq <- sum(summary(m)$df[1L:2L])
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)

  expect_equal(.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for binomial glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- suppressWarnings(
    as.lmitt(
      glm(round(exp(y) / (1 + exp(y))) ~ z, data = simdata, weights = ate(des),
          family = stats::binomial())
  ))

  nq <- sum(summary(m)$df[1L:2L])
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)

  expect_equal(.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})


test_that(".get_a11_inverse returns correct value for lm cmod", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  fim <- crossprod(stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), solve(fim))
})

test_that(".get_a11_inverse returns correct value for Gaussian glm cmod", {
  data(simdata)
  cmod <- glm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  dispersion <- sum((cmod$weights * cmod$residuals)^2) / sum(cmod$weights)
  fim <- crossprod(stats::model.matrix(cmod), stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), dispersion * solve(fim))
})

test_that(".get_a11_inverse returns correct value for poisson glm cmod", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::poisson())
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    glm(round(exp(y)) ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod),
        family = stats::poisson())
  )

  fim <- crossprod(stats::model.matrix(cmod) * exp(cmod$linear.predictors),
                   stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), solve(fim), tolerance = 1e-4) # tol due to diffs from chol2inv vs. solve(crossprod)
})

test_that(".get_a11_inverse returns correct value for quasipoisson glm cmod", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::quasipoisson())
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    glm(round(exp(y)) ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod),
        family = stats::poisson())
  )

  dispersion <- sum((cmod$weights * cmod$residuals)^2) / sum(cmod$weights)
  fim <- crossprod(stats::model.matrix(cmod) * exp(cmod$linear.predictors),
                   stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), dispersion * solve(fim), tolerance = 1e-4) # tol due to diffs from chol2inv vs. solve(crossprod)
})

test_that(".get_a11_inverse returns correct value for poisson glm cmod", {
  data(simdata)
  cmod <- suppressWarnings(
    glm(round(exp(y) / (1 + exp(y))) ~ x, data = simdata, family = stats::binomial())
  )
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- suppressWarnings(as.lmitt(
    glm(round(exp(y) / (1 + exp(y))) ~ z, data = simdata, weights = ate(des),
        offset = cov_adj(cmod), family = stats::binomial())
  ))

  fim <- crossprod(stats::model.matrix(cmod) * cmod$fitted.values * (1 - cmod$fitted.values),
                   stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), solve(fim), tolerance = 1e-6)
})

test_that(".get_b11 returns correct B_11 for one cluster column", {
  data(simdata)

  simdata[simdata$cid1 == 4, "z"] <- 0
  simdata[simdata$cid1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(cid1), data = simdata)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(simdata$cid1)
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )

  expect_equal(.get_b11(m), expected)
})

test_that(".get_b11 accepts custom cluster argument for meatCL", {
  data(simdata)

  simdata[simdata$cid1 == 4, "z"] <- 0
  simdata[simdata$cid1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(cid1), data = simdata)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(simdata$cid1)
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )

  expect_equal(.get_b11(m, cluster = uoas), expected)
})

test_that(".get_b11 returns correct B_11 for multiple cluster columns", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("cid1", "cid2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )

  expect_equal(.get_b11(m), expected)
})

test_that(".get_b11 returns correct B_11 for lm cmod (HC1)", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("cid1", "cid2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )

  expect_equal(.get_b11(m), expected)
})

test_that(".get_b11 returns correct B_11 for glm object (HC0)", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::poisson())
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    glm(round(exp(y)) ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod),
        family = stats::poisson()))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("cid1", "cid2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L)
  )

  expect_equal(.get_b11(m), expected)
})

test_that(paste(".get_b11 returns correct B_11 for experimental data that is a",
                "subset of cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  nc <- sum(summary(cmod)$df[1L:2L])

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata, subset = simdata$cid2 == 1)
  weighted_design <- ate(des, data = simdata[simdata$cid2 == 1,])
  m <- as.lmitt(
    lm(y ~ z, data = simdata[simdata$cid2 == 1,],
       weights = weighted_design, offset = cov_adj(cmod))
  )

  # replace NA's with distinct uoa values and recalculate nuoas for small-sample adjustment
  uoas <- Reduce(function(x, y) paste(x, y, sep = "_"), m$model$`(offset)`@keys)
  uoas[grepl("NA", uoas)] <- NA_integer_
  nuoas <- length(unique(uoas))
  nas <- is.na(uoas)
  uoas[nas] <- paste0(nuoas - 1 + seq_len(sum(nas)), "*")
  uoas <- factor(uoas)
  nuoas <- length(levels(uoas))

  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )

  expect_equal(.get_b11(m), expected)
})


test_that(paste(".get_b11 returns correct B_11 for experimental data that has",
                "no overlap with cov model data"), {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  C_no_cluster_ids <- simdata
  C_no_cluster_ids[, var_names(des, "u")] <- NA
  cmod <- lm(y ~ x, data = C_no_cluster_ids)
  nc <- sum(summary(cmod)$df[1L:2L])

  m <- as.lmitt(
    lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  # replace NA's with distinct uoa values and recalculate nuoas for small-sample adjustment
  uoas <- Reduce(function(x, y) paste(x, y, sep = "_"), m$model$`(offset)`@keys)
  uoas[grepl("NA", uoas)] <- NA_integer_
  nuoas <- length(unique(uoas))
  nas <- is.na(uoas)
  uoas[nas] <- paste0(nuoas - 1 + seq_len(sum(nas)), "*")
  uoas <- factor(uoas)

  # no adjustment for cases where each row is its own cluster
  expected <- crossprod(sandwich::estfun(cmod)) * (nc - 1L) / (nc - 2L)

  expect_equal(.get_b11(m, cadjust = FALSE), expected)
})

test_that(".get_b11 returns expected B_11 when cmod fit to one cluster", {
  data(simdata)
  cmod <- lm(y ~ x, simdata, subset = cid1 == 1)
  des <- rct_design(z ~ cluster(cid1), simdata, subset = simdata$cid1 %in% c(1, 5))

  msk <- simdata$cid1 %in% c(1, 5)
  m <- as.lmitt(
    lm(y ~ z, simdata[msk,],
       weights = ate(des, data = simdata[msk,]),
       offset = cov_adj(cmod, newdata = simdata[msk,]))
  )

  expect_equal(.get_b11(m, type = "HC0", cadjust = FALSE),
               crossprod(stats::model.matrix(cmod) * cmod$residuals))
})

test_that(".get_a21 returns correct matrix for lm cmod and lm damod w/ clustering", {
  data(simdata)

  cmod <- lm(y ~ x + force, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  m <- as.lmitt(lm(y ~ z, simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoanames <- var_names(des, "u")
  uoas <- factor(Reduce(function(...) paste(..., sep = "_"),
                        stats::expand.model.frame(m, uoanames)[, uoanames]))

  Qmat <- m$weights * stats::model.matrix(m)
  Cmat <- stats::model.matrix(cmod)

  a21 <- .get_a21(m)
  expect_equal(dim(a21), c(2, 3))
  expect_equal(a21, crossprod(Qmat, Cmat))
})

test_that(".get_a21 returns correct matrix for lm cmod and lm damod w/o clustering", {
  data(simdata)

  simdata$uid <- seq_len(nrow(simdata))
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ unitid(uid), data = simdata)
  m <- as.lmitt(lm(y ~ z, simdata, weights = ate(des), offset = cov_adj(cmod)))

  Qmat <- m$weights * stats::model.matrix(m)
  Cmat <- stats::model.matrix(cmod)

  a21 <- .get_a21(m)
  expect_equal(dim(a21), c(2, 2))
  expect_equal(a21, crossprod(Qmat, Cmat))
})

test_that(".get_a21 returns correct matrix for lm cmod and glm damod w/ clustering", {
  data(simdata)

  cmod <- lm(round(exp(y) / (1 + exp(y))) ~ x + force, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  m <- suppressWarnings(as.lmitt(
    glm(round(exp(y) / (1 + exp(y))) ~ z, simdata, weights = ate(des),
        offset = cov_adj(cmod), family = stats::binomial())
  ))

  uoanames <- var_names(des, "u")
  uoas <- factor(Reduce(function(...) paste(..., sep = "_"),
                        stats::expand.model.frame(m, uoanames)[, uoanames]))

  Qmat <- m$weights * m$family$mu.eta(m$linear.predictors) * stats::model.matrix(m)
  Cmat <- stats::model.matrix(cmod)

  a21 <- .get_a21(m)
  expect_equal(dim(a21), c(2, 3))
  expect_equal(a21, crossprod(Qmat, Cmat))
})

test_that("vcovDA returns px2 matrix", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(lm(y ~ z, simdata, weights = ate(des), offset = cov_adj(cmod)))

  vmat <- vcovDA(m)
  expect_equal(dim(vmat), c(2, 2))
})

test_that(paste("HC0 vcovDA lm w/o clustering",
                "no grp imbalance",
                "no cmod/damod data overlap", sep = ", "), {
  set.seed(50)
  nc <- 60
  nq <- 16
  beta <- c(0.5, -0.25)
  tau <- -2
  error_sd <- 0.5

  cmod_df <- data.frame(
    "x1" = rep(c(0, 1, 0, 1), each = nc/4),
    "x2" = rep(c(0, 1), each = nc/2)
  )
  cmod_df$y <- as.matrix(cmod_df) %*% beta + error_sd * rnorm(nc)
  cmod_form <- ~ x1 + x2
  X <- stats::model.matrix(cmod_form, stats::model.frame(cmod_form, cmod_df))
  cmod_df$uid <- NA_integer_
  cmod <- lm(y ~ x1 + x2, cmod_df)

  damod_df <- data.frame(
    "x1" = rep(c(0, 1, 0, 1), each = nq/4),
    "x2" = rep(rep(c(0, 1), 4), each = nq/8),
    "z" = rep(c(0, 1), each = nq/2)
  )
  damod_df$y <- as.matrix(damod_df) %*% c(beta, tau) + error_sd * rnorm(nq)
  damod_form <- ~ z
  C <- stats::model.matrix(damod_form, damod_df)
  Xstar <- stats::model.matrix(cmod_form, damod_df)

  damod_df$uid <- seq_len(nq)
  des <- rct_design(z ~ unitid(uid), damod_df)
  onemod <- lm(y ~ x1 + x2 + z, data = damod_df, weights = ate(des))
  O <- stats::model.matrix(onemod)

  damod <- as.lmitt(
    lm(y ~ z, data = damod_df, weights = ate(des), offset = cov_adj(cmod))
  )

  ## check blocks for our sandwich estimator are what we expect
  a11inv <- solve(crossprod(X))
  expect_equal(a11inv, flexida:::.get_a11_inverse(damod))
  b22 <- crossprod(C, diag((damod$weights * damod$residuals)^2)) %*% C
  expect_equal(b22, flexida:::.get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(C * damod$weights, C))
  expect_equal(a22inv, flexida:::.get_a22_inverse(damod))
  a21 <- crossprod(C, Xstar * damod$weights)
  expect_equal(a21, flexida:::.get_a21(damod))
  b12 <- matrix(0, nrow = dim(X)[2], ncol = dim(C)[2])
  expect_equal(b12, flexida:::.get_b12(damod))
  b11 <- crossprod(X * cmod$residuals^2, X)
  expect_equal(b11, flexida:::.get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(nq * a22inv, sandwich::bread(damod))

  ## meat is not equivalent to sandwich::meat! since there's no overlap of cmod
  ## and damod data, cov-adjusted meat term simplifies to:
  ## b22 + t(a21) %*% a11inv %*% b11 %*% a11inv %*% a21, which is a robust estimate
  ## of the cmod vcov matrix sandwiched by lm coefs of the covariates in the
  ## experimental sample on the trt assignment
  # we can compare to sandwich::meat when there's no clustering
  expect_equal(sandwich::meat(damod, adjust = FALSE),
               sandwich::meatCL(damod, cadjust = FALSE, type = "HC0"))

  expect_equal(b22 / nq, sandwich::meat(damod, adjust = FALSE))
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  # with perfect group balance, a22inv %*% a21 will have 0's in the trt row,
  # meaning the variance of the trt effect estimate will be the same as that
  # given by the sandwich package
  expect_true(all((a22inv %*% a21)[2,] == 0))

  ## variance of intercept estimate != sandwich package estimate, but variance
  ## of tau_hat == sandwich package estimate
  # confirm matrix multiplication for vcovDA is what we expect
  expect_equal(vcovDA(damod, type = "HC0", cadjust = FALSE),
               (1 / nq) * a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  expect_true(diag(vcovDA(damod, type = "HC0", cadjust = FALSE))[1] >
              (1 / nq) * diag(sandwich::sandwich(damod))[1])
  expect_equal(diag(vcovDA(damod, type = "HC0", cadjust = FALSE))[2],
               (1 / nq) * diag(sandwich::sandwich(damod))[2])

  ## DA model should have smaller variances than onemod
  expect_true(all(diag(vcovDA(damod, type = "HC0", cadjust = FALSE)) <
                  diag(vcov(onemod))[c(1, 4)]))
})

test_that(paste("HC0 vcovDA lm w/o clustering",
                "imbalanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  set.seed(50)
  nc <- 60
  nq <- 16
  beta <- c(0.5, -0.25)
  zeta <- c(-3, 2)
  tau <- -2
  error_sd <- 0.5

  cmod_df <- data.frame(
    "x1" = rep(c(0, 1, 0, 1), each = nc/4),
    "x2" = rep(c(0, 1), each = nc/2)
  )
  cmod_df$y <- as.matrix(cmod_df) %*% beta + error_sd * rnorm(nc)
  cmod_form <- ~ x1 + x2
  X <- stats::model.matrix(cmod_form, stats::model.frame(cmod_form, cmod_df))
  cmod_df$uid <- NA_integer_
  cmod <- lm(y ~ x1 + x2, cmod_df)

  damod_df <- data.frame(
    "x1" = rep(c(0, 1, 0, 1), each = nq/4),
    "x2" = rep(rep(c(0, 1), 4), each = nq/8)
  )
  damod_df$z <- rbinom(nq, 1,
                       exp(as.matrix(damod_df) %*% zeta) /
                        (1 + exp(as.matrix(damod_df) %*% zeta)))
  damod_df$y <- as.matrix(damod_df) %*% c(beta, tau) + error_sd * rnorm(nq)
  damod_form <- ~ z
  C <- stats::model.matrix(damod_form, damod_df)
  Xstar <- stats::model.matrix(cmod_form, damod_df)

  damod_df$uid <- seq_len(nq)
  des <- rct_design(z ~ unitid(uid), damod_df)
  onemod <- lm(y ~ x1 + x2 + z, data = damod_df, weights = ate(des))
  O <- stats::model.matrix(onemod)

  damod <- as.lmitt(
    lm(y ~ z, data = damod_df, weights = ate(des), offset = cov_adj(cmod))
  )

  ## check blocks for our sandwich estimator are what we expect
  a11inv <- solve(crossprod(X))
  expect_equal(a11inv, flexida:::.get_a11_inverse(damod))
  b22 <- crossprod(C, diag((damod$weights * damod$residuals)^2)) %*% C
  expect_equal(b22, flexida:::.get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(C * damod$weights, C))
  expect_equal(a22inv, flexida:::.get_a22_inverse(damod))
  a21 <- crossprod(C, Xstar * damod$weights)
  expect_equal(a21, flexida:::.get_a21(damod))
  b12 <- matrix(0, nrow = dim(X)[2], ncol = dim(C)[2])
  expect_equal(b12, flexida:::.get_b12(damod))
  b11 <- crossprod(X * cmod$residuals^2, X)
  expect_equal(b11, flexida:::.get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(nq * a22inv, sandwich::bread(damod))

  ## meat is not equivalent to sandwich::meat! since there's no overlap of cmod
  ## and damod data, cov-adjusted meat term simplifies to:
  ## b22 + t(a21) %*% a11inv %*% b11 %*% a11inv %*% a21, which is a robust estimate
  ## of the cmod vcov matrix sandwiched by lm coefs of the covariates in the
  ## experimental sample on the trt assignment
  # we can compare to sandwich::meat when there's no clustering
  expect_equal(sandwich::meat(damod, adjust = FALSE),
               sandwich::meatCL(damod, cadjust = FALSE, type = "HC0"))

  expect_equal(b22 / nq, sandwich::meat(damod, adjust = FALSE))
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  # with imperfect group balance, a22inv %*% a21 will not be 0 for the trt row,
  # so the variance of the trt effect estimate will increase from the sandwich
  # package's estimate
  expect_false(all((a22inv %*% a21)[2,] == 0))

  ## variance of vcovDA estimates != sandwich package estimates
  # confirm matrix multiplication for vcovDA is what we expect
  expect_equal(vcovDA(damod, type = "HC0", cadjust = FALSE),
               (1 / nq) * a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  expect_true(all(diag(vcovDA(damod, type = "HC0", cadjust = FALSE)) >
                  (1 / nq) * diag(sandwich::sandwich(damod))))

  ## DA model should have smaller variances than onemod
  expect_true(all(diag(vcovDA(damod, type = "HC0", cadjust = FALSE)) <
                    diag(vcov(onemod))[c(1, 4)]))
})

test_that(paste("HC0 vcovDA lm w/ clustering",
                "balanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  ## CLUSTERED DATA SIMULATION
  set.seed(50)
  nc <- nq <- 4
  mi <- 16

  # generate one categorical and one continuous covariate
  px1 <- round(stats::runif((nc + nq) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((nc + nq) / 2, 55, 90))

  df <- data.frame(
    "cid" = c(rep(NA_integer_, nc * mi), rep(seq_len(nq), each = mi)),
    "uid" = c(rep(NA_integer_, nc * mi), rep(seq_len(mi), nq)),
    "z" = c(rep(0, nc * mi), rep(c(0, 1), each = mi, nq / 2)),
    "x1" = c(
      mapply(function(x) {
        # perfectly balance categorical covariate within matched pairs
        rep(c(rep(0, round(mi * x)), rep(1, round(mi * (1 - x)))), 2)
      }, px1)
    ),
    "x2" = c(
      Reduce(
        cbind,
        Map(function(x) {
          # simulate dirichlet to perfectly balance continuous covariate
          # within matched pairs
          v <- stats::rgamma(mi * 2, 50)
          Reduce(cbind, by(v, rep(seq_len(2), each = mi),
                           function(vc) vc / sum(vc) * x * mi))
        }, mux2)
      )
    )
  )

  # generate clustered errors 
  error_sd <- round(stats::runif(nc + nq, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(nc + nq) - 1)) {
    msk <- (1 + i * mi):((i + 1) * mi)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = mi) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps
  
  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ z
  
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df[is.na(df$cid),])
  Xstar <- stats::model.matrix(as.formula(cmod_form[-2]), df[!is.na(df$cid),])
  C <- stats::model.matrix(as.formula(damod_form[-2]), df[!is.na(df$cid),])
  
  df$y <- c(
    X %*% theta[-length(theta)] + eps[is.na(df$cid)], # cmod sample y
    Xstar %*% theta[-length(theta)] + C[,2] * theta[length(theta)] +
      eps[!is.na(df$cid)] # damod sample y
  )
  cmod <- lm(cmod_form, df[is.na(df$cid),])
  des <- rct_design(z ~ cluster(cid), df[!is.na(df$cid),])
  damod <- as.lmitt(
    lm(damod_form, data = df[!is.na(df$cid),], weights = ate(des),
       offset = cov_adj(cmod))
  )

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  a11inv <- solve(crossprod(X))
  expect_equal(a11inv, .get_a11_inverse(damod))
  b22 <- crossprod(
    Reduce(rbind,
           by(C * damod$weights * damod$residuals, df$cid[!is.na(df$cid)],
              colSums))
  )
  expect_equal(b22, .get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(C * damod$weights, C))
  expect_equal(a22inv, .get_a22_inverse(damod))
  a21 <- crossprod(C, Xstar * damod$weights)
  expect_equal(a21, .get_a21(damod))
  b12 <- matrix(0, nrow = dim(X)[2], ncol = dim(C)[2])
  expect_equal(b12, .get_b12(damod))
  b11 <- crossprod(X * cmod$residuals^2, X)
  expect_equal(b11, .get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(nq * mi * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meatCL
  expect_equal(b22 / (nq * mi),
               sandwich::meatCL(damod, cluster = factor(df$cid[!is.na(df$cid)]),
                                cadjust = FALSE))
  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  ## HEURISTIC CHECKS
  # with perfect group balance, a22inv %*% a21 should have 0's in the trt row
  expect_equal((a22inv %*% a21)[2,],
               setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check vcovDA matches manual matrix multiplication
  expect_equal(vcovDA(damod, type = "HC0", cadjust = FALSE),
               (1 / (nq * mi)) * a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  # check that, since (a22inv %*% a21)["z",] == 0, flexida's var_hat(z) is the
  # same as that given by sandwich
  expect_equal(diag(vcovDA(damod, type = "HC0", cadjust = FALSE))[2],
               (1 / (nq * mi)) * diag(sandwich::sandwich(
                                       damod,
                                       meat. = sandwich::meatCL,
                                       cluster = factor(df$cid[!is.na(df$cid)]),
                                       cadjust = FALSE))[2])
})

test_that(paste("HC0 vcovDA lm w/ clustering",
                "imbalanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  ## CLUSTERED DATA SIMULATION
  set.seed(50)
  nc <- nq <- 4
  mi <- 16
  
  # generate one categorical and one continuous covariate
  px1 <- round(stats::runif((nc + nq) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((nc + nq) / 2, 55, 90))
  
  df <- data.frame(
    "cid" = c(rep(NA_integer_, nc * mi), rep(seq_len(nq), each = mi)),
    "uid" = c(rep(NA_integer_, nc * mi), rep(seq_len(mi), nq)),
    "z" = c(rep(0, nc * mi), rep(c(0, 1), each = mi, nq / 2)),
    "x1" = c(mapply(function(x) c(rbinom(mi, 1, x), rbinom(mi, 1, x)), px1)),
    "x2" = c(Reduce(cbind, Map(function(x) stats::rgamma(mi * 2, x), mux2)))
  )

  # generate clustered errors 
  error_sd <- round(stats::runif(nc + nq, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(nc + nq) - 1)) {
    msk <- (1 + i * mi):((i + 1) * mi)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = mi) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps
  
  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ z
  
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df[is.na(df$cid),])
  Xstar <- stats::model.matrix(as.formula(cmod_form[-2]), df[!is.na(df$cid),])
  C <- stats::model.matrix(as.formula(damod_form[-2]), df[!is.na(df$cid),])
  
  df$y <- c(
    X %*% theta[-length(theta)] + eps[is.na(df$cid)], # cmod sample y
    Xstar %*% theta[-length(theta)] + C[,2] * theta[length(theta)] +
      eps[!is.na(df$cid)] # damod sample y
  )
  cmod <- lm(cmod_form, df[is.na(df$cid),])
  des <- rct_design(z ~ cluster(cid), df[!is.na(df$cid),])
  damod <- as.lmitt(
    lm(damod_form, data = df[!is.na(df$cid),], weights = ate(des),
       offset = cov_adj(cmod))
  )
  
  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  a11inv <- solve(crossprod(X))
  expect_equal(a11inv, .get_a11_inverse(damod))
  b22 <- crossprod(
    Reduce(rbind,
           by(C * damod$weights * damod$residuals, df$cid[!is.na(df$cid)],
              colSums))
  )
  expect_equal(b22, .get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(C * damod$weights, C))
  expect_equal(a22inv, .get_a22_inverse(damod))
  a21 <- crossprod(C, Xstar * damod$weights)
  expect_equal(a21, .get_a21(damod))
  b12 <- matrix(0, nrow = dim(X)[2], ncol = dim(C)[2])
  expect_equal(b12, .get_b12(damod))
  b11 <- crossprod(X * cmod$residuals^2, X)
  expect_equal(b11, .get_b11(damod, cadjust = FALSE, type = "HC0"))
  
  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(nq * mi * a22inv, sandwich::bread(damod))
  
  # after scaling, b22 is equivalent to sandwich::meatCL
  expect_equal(b22 / (nq * mi),
               sandwich::meatCL(damod, cluster = factor(df$cid[!is.na(df$cid)]),
                                cadjust = FALSE))
  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))
  
  ## HEURISTIC CHECKS
  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))
  
  # check vcovDA matches manual matrix multiplication
  expect_equal(vcovDA(damod, type = "HC0", cadjust = FALSE),
               (1 / (nq * mi)) * a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)
  
  # check that, since (a22inv %*% a21)["z",] != 0, flexida's var_hat(z) should
  # be greater than that given by sandwich
  expect_true(diag(vcovDA(damod, type = "HC0", cadjust = FALSE))[2] >
               (1 / (nq * mi)) * diag(sandwich::sandwich(
                 damod,
                 meat. = sandwich::meatCL,
                 cluster = factor(df$cid[!is.na(df$cid)]),
                 cadjust = FALSE))[2])
})

test_that(paste("HC0 vcovDA lm w/o clustering",
                "imbalanced grps",
                "cmod is a strict subset of damod data", sep = ", "), {
  ## NO CLUSTERING DATA SIMULATION
  set.seed(50)
  nq <- 60
  pz <- 0.5

  # generate one categorical and one continuous covariate
  px1 <- round(stats::runif(1, 0.05, 0.94), 1)
  mux2 <- round(stats::runif(1, 55, 90))
  
  df <- data.frame(
    "uid" = seq_len(nq),
    "z" = c(rep(0, round(nq * pz)), rep(1, nq - round(nq * pz))),
    "x1" = rbinom(nq, 1, px1),
    "x2" = stats::rgamma(nq, mux2)
  )
  
  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + rnorm(nq)
  
  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ z
  cmod_idx <- df$z == 0
  cmod <- lm(cmod_form, df[cmod_idx,])
  des <- rct_design(z ~ uoa(uid), df)
  damod <- as.lmitt(
    lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
  )

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(as.formula(cmod_form[-2]), df[cmod_idx,])
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df)
  C <- stats::model.matrix(as.formula(damod_form[-2]), df)

  a11inv <- solve(crossprod(Xstar))
  expect_equal(a11inv, flexida:::.get_a11_inverse(damod))
  b22 <- crossprod(C, diag((damod$weights * damod$residuals)^2)) %*% C
  expect_equal(b22, flexida:::.get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(C * damod$weights, C))
  expect_equal(a22inv, flexida:::.get_a22_inverse(damod))
  a21 <- crossprod(C, X * damod$weights)
  expect_equal(a21, flexida:::.get_a21(damod))
  b12 <- crossprod(
    Reduce(rbind, by(Xstar * cmod$residuals, df$uid[cmod_idx], colSums)),
    Reduce(rbind,
           by(C[cmod_idx,] * (damod$residuals * damod$weights)[cmod_idx],
              df$uid[cmod_idx],
              colSums))
  )
  expect_equal(b12, flexida:::.get_b12(damod))
  b11 <- crossprod(Xstar * cmod$residuals^2, Xstar)
  expect_equal(b11, flexida:::.get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(nq * a22inv, sandwich::bread(damod))
  
  # after scaling, b22 is equivalent to sandwich::meat
  expect_equal(sandwich::meat(damod, adjust = FALSE),
               sandwich::meatCL(damod, cadjust = FALSE, type = "HC0"))
  expect_equal(b22 / nq, sandwich::meat(damod, adjust = FALSE))
  
  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))
  
  ## HEURISTIC CHECKS
  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))
  
  # check vcovDA matches manual matrix multiplication
  expect_equal(vcovDA(damod, type = "HC0", cadjust = FALSE),
               1 / nq * a22inv %*%
                 (b22 - a21 %*% a11inv %*% b12 - t(b12) %*% a11inv %*% t(a21) +
                    (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)
  
  # check that, since (a22inv %*% a21)["z",] != 0 and b21 != 0, flexida's
  # var_hat(z) should be greater than that given by sandwich
  expect_true(diag(vcovDA(damod, type = "HC0", cadjust = FALSE))[2] >
                1 / nq * diag(sandwich::sandwich(damod, adjust = FALSE))[2])
})