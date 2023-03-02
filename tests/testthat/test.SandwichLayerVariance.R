test_that("vcovDA errors when provided an invalid type", {
  data(simdata)
  cmod <- lm(y ~ z + x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  expect_error(vcovDA(damod, "not_a_type"))
})

test_that("vcovDA correctly dispatches", {
  data(simdata)
  cmod <- lm(y ~ z + x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  vmat <- suppressMessages(vcovDA(damod, type = "CR0"))
  expect_equal(vmat, suppressMessages(.vcovMB_CR0(damod, cluster = .make_uoa_ids(damod))))
})

test_that(paste("vcovDA produces correct calculations with valid `cluster` arugment",
                "when cluster ID's have no NA's"), {
  data(simdata)
  simdata$uid <- seq_len(nrow(simdata))
  uid <- factor(simdata$uid)
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2, uid) + block(bid), simdata)
  dmod <- lmitt(y ~ assigned(), data = simdata, design = des,
                weights = ate(des), offset = cov_adj(cmod))

  # check default clustering level is the same when specified using cluster arg
  expect_equal(suppressMessages(vcovDA(dmod)),
               suppressMessages(vcovDA(dmod, cluster = c("uid"))))
})

test_that(paste("vcovDA produces correct calculations with valid `cluster` arugment",
                "when cluster ID's have NA's (must be via column name)"), {
  data(simdata)
  df <- rbind(simdata, simdata)
  df[1:50, c("cid1", "cid2", "bid", "z")] <- NA
  cmod <- lm(y ~ x, df[1:50,])
  des <- rct_design(z ~ cluster(cid1, cid2), df[51:100,])
  dmod <- lmitt(y ~ assigned(), data = df[51:100,], design = des,
                weights = ate(des), offset = cov_adj(cmod))

  expect_warning(vmat <- suppressMessages(vcovDA(dmod, cluster = c("cid1", "cid2"))),
                 "these observations should be treated as IID")
  expect_warning(expected <- suppressMessages(vcovDA(dmod)),
                 "these observations should be treated as IID")
  expect_equal(vmat, expected)
})

test_that(".make_uoa_ids fails without cluster argument or DirectAdjusted model", {
  data(simdata)
  mod <- lm(y ~ z, data = simdata)
  expect_error(.make_uoa_ids(mod), "Cannot deduce")
})

test_that(".make_uoa_ids fails with invalid cluster argument", {
  data(simdata)
  mod <- lm(y ~ z, data = simdata)
  expect_error(.make_uoa_ids(mod, cluster = "not_uoas"),
               "columns not_uoas in ITT effect model data")

  invalid_ids <- apply(simdata[, c("cid1", "cid2")], 1,
                       function(...) paste(..., collapse = "_"))
  expect_error(.make_uoa_ids(mod, cluster = invalid_ids),
               ", 5_2 in ITT effect model data")
})

test_that(".make_uoa_ids produces expected warnings for NA uoas in C", {
  data(simdata)
  cmod_data <- data.frame("x" = rnorm(50), "y" = rnorm(50),
                          "cid1" = rep(c(1, NA), each = 25),  "cid2" = NA)

  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod <- lmitt(y ~ assigned(), data = simdata, design = des,
                offset = cov_adj(cmod))
  Q_ids <- apply(simdata[, c("cid1", "cid2")], 1,
                 function(...) paste(..., collapse = "_"))
  C_ids <- c(rep("1_NA", 25), paste0(length(unique(Q_ids)) + seq_len(25), "*"))
  expected_ids <- factor(c(Q_ids, C_ids), levels = unique(c(Q_ids, C_ids)))

  expect_warning(
    expect_warning(ids <- .make_uoa_ids(dmod), "ID's will be clustered"),
    "should be treated as IID"
  )
  expect_equal(length(ids), nrow(cmod_data) + nrow(simdata))
  expect_equal(ids, expected_ids)
})

test_that(".make_uoa_ids returns correct ID's for non-DirectAdjusted model", {
  data(simdata)

  mod <- lm(y ~ z, data = simdata)
  expected_out <- factor(
    apply(simdata[, "cid1", drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  )

  expect_equal(.make_uoa_ids(mod, cluster = "cid1"), expected_out)
})

test_that(".make_uoa_ids returns correct ID's for non-SandwichLayer offset", {
  data(simdata)

  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  mod <- lmitt(y ~ assigned(), data = simdata, design = des)

  expected_out <- factor(
    apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  )
  expect_equal(.make_uoa_ids(mod), expected_out)
})

test_that(".make_uoa_ids returns correct ID's for full overlap of C and Q", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod <- lmitt(y ~ assigned(), data = simdata, design = des, offset = cov_adj(cmod))

  expected_out <- factor(
    apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1, function(...) paste(..., collapse = "_"))
  )
  expect_equal(.make_uoa_ids(dmod), expected_out)
})

test_that(".make_uoa_ids returns correct ID's for no overlap of C and Q", {
  data(simdata)
  set.seed(300)
  cmod_data <- data.frame("y" = rnorm(50), "x" = rnorm(50), "cid1" = NA, "cid2" = NA)

  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ uoa(cid1, cid2), simdata)
  dmod <- lmitt(y ~ assigned(), data = simdata, design = des, offset = cov_adj(cmod))

  Q_uoas <- apply(simdata[, c("cid1", "cid2"), drop = FALSE], 1,
                  function(...) paste(..., collapse = "_"))
  n_Q_uoas <- length(unique(Q_uoas))
  all_uoas <- c(Q_uoas, paste0(n_Q_uoas + seq_len(nrow(cmod_data)), "*"))
  expected_out <- factor(all_uoas, levels = unique(all_uoas))

  expect_warning(ids <- .make_uoa_ids(dmod), "treated as IID")
  expect_equal(ids, expected_out)
})

test_that("variance helper functions fail without a DirectAdjusted model", {
  data(simdata)
  cmod <- lm(y ~ z, data = simdata)

  expect_error(.vcovMB_CR0(cmod), "must be a DirectAdjusted")
  expect_error(.get_a22_inverse(cmod), "must be a DirectAdjusted")
  expect_error(.get_a11_inverse(cmod), "must be a DirectAdjusted")
  expect_error(.get_a21(cmod), "must be a DirectAdjusted")
})

test_that(paste(".get_a11_inverse, .get_a21 used with DirectAdjusted model",
                "without SandwichLayer offset"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  offset <- stats::predict(cmod, simdata)
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = offset)
  )

  expect_error(.get_a11_inverse(m), "must have an offset of class")
  expect_error(.get_a21(m), "must have an offset of class")
})

test_that(".get_a22_inverse returns correct value for lm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))

  mm <- stats::model.matrix(m)
  fim <- solve(crossprod(mm * m$weights, mm))

  expect_equal(.get_a22_inverse(m), fim)
})

test_that(".get_a22_inverse returns correct value for glm fit with Gaussian family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(glm(y ~ assigned(), data = simdata, weights = ate(des)))

  mm <- stats::model.matrix(m)
  fim <- solve(crossprod(mm * sqrt(m$weights)))

  expect_equal(.get_a22_inverse(m), fim)
})

test_that(".get_a22_inverse returns correct value for glm fit with poisson family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des),
        family = stats::poisson())
  )

  mm <- stats::model.matrix(m)
  fim <- solve(crossprod(mm * sqrt(m$weights)))

  expect_equal(.get_a22_inverse(m), fim)
})

test_that(".get_a22_inverse returns correct value for glm fit with quasipoisson family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des),
        family = stats::quasipoisson())
  )

  dispersion <- sum((m$weights * m$residuals)^2) / sum(m$weights)
  mm <- stats::model.matrix(m)
  fim <- solve(crossprod(mm * sqrt(m$weights)))

  expect_equal(.get_a22_inverse(m), fim)
})

test_that(".get_a11_inverse returns correct value for lm cmod", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  fim <- crossprod(stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), solve(fim))
})

test_that(".get_a11_inverse returns correct value for Gaussian glm cmod", {
  data(simdata)
  cmod <- glm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
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
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod),
        family = stats::poisson())
  )

  fim <- crossprod(stats::model.matrix(cmod) * exp(cmod$linear.predictors),
                   stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), solve(fim),
               tolerance = 1e-4) # tol due to diffs from chol2inv vs. solve(crossprod)
})

test_that(".get_a11_inverse returns correct value for quasipoisson glm cmod", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::quasipoisson())
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod),
        family = stats::poisson())
  )

  dispersion <- sum((cmod$weights * cmod$residuals)^2) / sum(cmod$weights)
  fim <- crossprod(stats::model.matrix(cmod) * exp(cmod$linear.predictors),
                   stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), dispersion * solve(fim),
               tolerance = 1e-4) # tol due to diffs from chol2inv vs. solve(crossprod)
})

test_that(".get_a11_inverse returns correct value for poisson glm cmod", {
  data(simdata)
  cmod <- suppressWarnings(
    glm(round(exp(y) / (1 + exp(y))) ~ x, data = simdata, family = stats::binomial())
  )
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- suppressWarnings(as.lmitt(
    glm(round(exp(y) / (1 + exp(y))) ~ assigned(), data = simdata, weights = ate(des),
        offset = cov_adj(cmod), family = stats::binomial())
  ))

  fim <- crossprod(stats::model.matrix(cmod) * cmod$fitted.values * (1 - cmod$fitted.values),
                   stats::model.matrix(cmod))
  expect_equal(.get_a11_inverse(m), solve(fim), tolerance = 1e-6)
})

test_that(".get_a21 returns correct matrix for lm cmod and lm damod w/ clustering", {
  data(simdata)

  cmod <- lm(y ~ x + force, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  m <- as.lmitt(lm(y ~ assigned(), simdata, weights = ate(des), offset = cov_adj(cmod)))

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
  m <- as.lmitt(lm(y ~ assigned(), simdata, weights = ate(des), offset = cov_adj(cmod)))

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
    glm(round(exp(y) / (1 + exp(y))) ~ assigned(), simdata, weights = ate(des),
        offset = cov_adj(cmod), family = stats::binomial())
  ))

  uoanames <- var_names(des, "u")
  uoas <- factor(Reduce(function(...) paste(..., sep = "_"),
                        stats::expand.model.frame(m, uoanames)[, uoanames]))

  Qmat <- m$prior.weights * m$family$mu.eta(m$linear.predictors) * stats::model.matrix(m)
  Cmat <- stats::model.matrix(cmod)

  a21 <- .get_a21(m)
  expect_equal(dim(a21), c(2, 3))
  expect_equal(a21, crossprod(Qmat, Cmat))
})

test_that(paste(".get_a21 returns correct matrix when data input for lmitt",
                "has NA values for some covariate values"), {
  data(simdata)
  simdata[simdata$cid1 == 1, "x"] <- NA_real_
  cmod_data <- subset(simdata, cid1 %in% c(2, 4))
  m_data <- subset(simdata, cid1 %in% c(1, 3, 5))
  cmod <- lm(y ~ x, data = cmod_data)
  des <- rct_design(z ~ cluster(cid1, cid2), m_data)

  expect_warning(expect_warning(
    lmitt(lm(y ~ assigned(), m_data, offset = cov_adj(cmod, design = des))),
    "covariance adjustments are NA"), "covariance adjustments are NA")
  # warning procs twice
  m <- suppressWarnings(
    lmitt(lm(y ~ assigned(), m_data, weights = ate(des),
             offset = cov_adj(cmod, design = des)))
  )

  damod_mm <- suppressWarnings(
    m$weights * stats::model.matrix(
      formula(m), stats::model.frame(m, na.action = na.pass))[!is.na(m_data$x),]
  )
  pg <- stats::model.matrix(formula(cmod), m_data)

  a21 <- suppressWarnings(.get_a21(m))
  expect_equal(a21, crossprod(damod_mm, pg))
})

test_that(paste(".get_a21 returns correct matrix when data input for lmitt has
                NA values for some treatment assignments"), {
  data(simdata)
  simdata[simdata$cid1 %in% c(2, 4), "z"] <- NA_integer_
  cmod <- lm(y ~ x, data = simdata, subset = cid1 %in% c(2, 4))
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  m <- lmitt(lm(y ~ assigned(), simdata, weights = ate(des), offset = cov_adj(cmod)))

  damod_mm <- m$weights * stats::model.matrix(m)
  pg <- stats::model.matrix(formula(cmod), simdata)[!is.na(simdata$z),]

  a21 <- .get_a21(m)
  expect_equal(a21, crossprod(damod_mm, pg))
})

test_that(".vcovMB_CR0 returns px2 matrix", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), simdata, weights = ate(des), offset = cov_adj(cmod)))

  vmat <- .vcovMB_CR0(m)
  expect_equal(dim(vmat), c(2, 2))
})

test_that(".vcovMB_CR0 doesn't accept `type` argument", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), simdata, weights = ate(des), offset = cov_adj(cmod)))

  expect_error(.vcovMB_CR0(m, type = "CR0"), "Cannot override the `type`")
})

test_that(".vcovMB_CR0 returns DA model sandwich if it has no SandwichLayer", {
  data(simdata)

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  m <- lmitt(y ~ assigned(), data = simdata, design = des, weights = ate(des))

  uoas <- apply(simdata[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_"))
  expect_equal(vcovDA(m, type = "CR0"),
               sandwich::sandwich(m, meat. = sandwich::meatCL, cluster = uoas))
})

test_that(paste("HC0 .vcovMB_CR0 lm w/o clustering",
                "balanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  set.seed(50)
  N <- 60

  # trt variable
  z <- rep(c(0, 1), N / 2)

  # p and mu for one categorical and one continuous cov
  px1 <- round(stats::runif(1, 0.05, 0.94), 1)
  mux2 <- round(stats::runif(1, 55, 90))
  x1 <- rep(c(0, 1), each = N / 2)
  v <- stats::rgamma(N, 50) # simulate dirichlet to balance continuous cov
  x2 <- c(t(
    Reduce(cbind, by(v, z, function(vc) vc / sum(vc) * mux2 * N / 2))
  ))

  df1 <- data.frame("z" = z, "x1" = x1, "x2" = x2, "uid" = seq_len(N))

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df1$y <- stats::model.matrix(~ x1 + x2 + z, df1) %*% theta + rnorm(N)
  df2 <- df1
  df2$uid <- NA_integer_
  df <- rbind(df1, df2)

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()
  cmod <- lm(cmod_form, df, subset = is.na(uid))
  des <- rct_design(z ~ uoa(uid), df, subset = !is.na(df$uid))
  damod <- as.lmitt(
    lm(damod_form, data = df, subset = !is.na(uid), weights = ate(des),
       offset = cov_adj(cmod))
  )
  onemod <- lm(y ~ x1 + x2 + z, data = df, subset = !is.na(uid),
               weights = ate(des))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(damod)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$uid),])
  nc <- nrow(Xstar)
  n <- nc + nrow(Z)

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(damod), solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(damod), crossprod(Z, X * damod$weights))
  expect_equal(a22inv <- sandwich::bread(damod),
               n * solve(crossprod(Z * damod$weights, Z)))

  ef_damod <- utils::getS3method("estfun", "lm")(damod)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = nc, ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat <- crossprod(estfun(damod)) / n,
               (crossprod(ef_damod) -
                  crossprod(ef_damod, ef_cmod) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(ef_cmod, ef_damod) +
                  a21 %*% a11inv %*% crossprod(ef_cmod) %*% a11inv %*% t(a21)
               ) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat, sandwich::meat(damod, adjust = FALSE))

  # with perfect group balance, a22inv %*% a21 ((Z'WZ)^(-1)Z'WX) should be 0
  # for the trt row
  expect_equal((a22inv %*% a21)[2,], setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_equal(vmat <- .vcovMB_CR0(damod, cluster = seq_len(n), cadjust = FALSE),
               (1/n) * a22inv %*% meat %*% a22inv)

  # vmat should be equal to the outputs of sandwich
  expect_equal(vmat, sandwich::sandwich(damod))

  # var_hat(z) should be smaller than var_hat(z) from onemod
  expect_true(all(diag(vmat) <
                  diag(vcov(onemod))[c(1, 4)]))
})

test_that(paste("HC0 .vcovMB_CR0 lm w/o clustering",
                "imbalanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  set.seed(50)
  N <- 60

  # trt variable
  z <- rep(c(0, 1), N / 2)

  # p and mu for one categorical and one continuous cov
  px1 <- round(stats::runif(1, 0.05, 0.94), 1)
  mux2 <- round(stats::runif(1, 55, 90))
  x1 <- rbinom(N, 1, px1)
  x2 <- stats::rgamma(N, mux2)

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2, "uid" = seq_len(N))

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + rnorm(N)
  df$uid[seq_len(3 * N / 4)] <- NA_integer_

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()
  cmod <- lm(cmod_form, df, subset = is.na(uid))
  des <- rct_design(z ~ uoa(uid), df, subset = !is.na(df$uid))
  damod <- as.lmitt(
    lm(damod_form, data = df, subset = !is.na(uid), weights = ate(des),
       offset = cov_adj(cmod))
  )
  onemod <- lm(y ~ x1 + x2 + z, data = df, subset = !is.na(uid),
               weights = ate(des))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(damod)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$uid),])
  nc <- nrow(Xstar)
  n <- nc + nrow(Z)

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(damod), solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(damod), crossprod(Z, X * damod$weights))
  expect_equal(a22inv <- sandwich::bread(damod),
               n * solve(crossprod(Z * damod$weights, Z)))

  ef_damod <- utils::getS3method("estfun", "lm")(damod)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = nc, ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat <- crossprod(estfun(damod)) / n,
               (crossprod(ef_damod) -
                  crossprod(ef_damod, ef_cmod) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(ef_cmod, ef_damod) +
                  a21 %*% a11inv %*% crossprod(ef_cmod) %*% a11inv %*% t(a21)
               ) / n)


  # meat should be the same as the output of sandwich::meat
  expect_equal(meat, sandwich::meat(damod, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_equal(vmat <- .vcovMB_CR0(damod, cluster = seq_len(n), cadjust = FALSE),
               (1/n) * a22inv %*% meat %*% a22inv)

  # vmat should be equal the outputs of sandwich
  expect_equal(vmat, sandwich::sandwich(damod))

  # var_hat(z) should be smaller than var_hat(z) from onemod
  expect_true(all(diag(vmat) <
                    diag(vcov(onemod))[c(1, 4)]))
})

test_that(paste("HC0 .vcovMB_CR0 lm w/ clustering",
                "balanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  set.seed(50)
  NC <- NQ <- 4
  MI <- 16

  # trt variable
  z <- c(rep(rep(c(0, 1), each = MI), NC / 2),
         rep(rep(c(0, 1), each = MI), NQ / 2))

  # p and mu for each matched pair's categorical and continuous cov
  px1 <- round(stats::runif((NC + NQ) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((NC + NQ) / 2, 55, 90))
  x1 <- c(
    mapply(function(x) {
      rep(c(rep(0, round(MI * x)), rep(1, round(MI * (1 - x)))), 2)
    }, px1)
  )
  x2 <- c(
    Reduce(
      cbind,
      Map(function(x) {
        # simulate dirichlet to balance continuous cov
        v <- stats::rgamma(MI * 2, 50)
        Reduce(cbind, by(v, rep(seq_len(2), each = MI),
                         function(vc) vc / sum(vc) * x * MI))
      }, mux2)
    )
  )

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2)

  # generate clustered errors
  error_sd <- round(stats::runif(NC + NQ, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(NC + NQ) - 1)) {
    msk <- (1 + i * MI):((i + 1) * MI)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = MI) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps
  df$cid <- c(rep(NA_integer_, NC * MI), rep(seq_len(NQ), each = MI))

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()

  cmod <- lm(cmod_form, df, subset = is.na(cid))
  des <- rct_design(z ~ cluster(cid), df, subset = !is.na(df$cid))
  damod <- as.lmitt(
    lm(damod_form, data = df, subset = !is.na(cid), weights = ate(des),
       offset = cov_adj(cmod))
  )

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(damod)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$cid),])
  nc <- nrow(Xstar)
  n <- nc + nrow(Z)


  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(damod), solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(damod), crossprod(Z, X * damod$weights))
  expect_equal(a22inv <- sandwich::bread(damod),
               n * solve(crossprod(Z * damod$weights, Z)))

  ids <- c(df$cid[!is.na(df$cid)],
           paste0(length(df$cid[!is.na(df$cid)]) + seq_len(length(df$cid[is.na(df$cid)])),
                  "*"))
  ef_damod <- utils::getS3method("estfun", "lm")(damod)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = nc, ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat <- crossprod(Reduce(rbind, by(estfun(damod), ids, colSums))) / n,
               (crossprod(Reduce(rbind, by(ef_damod, ids, colSums))) -
                  crossprod(Reduce(rbind, by(ef_damod, ids, colSums)),
                            Reduce(rbind, by(ef_cmod, ids, colSums))) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(ef_cmod, ids, colSums)),
                                               Reduce(rbind, by(ef_damod, ids, colSums))) +
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(ef_cmod, ids, colSums))) %*%
                  a11inv %*% t(a21)
               ) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat, sandwich::meatCL(damod, cluster = ids, cadjust = FALSE))

  # with perfect group balance, a22inv %*% a21 should have 0's in the trt row
  expect_equal((a22inv %*% a21)[2,], setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_equal(vmat <- .vcovMB_CR0(damod, cluster = ids, cadjust = FALSE),
               (1 / n) * a22inv %*% meat %*% a22inv)

  # vmat should be equal the outputs of sandwich
  expect_equal(
    vmat,
    sandwich::sandwich(damod, meat. = sandwich::meatCL, cluster = ids, cadjust = FALSE)
  )
})

test_that(paste("HC0 .vcovMB_CR0 lm w/ clustering",
                "imbalanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  set.seed(50)
  NC <- NQ <- 4
  MI <- 16

  # trt variable
  z <- c(rep(rep(c(0, 1), each = MI), NC / 2),
         rep(rep(c(0, 1), each = MI), NQ / 2))

  # p and mu for each matched pair's categorical and continuous cov
  px1 <- round(stats::runif((NC + NQ) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((NC + NQ) / 2, 55, 90))
  x1 <- c(mapply(function(x) c(rbinom(MI, 1, x), rbinom(MI, 1, x)), px1))
  x2 <- c(Reduce(cbind, Map(function(x) stats::rgamma(MI * 2, x), mux2)))

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2)

  # generate clustered errors
  error_sd <- round(stats::runif(NC + NQ, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(NC + NQ) - 1)) {
    msk <- (1 + i * MI):((i + 1) * MI)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = MI) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps
  df$cid <- c(rep(NA_integer_, NC * MI), rep(seq_len(NQ), each = MI))

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()

  cmod <- lm(cmod_form, df, subset = is.na(cid))
  des <- rct_design(z ~ cluster(cid), df, subset = !is.na(df$cid))
  damod <- as.lmitt(
    lm(damod_form, data = df, subset = !is.na(cid), weights = ate(des),
       offset = cov_adj(cmod))
  )

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(damod)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$cid),])
  nc <- nrow(Xstar)
  n <- nc + nrow(Z)


  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(damod), solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(damod), crossprod(Z, X * damod$weights))
  expect_equal(a22inv <- sandwich::bread(damod),
               n * solve(crossprod(Z * damod$weights, Z)))

  ids <- c(df$cid[!is.na(df$cid)],
           paste0(length(df$cid[!is.na(df$cid)]) + seq_len(length(df$cid[is.na(df$cid)])),
                  "*"))
  ef_damod <- utils::getS3method("estfun", "lm")(damod)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = nc, ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat <- crossprod(Reduce(rbind, by(estfun(damod), ids, colSums))) / n,
               (crossprod(Reduce(rbind, by(ef_damod, ids, colSums))) -
                  crossprod(Reduce(rbind, by(ef_damod, ids, colSums)),
                            Reduce(rbind, by(ef_cmod, ids, colSums))) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(ef_cmod, ids, colSums)),
                                               Reduce(rbind, by(ef_damod, ids, colSums))) +
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(ef_cmod, ids, colSums))) %*%
                  a11inv %*% t(a21)
               ) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat, sandwich::meatCL(damod, cluster = ids, cadjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_equal(vmat <- .vcovMB_CR0(damod, cluster = ids, cadjust = FALSE),
               (1 / n) * a22inv %*% meat %*% a22inv)

  # vmat should be the same as the outputs of sandwich
  expect_equal(
    vmat,
    sandwich::sandwich(damod, meat. = sandwich::meatCL, cluster = ids, cadjust = FALSE)
  )
})

test_that(paste("HC0 .vcovMB_CR0 lm w/o clustering",
                "imbalanced grps",
                "cmod is a strict subset of damod data", sep = ", "), {
  set.seed(50)
  N <- 60

  # trt variable
  z <- rep(c(0, 1), N / 2)

  # p and mu for one categorical and one continuous cov
  px1 <- round(stats::runif(1, 0.05, 0.94), 1)
  mux2 <- round(stats::runif(1, 55, 90))
  x1 <- rbinom(N, 1, px1)
  x2 <- stats::rgamma(N, mux2)

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2, "uid" = seq_len(N))

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + rnorm(N)

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()
  cmod_idx <- df$z == 0
  cmod <- lm(cmod_form, df[cmod_idx,])
  des <- rct_design(z ~ uoa(uid), df)
  damod <- as.lmitt(
    lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
  )

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(damod)
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df)
  nc <- nrow(Xstar)
  n <- nrow(Z)

  expect_equal(a11inv <- .get_a11_inverse(damod), solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(damod), crossprod(Z, X * damod$weights))
  expect_equal(a22inv <- sandwich::bread(damod),
               n * solve(crossprod(Z * damod$weights, Z)))

  ef_damod <- utils::getS3method("estfun", "lm")(damod)
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = n - nc, ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat <- crossprod(estfun(damod)) / n,
               (crossprod(ef_damod) -
                  crossprod(ef_damod, ef_cmod) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(ef_cmod, ef_damod) +
                  a21 %*% a11inv %*% crossprod(ef_cmod) %*% a11inv %*% t(a21)
               ) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat, sandwich::meat(damod, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_equal(vmat <- .vcovMB_CR0(damod, cluster = df$uid, cadjust = FALSE),
               (1 / n) * a22inv %*% meat %*% a22inv)

  # vmat should be the same as the outputs of sandwich
  expect_equal(vmat, sandwich::sandwich(damod, adjust = FALSE))
})

test_that(paste("HC0 .vcovMB_CR0 lm w/ clustering",
                "imbalanced grps",
                "cmod is a strict subset of damod data", sep = ", "), {
  set.seed(50)
  NC <- NQ <- 4
  MI <- 16

  # trt variable
  z <- c(rep(rep(c(0, 1), each = MI), NC / 2),
         rep(rep(c(0, 1), each = MI), NQ / 2))

  # p and mu for each matched pair's categorical and continuous cov
  px1 <- round(stats::runif((NC + NQ) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((NC + NQ) / 2, 55, 90))
  x1 <- c(mapply(function(x) c(rbinom(MI, 1, x), rbinom(MI, 1, x)), px1))
  x2 <- c(Reduce(cbind, Map(function(x) stats::rgamma(MI * 2, x), mux2)))

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2)

  # generate clustered errors
  error_sd <- round(stats::runif(NC + NQ, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(NC + NQ) - 1)) {
    msk <- (1 + i * MI):((i + 1) * MI)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = MI) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps
  df$cid <- rep(seq_len(NC + NQ), each = MI)

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()
  cmod_idx <- df$z == 0

  cmod <- lm(cmod_form, df[cmod_idx,])
  des <- rct_design(z ~ cluster(cid), df)
  damod <- as.lmitt(
    lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
  )

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(damod)
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df)
  nc <- nrow(Xstar)
  n <- nrow(Z)

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(damod), solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(damod), crossprod(Z, X * damod$weights))
  expect_equal(a22inv <- sandwich::bread(damod), n * solve(crossprod(Z * damod$weights, Z)))

  ids <- df$cid
  ef_damod <- utils::getS3method("estfun", "lm")(damod)
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = n - nc, ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat <- crossprod(Reduce(rbind, by(estfun(damod), ids, colSums))) / n,
               (crossprod(Reduce(rbind, by(ef_damod, ids, colSums))) -
                  crossprod(Reduce(rbind, by(ef_damod, ids, colSums)),
                            Reduce(rbind, by(ef_cmod, ids, colSums))) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(ef_cmod, ids, colSums)),
                                               Reduce(rbind, by(ef_damod, ids, colSums))) +
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(ef_cmod, ids, colSums))) %*%
                  a11inv %*% t(a21)
               ) / n)

  # meat should be the same as the output of sandwich::meatCL
  expect_equal(meat, sandwich::meatCL(damod, cluster = ids, cadjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_equal(vmat <- .vcovMB_CR0(damod, cluster = ids, cadjust = FALSE),
               (1 / n) * a22inv %*% meat %*% a22inv)

  # vmat should be the same as the outputs of sandwich
  expect_equal(
    vmat,
    sandwich::sandwich(damod, meat. = sandwich::meatCL, cluster = ids, cadjust = FALSE)
  )
})

test_that(paste("HC0 .vcovMB_CR0 binomial glm cmod",
                "w/o clustering",
                "imbalanced grps",
                "cmod is a strict subset of damod data", sep = ", "), {
  set.seed(50)
  N <- 60

  # trt variable
  z <- rep(c(0, 1), N / 2)

  # p and mu for one categorical and one continuous cov
  px1 <- round(stats::runif(1, 0.05, 0.94), 1)
  mux2 <- round(stats::runif(1, 55, 90))
  x1 <- rbinom(N, 1, px1)
  x2 <- stats::rgamma(N, mux2)

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2, "uid" = seq_len(N))

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  py <- 1 / (1 + exp(-scale(stats::model.matrix(~ x1 + x2 + z, df),
                            scale = FALSE) %*% theta))
  df$y <- rbinom(N, 1, py)

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()
  cmod <- glm(cmod_form, data = df, subset = z == 0, family = stats::binomial())
  des <- rct_design(z ~ uoa(uid), df)
  damod <- lmitt(
    lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
  )

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(cmod)
  fit_wstar <- cmod$weights
  rstar <- cmod$residuals
  mu_etastar <- cmod$family$mu.eta(cmod$linear.predictors)
  Z <- stats::model.matrix(damod)
  X <- stats::model.matrix(formula(stats::delete.response(terms(cmod))), df)
  w <- damod$weights
  r <- damod$residuals
  mu_eta <- cmod$family$mu.eta(drop(X %*% cmod$coefficients))
  nc <- nrow(Xstar)
  n <- nrow(Z)

  expect_equal(a11inv <- .get_a11_inverse(damod), solve(crossprod(Xstar * sqrt(fit_wstar))))
  expect_equal(a21 <- .get_a21(damod), crossprod(Z, X * w * mu_eta))
  expect_equal(a22inv <- sandwich::bread(damod), n * solve(crossprod(Z * w, Z)))

  ef_damod <- utils::getS3method("estfun", "lm")(damod)
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = n - nc, ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat <- crossprod(estfun(damod)) / n,
               (crossprod(ef_damod) -
                  crossprod(ef_damod, ef_cmod) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(ef_cmod, ef_damod) +
                  a21 %*% a11inv %*% crossprod(ef_cmod) %*% a11inv %*% t(a21)
               ) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat, sandwich::meat(damod, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_equal(vmat <- .vcovMB_CR0(damod, cluster = seq_len(n), cadjust = FALSE),
               (1 / n) * a22inv %*% meat %*% a22inv)

  # vmat should be the same as the outputs of sandwich
  expect_equal(vmat, sandwich::sandwich(damod, adjust = FALSE))
})

test_that(paste("HC0 .vcovMB_CR0 binomial glm cmod",
                "w/ clustering",
                "imbalanced grps",
                "damod is a strict subset of cmod data", sep = ", "), {
  set.seed(50)
  NC <- NQ <- 4
  MI <- 16

  # trt variable
  z <- c(rep(rep(0, MI * 2), NC / 2),
         rep(rep(c(0, 1), each = MI), NQ / 2))

  # p and mu for each matched pair's categorical and continuous cov
  px1 <- round(stats::runif((NC + NQ) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((NC + NQ) / 2, 55, 90))
  x1 <- c(mapply(function(x) c(rbinom(MI, 1, x), rbinom(MI, 1, x)), px1))
  x2 <- c(Reduce(cbind, Map(function(x) stats::rgamma(MI * 2, x), mux2)))

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2)

  # generate clustered errors
  error_sd <- round(stats::runif(NC + NQ, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(NC + NQ) - 1)) {
    msk <- (1 + i * MI):((i + 1) * MI)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = MI) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  py <- 1 / (1 + exp(-scale(stats::model.matrix(~ x1 + x2 + z, df),
                            scale = FALSE) %*% theta))
  df$y <- rbinom((NC + NQ) * MI, 1, py)
  df$cid <- rep(seq_len(NC + NQ), each = MI)

  # set trt to NA for cmod rows
  df[1:(NC * MI), "z"] <- NA_integer_

  # model
  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()
  cmod <- glm(cmod_form, data = df, family = stats::binomial())
  des <- rct_design(z ~ uoa(cid), df[!is.na(df$z),])
  damod <- lmitt(
    lm(damod_form, data = df[!is.na(df$z),], weights = ate(des), offset = cov_adj(cmod))
  )

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(cmod)
  fit_wstar <- cmod$weights
  rstar <- cmod$residuals
  mu_etastar <- cmod$family$mu.eta(cmod$linear.predictors)
  Z <- stats::model.matrix(damod)
  X <- stats::model.matrix(formula(stats::delete.response(terms(cmod))), df[!is.na(df$z),])
  w <- damod$weights
  r <- damod$residuals
  mu_eta <- cmod$family$mu.eta(drop(X %*% cmod$coefficients))
  nc <- nrow(Xstar)
  n <- nrow(Xstar)

  expect_equal(a11inv <- .get_a11_inverse(damod), solve(crossprod(Xstar * sqrt(fit_wstar))))
  expect_equal(a21 <- .get_a21(damod), crossprod(Z, X * w * mu_eta))
  expect_equal(a22inv <- sandwich::bread(damod), n * solve(crossprod(Z * w, Z)))

  ids <- c(df$cid[!is.na(df$z)], df$cid[is.na(df$z)])
  ef_damod <- utils::getS3method("estfun", "lm")(damod)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = n - nrow(X), ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  expect_equal(meat <- crossprod(Reduce(rbind, by(estfun(damod), ids, colSums))) / n,
               (crossprod(Reduce(rbind, by(ef_damod, ids, colSums))) -
                  crossprod(Reduce(rbind, by(ef_damod, ids, colSums)),
                            Reduce(rbind, by(ef_cmod, ids, colSums))) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(ef_cmod, ids, colSums)),
                                               Reduce(rbind, by(ef_damod, ids, colSums))) +
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(ef_cmod, ids, colSums))) %*%
                  a11inv %*% t(a21)
               ) / n)

  # meat should be the same as the output of sandwich::meatCL
  expect_equal(meat, sandwich::meatCL(damod, cluster = ids, cadjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_equal(vmat <- .vcovMB_CR0(damod, cluster = ids, cadjust = FALSE),
               (1 / n) * a22inv %*% meat %*% a22inv)

  # vmat should be the same as the outputs of sandwich
  expect_equal(
    vmat,
    sandwich::sandwich(damod, meat. = sandwich::meatCL, cluster = ids, cadjust = FALSE)
  )
})
