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

  vmat1 <- suppressMessages(vcovDA(damod, type = "CR0"))
  expect_equal(vmat1, suppressMessages(.vcovMB_CR0(damod, cluster = .make_uoa_ids(damod))))

  vmat2 <- suppressMessages(vcovDA(damod, type = "MB0"))
  expect_true(all.equal(vmat2, vmat1, check.attributes = FALSE))
  expect_true(attr(vmat1, "type") != attr(vmat2, "type"))

  vmat3 <- suppressMessages(vcovDA(damod, type = "HC0"))
  expect_true(all.equal(vmat3, vmat1, check.attributes = FALSE))
  expect_true(attr(vmat1, "type") != attr(vmat3, "type"))
})

test_that(paste("vcovDA produces correct calculations with valid `cluster` arugment",
                "when cluster ID's have no NA's"), {
  data(simdata)
  simdata$uid <- seq_len(nrow(simdata))
  uid <- factor(simdata$uid)
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2, uid) + block(bid), simdata)
  dmod <- lmitt(y ~ 1, data = simdata, design = des,
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
  dmod <- lmitt(y ~ 1, data = df[51:100,], design = des,
                weights = ate(des), offset = cov_adj(cmod))

  expect_warning(vmat <- suppressMessages(vcovDA(dmod, cluster = c("cid1", "cid2"))),
                 "these observations should be treated as independent")
  expect_warning(expected <- suppressMessages(vcovDA(dmod)),
                 "these observations should be treated as independent")
  expect_equal(vmat, expected)
})

test_that("variance helper functions fail without a DirectAdjusted model", {
  data(simdata)
  cmod <- lm(y ~ z, data = simdata)

  expect_error(flexida:::.vcovMB_CR0(cmod), "must be a DirectAdjusted")
  expect_error(flexida:::.get_b12(cmod), "must be a DirectAdjusted")
  expect_error(flexida:::.get_a22_inverse(cmod), "must be a DirectAdjusted")
  expect_error(flexida:::.get_b22(cmod), "must be a DirectAdjusted")
  expect_error(flexida:::.get_a22_inverse(cmod), "must be a DirectAdjusted")
  expect_error(flexida:::.get_a11_inverse(cmod), "must be a DirectAdjusted")
  expect_error(flexida:::.get_b11(cmod), "must be a DirectAdjusted")
  expect_error(flexida:::.get_a21(cmod), "must be a DirectAdjusted")

})

test_that(paste(".get_b12, .get_a11_inverse, .get_b11, .get_a21 used with DirectAdjusted model",
                "without SandwichLayer offset"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  offset <- stats::predict(cmod, simdata)
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = offset)
  )

  expect_error(.get_b12(m), "must have an offset of class")
  expect_error(.get_a11_inverse(m), "must have an offset of class")
  expect_error(.get_b11(m), "must have an offset of class")
  expect_error(.get_a21(m), "must have an offset of class")
})

test_that(".get_b12 fails if ITT model isn't strictl an `lm`", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  offset <- stats::predict(cmod, simdata)
  m <- as.lmitt(
    glm(y ~ assigned(), data = simdata, weights = ate(des), offset = offset)
  )

  expect_error(.get_b12(m), "x must be an `lm` object")
})

test_that(paste(".get_b12 returns expected B_12 for individual-level",
                "experimental data identical to cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

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
  msk <- !is.na(flexida:::.merge_preserve_order(
    Q,
    merge(unique(m$model$`(offset)`@keys), m@Design@structure),
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))
  expect_message(b12 <- flexida:::.get_b12(m), paste(nrow(simdata), "rows"))
  expect_equal(b12,
               t(cmod_eqns) %*% m_eqns)
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
    lm(y ~ assigned(), data = Q_cluster, weights = ate(des), offset = cov_adj(cmod))
  )

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
  msk <- !is.na(flexida:::.merge_preserve_order(
    Q,
    merge(unique(m$model$`(offset)`@keys), m@Design@structure),
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))
  expect_message(b12 <- flexida:::.get_b12(m), paste(nrow(simdata), "rows"))
  expect_equal(b12, t(cmod_eqns) %*% m_eqns)
})

test_that(".get_b12 fails with invalid custom cluster argument", {
  data(simdata)
  cmod_data <- cbind(simdata, new_col = factor(rbinom(nrow(simdata), 1, 0.5)))
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  expect_error(flexida:::.get_b12(m, cluster = c("cid3")),
               "cid3 are missing from the covariance adjustment model dataset")
  expect_error(flexida:::.get_b12(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
  expect_error(.get_b12(m, cluster = "new_col"),
               "in the DirectAdjusted object's Design: new_col")
})

test_that(".get_b12 produces correct estimates with valid custom cluster argument", {
  data(simdata)
  simdata[simdata$cid2 == 1, "z"] <- 0
  simdata[simdata$cid2 == 2, "z"] <- 1

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid2), data = simdata)

  m <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))

  cmod_eqns <- Reduce(
    rbind,
    by(estfun(cmod), list(simdata$cid2), colSums)
  )
  dmod_eqns <- Reduce(
    rbind,
    by(estfun(as(m, "lm")), list(simdata$cid2), colSums)
  )
  expected <- crossprod(cmod_eqns, dmod_eqns)

  # default (columns specified in `cluster` argument of Design) matches expected
  expect_equal(suppressMessages(flexida:::.get_b12(m)), expected)
  expect_equal(suppressMessages(flexida:::.get_b12(m, cluster = "cid2")), expected)
})

test_that("get_b12 handles custom cluster columns with only NA's", {
  data(simdata)
  set.seed(200)

  cmod_data <- data.frame("y" = rnorm(100), "x" = rnorm(100),
                          "cid1" = NA_integer_, "cid2" = NA_integer_)
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_message(b12 <- flexida:::.get_b12(dmod, cluster = c("cid1", "cid2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))
})

test_that(paste("get_b12 handles multiple custom cluster columns where one is",
                "only NA's and the other has no NA's"), {
  # first test is where the non-NA cluster ID's are distinct
  cmod_data <- data.frame("y" = rnorm(100), "x" = rnorm(100),
                          "cid1" = rep(seq(6, 10), each = 20),
                          "cid2" = NA_integer_)
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_message(b12 <- flexida:::.get_b12(dmod, cluster = c("cid1", "cid2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))

  # second test is where the non-NA cluster ID's overlap with design
  cmod_data$cid1 <- rep(seq(1, 5), each = 20)
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_message(b12 <- flexida:::.get_b12(dmod, cluster = c("cid1", "cid2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))
})

test_that(paste(".get_b12 handles multiple custom cluster columns where both",
                "have a mix of NA's and non-NA's"), {
  data(simdata)
  cmod_data <- rbind(simdata, simdata)
  cmod_data[(nrow(simdata)+1):(2*nrow(simdata)), c("cid1", "cid2")] <- NA_integer_
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des),
                      offset = cov_adj(cmod)))

  cmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(cmod), list(cmod_data$cid1, cmod_data$cid2), colSums)
  )
  dmod_eqns <- Reduce(
    rbind,
    by(sandwich:::estfun.lm(dmod), list(simdata$cid1, simdata$cid2), colSums)
  )
  expected <- crossprod(cmod_eqns, dmod_eqns)

  expect_message(b12 <- flexida:::.get_b12(dmod, cluster = c("cid1", "cid2")),
                 paste(nrow(simdata), "rows"))
  expect_equal(b12,
               expected)
})

test_that(paste(".get_b12 returns expected B_12 for individual-level",
                "experimental data that is a subset of cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata, subset = simdata$cid2 == 1)
  weighted_design <- ate(des, data = simdata[simdata$cid2 == 1,])
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata[simdata$cid2 == 1,],
       weights = weighted_design, offset = cov_adj(cmod))
  )

  uoanames <- var_names(m@Design, "u")
  zname <- var_names(m@Design, "t")

  # est eqns for lm are wts * resids * dmatrix
  cmod <- m$model$`(offset)`@fitted_covariance_model
  cmod_eqns <- Reduce(
    rbind,
    by((cmod$residuals * model.matrix(cmod))[m$model$`(offset)`@keys$in_Q,],
       lapply(uoanames, function(col) m$model$`(offset)`@keys[m$model$`(offset)`@keys$in_Q,col]),
       colSums))

  Q <- stats::expand.model.frame(m, uoanames)
  msk <- !is.na(.merge_preserve_order(
    Q,
    merge(unique(m$model$`(offset)`@keys[, uoanames, drop = FALSE]), m@Design@structure),
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))
  expect_message(b12 <- flexida:::.get_b12(m), paste(nrow(simdata[simdata$cid2 == 1,]), "rows"))
  expect_equal(b12,
               crossprod(cmod_eqns, m_eqns))
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
    lm(y ~ assigned(), data = Q_cluster_subset,
       weights = weighted_design, offset = cov_adj(cmod))
  )

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
  msk <- !is.na(flexida:::.merge_preserve_order(
    Q,
    merge(unique(m$model$`(offset)`@keys), m@Design@structure),
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))
  expect_equal(suppressMessages(flexida:::.get_b12(m)),
               t(cmod_eqns) %*% m_eqns)
})

test_that(paste(".get_b12 returns expected B_12 for experimental",
                "data that has no overlap with cov model data"), {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  C_no_cluster_ids <- simdata
  C_no_cluster_ids[, var_names(des, "u")] <- NA
  cmod <- lm(y ~ x, data = C_no_cluster_ids)

  m <- as.lmitt(
    lm(y ~ assigned() + force, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  expect_message(b12 <- flexida:::.get_b12(m), "0 rows")
  expect_equal(dim(b12), c(2, 3))
  expect_true(all(b12 == 0))
})

test_that(paste(".get_b12 returns B_12 with correct dimensions when only one",
                "cluster overlapps between the covariance and direct adjustment",
                "samples"), {
  data(simdata)
  set.seed(200)

  cmod_data <- data.frame("y" = rnorm(10), "x" = rnorm(10),
                          "cid1" = rep(1, 10), "cid2" = rep(1, 10))
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(cid1), simdata, subset = simdata$cid1 %in% c(1, 5))

  msk <- simdata$cid1 %in% c(1, 5)
  m <- as.lmitt(
    lm(y ~ assigned(), simdata[msk,],
       weights = ate(des, data = simdata[msk,]),
       offset = cov_adj(cmod, newdata = simdata[msk,]))
  )

  expect_warning(
    expect_message(b12 <- flexida:::.get_b12(m), paste(nrow(cmod_data), "rows")),
    "numerically indistinguishable from 0"
  )
  expect_equal(dim(b12), c(2, 2))
  expect_equal(b12, matrix(0, nrow = 2, ncol = 2))
})

test_that(paste(".get_b12 returns expected matrix when some rows in cmod data",
                "have NA values for covariates"), {
  data(simdata)
  set.seed(96)

  # random rows in cmod data have NA covariate values
  rvals <- runif(n = length(which(simdata$cid1 %in% c(2, 4))))
  names(rvals) <- which(simdata$cid1 %in% c(2, 4))
  rvals <- sort(rvals)
  simdata[as.numeric(names(rvals)[seq_len(round(length(rvals) / 4))]),
          "x"] <- NA_real_

  cmod <- lm(y ~ x, simdata, subset = cid1 %in% c(2, 4))
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  expect_warning(expect_warning(lmitt(lm(y ~ assigned(), simdata,
                                         offset = cov_adj(cmod, design = des))),
                                "adjustments are NA"), "adjustments are NA")
  ## warning procs twice
  m <- suppressWarnings(
    lmitt(lm(y ~ assigned(), simdata, offset = cov_adj(cmod, design = des)))
  )

  cmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(cmod),
       list(cid1 = simdata$cid1[!is.na(simdata$x) & simdata$cid1 %in% c(2, 4)],
            cid2 = simdata$cid2[!is.na(simdata$x) & simdata$cid1 %in% c(2, 4)]),
       colSums))
  msk <- which(simdata$cid1[!is.na(simdata$x)] %in% c(2, 4))
  dmod_eqns <- Reduce(
    rbind,
    by(estfun(as(m, "lm"))[msk, , drop = FALSE],
       list(cid1 = simdata$cid1[!is.na(simdata$x)][msk],
            cid2 = simdata$cid2[!is.na(simdata$x)][msk]),
       colSums))
  b12 <- suppressMessages(flexida:::.get_b12(m))
  expect_equal(b12, crossprod(cmod_eqns, dmod_eqns))
})

test_that(paste(".get_b12 returns expected value for B12 when no intercept is",
                "included in the direct adjustment model"), {
  data(simdata)
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)

  m <- as.lmitt(lm(y ~ assigned() - 1, simdata, weights = ate(des),
                   offset = cov_adj(cmod)))

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
  msk <- !is.na(flexida:::.merge_preserve_order(
    Q,
    merge(unique(m$model$`(offset)`@keys), m@Design@structure),
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))
  expect_equal(suppressMessages(flexida:::.get_b12(m)),
               t(cmod_eqns) %*% m_eqns)
})

test_that(".get_a22_inverse returns correct value for lm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))

  mm <- stats::model.matrix(m)
  fim <- solve(crossprod(mm * m$weights, mm))

  expect_equal(flexida:::.get_a22_inverse(m), fim)
})

test_that(".get_a22_inverse returns correct value for glm fit with Gaussian family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(glm(y ~ assigned(), data = simdata, weights = ate(des)))

  mm <- stats::model.matrix(m)
  fim <- solve(crossprod(mm * sqrt(m$weights)))

  expect_equal(flexida:::.get_a22_inverse(m), fim)
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

  expect_equal(flexida:::.get_a22_inverse(m), fim)
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

  expect_equal(flexida:::.get_a22_inverse(m), fim)
})

test_that(".get_b22 returns correct value for lm object w/o offset", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))
  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(flexida:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
  expect_equal(flexida:::.get_b22(m, type = "HC1"),
               vmat * nuoas / (nuoas - 1L) * (nq - 1L) / (nq - 2L))
})

test_that(".get_b22 returns correct value for lm object w offset", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  nq <- nrow(simdata)
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(flexida:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
  expect_equal(flexida:::.get_b22(m, type = "HC1"),
               vmat * nuoas / (nuoas - 1L) * (nq - 1L) / (nq - 2L))
})

test_that(".get_b22 fails with invalid custom cluster argument", {
  data(simdata)
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  expect_error(flexida:::.get_b22(m, cluster = c("cid3")),
               "cid3 are missing from the ITT effect model dataset")
  expect_error(flexida:::.get_b22(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
})

test_that(".get_b22 produces correct estimates with valid custom cluster argument", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  nq <- nrow(simdata)
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")

  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(flexida:::.get_b22(m, cluster = uoanames, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
})

test_that(".get_b22 with one clustering column", {
  data(simdata)

  simdata[simdata$cid1 == 4, "z"] <- 0
  simdata[simdata$cid1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(cid1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(simdata$cid1)
  nuoas <- length(levels(uoas))

  nq <- nrow(simdata)
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoa_matrix <- stats::model.matrix(as.formula("~ -1 + as.factor(cid1)"),
                                    stats::expand.model.frame(m, "cid1")[, "cid1", drop = FALSE])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(flexida:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
})

test_that(".get_b22 returns corrrect value for glm fit with Gaussian family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(glm(y ~ assigned(), data = simdata, weights = ate(des)))
  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(flexida:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for poisson glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des),
        family = stats::poisson())
  )
  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(flexida:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for quasipoisson glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des),
        family = stats::quasipoisson())
  )
  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(flexida:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for binomial glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- suppressWarnings(
    as.lmitt(
      glm(round(exp(y) / (1 + exp(y))) ~ assigned(), data = simdata, weights = ate(des),
          family = stats::binomial())
    ))

  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@Design, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(flexida:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_a11_inverse returns correct value for lm cmod", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  fim <- crossprod(stats::model.matrix(cmod))
  expect_equal(flexida:::.get_a11_inverse(m), solve(fim))
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
  expect_equal(flexida:::.get_a11_inverse(m), dispersion * solve(fim))
})

test_that(".get_a11_inverse returns correct value for poisson glm cmod", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::poisson())
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des),
        offset = cov_adj(cmod), family = stats::poisson())
  )

  fim <- crossprod(stats::model.matrix(cmod) * exp(cmod$linear.predictors),
                   stats::model.matrix(cmod))
  expect_equal(flexida:::.get_a11_inverse(m), solve(fim), tolerance = 1e-4) # tol due to diffs from chol2inv vs. solve(crossprod)
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
  expect_equal(flexida:::.get_a11_inverse(m), dispersion * solve(fim),
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
  expect_equal(flexida:::.get_a11_inverse(m), solve(fim), tolerance = 1e-6)
})

test_that(".get_b11 returns correct B_11 for one cluster column", {
  data(simdata)

  simdata[simdata$cid1 == 4, "z"] <- 0
  simdata[simdata$cid1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(cid1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des),
       offset = cov_adj(cmod)))

  uoas <- factor(simdata$cid1)
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(.get_b11(m), expected)
})

test_that(".get_b11 fails with invalid custom cluster argument", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des),
       offset = cov_adj(cmod))
  )

  expect_error(flexida:::.get_b11(m, cluster = c("cid3")),
               "cid3 are missing from the covariance adjustment model dataset")
  expect_error(flexida:::.get_b11(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
})

test_that(".get_b11 produces correct estimates with valid custom cluster argument", {
  data(simdata)

  simdata[simdata$cid1 == 4, "z"] <- 0
  simdata[simdata$cid1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(cid1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(simdata$cid1)
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )

  expect_equal(flexida:::.get_b11(m, cluster = "cid1"), expected)

  # test different clustering level
  bids <- factor(simdata[, "bid"])
  nbids <- length(levels(bids))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), simdata[, "bid"], colSums))) *
      nbids / (nbids - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(flexida:::.get_b11(m, cluster = "bid"), expected)
})

test_that(".get_b11 handles NA's correctly in custom clustering columns", {
  data(simdata)
  set.seed(200)

  # check case where all clustering columns only have NA's
  cmod_data <- data.frame("y" = rnorm(100), "x" = rnorm(100),
                          "cid1" = NA_integer_, "cid2" = NA_integer_)
  cmod <- lm(y ~ x, cmod_data)
  nc <- sum(summary(cmod)$df[1L:2L])

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_warning(flexida:::.get_b11(dmod, cluster = c("cid1", "cid2")),
                 "are found to have NA's")
  expect_equal(suppressWarnings(flexida:::.get_b11(dmod, cluster = c("cid1", "cid2"),
                                         type = "HC0", cadjust = FALSE)),
               crossprod(sandwich::estfun(cmod))) # there should be no clustering

  # check case where one clustering column doesn't only have NA's
  cmod_data$cid1 <- rep(seq(6, 10), each = 20)
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_warning(flexida:::.get_b11(dmod, cluster = c("cid1", "cid2")),
                 "have NA's for some but not all")
})

test_that(".get_b11 returns correct B_11 for multiple cluster columns", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("cid1", "cid2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(flexida:::.get_b11(m), expected)
})

test_that(".get_b11 returns correct B_11 for lm cmod (HC1)", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("cid1", "cid2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(flexida:::.get_b11(m), expected)
})

test_that(".get_b11 returns correct B_11 for glm object (HC0)", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::poisson())
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod),
        family = stats::poisson()))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("cid1", "cid2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L)
  )
  expect_equal(flexida:::.get_b11(m), expected)
})

test_that(paste(".get_b11 returns correct B_11 for experimental data that is a",
                "subset of cov model data (also tests NA's in some cluster",
                "but not all cluster columns)"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  nc <- sum(summary(cmod)$df[1L:2L])

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata, subset = simdata$cid2 == 1)
  weighted_design <- ate(des, data = simdata[simdata$cid2 == 1,])
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata[simdata$cid2 == 1,],
       weights = weighted_design, offset = cov_adj(cmod))
  )

  # replace NA's with distinct uoa values and recalculate nuoas for small-sample adjustment
  uoas <- Reduce(function(x, y) paste(x, y, sep = "_"), m$model$`(offset)`@keys)
  uoas <- factor(uoas)
  nuoas <- length(levels(uoas))

  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(flexida:::.get_b11(m), expected)
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
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
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
  expect_equal(flexida:::.get_b11(m, cadjust = FALSE), expected)
})

test_that(".get_b11 returns expected B_11 when cmod fit to one cluster", {
  data(simdata)
  cmod <- lm(y ~ x, simdata, subset = cid1 == 1)
  des <- rct_design(z ~ cluster(cid1), simdata, subset = simdata$cid1 %in% c(1, 5))

  msk <- simdata$cid1 %in% c(1, 5)
  m <- as.lmitt(
    lm(y ~ assigned(), simdata[msk,],
       weights = ate(des, data = simdata[msk,]),
       offset = cov_adj(cmod, newdata = simdata[msk,]))
  )
  expect_warning(flexida:::.get_b11(m),
                 "meat matrix numerically indistinguishable")
  expect_equal(suppressWarnings(flexida:::.get_b11(m)),
               matrix(0, nrow = 2, ncol = 2))
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

  a21 <- flexida:::.get_a21(m)
  expect_equal(dim(a21), c(2, 3))
  expect_equal(a21, crossprod(Qmat, Cmat))
})

test_that(".get_a21 returns correct matrix for lm cmod and lm damod w/o clustering", {
  data(simdata)

  new_df <- simdata
  new_df$uid <- seq_len(nrow(new_df))
  cmod <- lm(y ~ x, new_df)
  des <- rct_design(z ~ unitid(uid), data = new_df)
  m <- as.lmitt(lm(y ~ assigned(), new_df, weights = ate(des), offset = cov_adj(cmod)))

  Qmat <- m$weights * stats::model.matrix(m)
  Cmat <- stats::model.matrix(cmod)

  a21 <- flexida:::.get_a21(m)
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

  a21 <- flexida:::.get_a21(m)
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

  a21 <- suppressWarnings(flexida:::.get_a21(m))
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

  a21 <- flexida:::.get_a21(m)
  expect_equal(a21, crossprod(damod_mm, pg))
})

test_that("flexida:::.vcovMB_CR0 returns px2 matrix", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), simdata, weights = ate(des), offset = cov_adj(cmod)))

  vmat <- flexida:::.vcovMB_CR0(m)
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
  m <- lmitt(y ~ 1, data = simdata, design = des, weights = ate(des))

  uoas <- apply(simdata[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_"))
  expect_true(all.equal(vcovDA(m, type = "CR0"),
                        sandwich::sandwich(m,
                                           meat. = sandwich::meatCL,
                                           cluster = uoas),
                        check.attributes = FALSE))
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
  expect_equal((a22inv %*% a21)[2,],
               setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcovMB_CR0 matches manual matrix multiplication
  vmat <- flexida:::.vcovMB_CR0(damod, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1/n) * a22inv %*% meat %*% a22inv,
                        check.attributes = FALSE))

  # vmat should be equal to the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(damod),
                        check.attributes = FALSE))

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
  vmat <- flexida:::.vcovMB_CR0(damod, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1/n) * a22inv %*% meat %*% a22inv,
                        check.attributes = FALSE))

  # vmat should be equal the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(damod),
                       check.attributes = FALSE))

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
  expect_equal(a11inv <- flexida:::.get_a11_inverse(damod), solve(crossprod(Xstar)))
  expect_equal(a21 <- flexida:::.get_a21(damod), crossprod(Z, X * damod$weights))
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
  vmat <- flexida:::.vcovMB_CR0(damod, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * a22inv %*% meat %*% a22inv,
                        check.attributes = FALSE))

  # vmat should be equal the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(damod,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
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
  vmat <- .vcovMB_CR0(damod, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * a22inv %*% meat %*% a22inv,
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(damod,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
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
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_damod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
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
  vmat <- .vcovMB_CR0(damod, cluster = df$uid, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * a22inv %*% meat %*% a22inv,
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(damod, adjust = FALSE),
                        check.attributes = FALSE))
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
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_damod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
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
  vmat <- .vcovMB_CR0(damod, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * a22inv %*% meat %*% a22inv,
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(damod,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
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
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_damod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
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
  vmat <- .vcovMB_CR0(damod, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * a22inv %*% meat %*% a22inv,
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(damod, adjust = FALSE),
                        check.attributes = FALSE))
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
  ef_cmod <- estfun(cmod)
  nonzero_ef_damod <- estfun(as(damod, "lm"))
  ef_damod <- matrix(0, nrow = nrow(ef_cmod), ncol = ncol(nonzero_ef_damod))
  colnames(ef_damod) <- colnames(nonzero_ef_damod)
  ef_damod[!is.na(df$z), ] <- nonzero_ef_damod
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
  vmat <- .vcovMB_CR0(damod, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * a22inv %*% meat %*% a22inv,
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(damod,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
})

test_that("type attribute", {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  damod <- lmitt(y ~ 1, data = simdata, weights = "ate", design = des)
  expect_identical(attr(vcov(damod), "type"), "CR0")
  expect_identical(attr(vcov(damod, type = "CR0"), "type"), "CR0")
  expect_identical(attr(vcov(damod, type = "MB0"), "type"), "MB0")
  expect_identical(attr(vcov(damod, type = "HC0"), "type"), "HC0")
})

test_that("#119 flagging vcovDA entries as NA", {
  ### factor moderator variable
  data(simdata)
  copy_simdata <- simdata
  copy_simdata$o_fac <- as.factor(copy_simdata$o)
  des <- rct_design(z ~ cluster(cid1, cid2), copy_simdata)

  #### lmitt.formula
  damod <- lmitt(y ~ o_fac, data = copy_simdata, design = des)
  expect_warning(vc <- vcov(damod), "will be returned as NA: o_fac1, o_fac3")

  # Issue is in subgroup o_fac=1, so *find that* entry in the vcov matrix
  na_dim <- which(grepl("z._o_fac1", colnames(vc)))
  expect_true(all(is.nan(vc[na_dim, ])))
  expect_true(all(is.nan(vc[, na_dim])))
  expect_true(all(!is.nan(vc[-na_dim, -na_dim])))

  #### lmitt.lm
  damod <- lmitt(lm(y ~ o_fac + o_fac:assigned(des), data = copy_simdata), design = des)
  vc <- vcov(damod)[5:7, 5:7]

  #****************************************
  ### Setting these to NA manually only for testing purposes!
  vc[1, ] <- NA
  vc[, 1] <- NA
  ### Remove these once #119 is addressed!!!!!
  #****************************************

  # Issue is in subgroup o_fac=1, so the first entry in the vcov matrix
  expect_true(all(is.na(vc[1, ])))
  expect_true(all(is.na(vc[, 1])))
  expect_true(all(!is.na(vc[-1, -1])))
  
  ### valid continuous moderator variable
  damod <- lmitt(y ~ o, data = copy_simdata, design = des)
  vc <- vcov(damod)
  expect_true(all(!is.na(vc)))
  
  ### invalid continuous moderator variable
  copy_simdata$invalid_o <- 0
  copy_simdata$invalid_o[(copy_simdata$cid1 == 2 & copy_simdata$cid2 == 2) |
                           (copy_simdata$cid1 == 2 & copy_simdata$cid2 == 1)] <- 1
  damod <- lmitt(y ~ invalid_o, data = copy_simdata, design = des)
  expect_warning(vc <- vcov(damod), "will be returned as NA: invalid_o")
  na_dim <- which(grepl("z._invalid_o", colnames(vc)))
  expect_true(all(is.na(vc[na_dim, ])))
  expect_true(all(is.na(vc[, na_dim])))
  expect_true(all(!is.na(vc[-na_dim, -na_dim])))
})

test_that(".check_df_moderator_estimates other warnings", {
  data(simdata)
  
  # fail without a DirectAdjusted model
  nodamod <- lm(y ~ x, simdata)
  vmat <- vcov(nodamod)
  cluster_ids <- apply(simdata[, c("cid1", "cid2")], 1, function(...) paste(..., collapse = "_"))
  expect_error(.check_df_moderator_estimates(vmat, nodamod, cluster_ids),
               "must be a DirectAdjusted")
  
  # invalid `data` arg
  des <- rct_design(z ~ cluster(cid1, cid2), simdata)
  damod <- lmitt(y ~ factor(o), design = des, data = simdata)
  expect_error(.check_df_moderator_estimates(vmat, damod, cluster_ids,
                                             cbind(y = simdata$y, z = simdata$z),
                                             envir = parent.frame()),
               "`data` must be a dataframe")
  
  # column names of vmat don't match moderator variable of model
  expect_warning(
    expect_warning(.check_df_moderator_estimates(vmat, damod, cluster_ids),
                   "will be returned as NA"),
    "will not be returned as NA")
})

test_that("#123 ensure PreSandwich are converted to Sandwich", {
  data(simdata)
  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  cmod <- lm(y ~ x, data = simdata)
  # Make sure its PreSandwich prior
  ca <- cov_adj(cmod, newdata = simdata)
  expect_false(is(ca, "SandwichLayer"))
  damod <- as.lmitt(lm(y ~ a.(des), data = simdata, offset = ca),
                    design = des)
  expect_true(is(damod$model$`(offset)`, "SandwichLayer"))


  copy_simdata <- simdata
  copy_simdata$schoolid <- seq(900, 949)
  C_df <- rbind(copy_simdata[, c("schoolid", "x", "y")],
                data.frame(schoolid = seq(1000, 1049),
                           x = rnorm(50),
                           y = rnorm(50)))
  cmod <- lm(y ~ x, C_df)
  des <- rct_design(z ~ uoa(schoolid), copy_simdata)
  lm1 <- lm(y ~ assigned(), copy_simdata, weights = ett(design = des),
            offset = cov_adj(cmod, NULL, NULL, "schoolid"))
  dmod1 <- lmitt(lm1, design = des)
  expect_true(is(damod$model$`(offset)`, "SandwichLayer"))
  v1 <- vcovDA(dmod1)

  wts2 <- ett(des, data = copy_simdata)
  offst2 <- cov_adj(cmod, copy_simdata, by = "schoolid")
  lm2 <- lm(y ~ assigned(), copy_simdata, weights = wts2, offset = offst2)
  dmod2 <- lmitt(lm2, design = des)
  expect_true(is(damod$model$`(offset)`, "SandwichLayer"))
  v2 <- vcovDA(dmod2)

  expect_true(all.equal(v1, v2))

})
