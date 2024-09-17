test_that("vcov_tee errors when provided an invalid type", {
  data(simdata)
  cmod <- lm(y ~ z + x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  expect_error(vcov_tee(damod, "not_a_type"))
})

test_that("vcov_tee correctly dispatches", {
  data(simdata)
  cmod <- lm(y ~ z + x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
  damod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  vmat1 <- suppressMessages(vcov_tee(damod, type = "CR0"))
  expect_equal(vmat1, suppressMessages(.vcov_CR0(damod, cluster = .make_uoa_ids(damod, "CR"))))

  vmat2 <- suppressMessages(vcov_tee(damod, type = "MB0"))
  expect_true(all.equal(vmat2, vmat1, check.attributes = FALSE))
  expect_true(attr(vmat1, "type") != attr(vmat2, "type"))

  vmat3 <- suppressMessages(vcov_tee(damod, type = "HC0"))
  expect_true(all.equal(vmat3, vmat1, check.attributes = FALSE))
  expect_true(attr(vmat1, "type") != attr(vmat3, "type"))
})

test_that(paste("vcov_tee produces correct calculations with valid `cluster` arugment",
                "when cluster ID's have no NA's"), {
  data(simdata)
  simdata_copy <- simdata
  simdata_copy$uid <- seq_len(nrow(simdata_copy))
  cmod <- lm(y ~ x, simdata_copy)
  des <- rct_design(z ~ cluster(uoa1, uoa2, uid) + block(bid), simdata_copy)
  dmod <- lmitt(y ~ 1, data = simdata_copy, design = des,
                weights = ate(des), offset = cov_adj(cmod))

  # check default clustering level is the same when specified using cluster arg
  expect_equal(suppressMessages(vcov_tee(dmod)),
               suppressMessages(vcov_tee(dmod, cluster = c("uid"))))
})

test_that(paste("vcov_tee produces correct calculations with valid `cluster` arugment",
                "when cluster ID's have NA's (must be via column name)"), {
  data(simdata)
  df <- rbind(simdata, simdata)
  df[1:50, c("uoa1", "uoa2", "bid", "z")] <- NA
  cmod <- lm(y ~ x, df[1:50,])
  des <- rct_design(z ~ cluster(uoa1, uoa2), df[51:100,])
  dmod <- lmitt(y ~ 1, data = df[51:100,], design = des,
                weights = ate(des), offset = cov_adj(cmod))

  expect_equal(vcov_tee(dmod, cluster = c("uoa1", "uoa2")), vcov_tee(dmod))
})

test_that("variance helper functions fail without a teeMod model", {
  data(simdata)
  cmod <- lm(y ~ z, data = simdata)

  expect_error(propertee:::.vcov_CR0(cmod), "must be a teeMod")
  expect_error(propertee:::.get_b12(cmod), "must be a teeMod")
  expect_error(propertee:::.get_a22_inverse(cmod), "must be a teeMod")
  expect_error(propertee:::.get_b22(cmod), "must be a teeMod")
  expect_error(propertee:::.get_a22_inverse(cmod), "must be a teeMod")
  expect_error(propertee:::.get_a11_inverse(cmod), "must be a teeMod")
  expect_error(propertee:::.get_b11(cmod), "must be a teeMod")
  expect_error(propertee:::.get_a21(cmod), "must be a teeMod")

})

test_that(paste(".get_b12, .get_a11_inverse, .get_b11, .get_a21 used with teeMod model",
                "without SandwichLayer offset"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

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
  msk <- !is.na(propertee:::.merge_preserve_order(
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
  expect_message(b12 <- propertee:::.get_b12(m), paste(nrow(simdata), "rows"))
  expect_equal(b12,
               t(cmod_eqns) %*% m_eqns)
})

test_that(paste(".get_b12 returns expected B_12 for cluster-level",
                "experimental data whose rows fully overlap with cov model data"), {
  data(simdata)
  cluster_ids <- unique(simdata[, c("uoa1", "uoa2")])
  Q_cluster <- data.frame(Reduce(
    rbind,
    by(simdata,
       list(simdata$uoa1, simdata$uoa2),
       function(x) {colMeans(x[, c("uoa1", "uoa2", "x", "y", "z")])}),
  ), row.names = NULL)

  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = Q_cluster)

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
  msk <- !is.na(propertee:::.merge_preserve_order(
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
  expect_message(b12 <- propertee:::.get_b12(m), paste(nrow(simdata), "rows"))
  expect_equal(b12, t(cmod_eqns) %*% m_eqns)
})

test_that(".get_b12 fails with invalid custom cluster argument", {
  data(simdata)
  cmod_data <- cbind(simdata, new_col = factor(rbinom(nrow(simdata), 1, 0.5)))
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  expect_error(propertee:::.get_b12(m, cluster = c("cid3")),
               "cid3 are missing from the covariance adjustment model dataset")
  expect_error(propertee:::.get_b12(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
  expect_error(.get_b12(m, cluster = "new_col"),
               "in the teeMod object's Design: new_col")
})

test_that(".get_b12 produces correct estimates with valid custom cluster argument", {
  data(simdata)
  simdata[simdata$uoa2 == 1, "z"] <- 0
  simdata[simdata$uoa2 == 2, "z"] <- 1

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(uoa2), data = simdata)

  m <- lmitt(y ~ 1, data = simdata, design = des, offset = cov_adj(cmod))

  cmod_eqns <- Reduce(
    rbind,
    by(estfun(cmod), list(simdata$uoa2), colSums)
  )
  dmod_eqns <- Reduce(
    rbind,
    by(estfun(as(m, "lm")), list(simdata$uoa2), colSums)
  )
  expected <- crossprod(cmod_eqns, dmod_eqns)

  # default (columns specified in `cluster` argument of Design) matches expected
  expect_equal(suppressMessages(propertee:::.get_b12(m)), expected)
  expect_equal(suppressMessages(propertee:::.get_b12(m, cluster = "uoa2")), expected)
})

test_that("get_b12 handles custom cluster columns with only NA's", {
  data(simdata)
  set.seed(200)

  cmod_data <- data.frame("y" = rnorm(100), "x" = rnorm(100),
                          "uoa1" = NA_integer_, "uoa2" = NA_integer_)
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_message(b12 <- propertee:::.get_b12(dmod, cluster = c("uoa1", "uoa2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))
})

test_that(paste("get_b12 handles multiple custom cluster columns where one is",
                "only NA's and the other has no NA's"), {
  # first test is where the non-NA cluster ID's are distinct
  cmod_data <- data.frame("y" = rnorm(100), "x" = rnorm(100),
                          "uoa1" = rep(seq(6, 10), each = 20),
                          "uoa2" = NA_integer_)
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_message(b12 <- propertee:::.get_b12(dmod, cluster = c("uoa1", "uoa2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))

  # second test is where the non-NA cluster ID's overlap with design
  cmod_data$uoa1 <- rep(seq(1, 5), each = 20)
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_message(b12 <- propertee:::.get_b12(dmod, cluster = c("uoa1", "uoa2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))
})

test_that(paste(".get_b12 handles multiple custom cluster columns where both",
                "have a mix of NA's and non-NA's"), {
  data(simdata)
  cmod_data <- rbind(simdata, simdata)
  cmod_data[(nrow(simdata)+1):(2*nrow(simdata)), c("uoa1", "uoa2")] <- NA_integer_
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des),
                      offset = cov_adj(cmod)))

  cmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(cmod), list(cmod_data$uoa1, cmod_data$uoa2), colSums)
  )
  dmod_eqns <- Reduce(
    rbind,
    by(sandwich:::estfun.lm(dmod), list(simdata$uoa1, simdata$uoa2), colSums)
  )
  expected <- crossprod(cmod_eqns, dmod_eqns)

  expect_message(b12 <- propertee:::.get_b12(dmod, cluster = c("uoa1", "uoa2")),
                 paste(nrow(simdata), "rows"))
  expect_equal(b12,
               expected)
})

test_that(paste(".get_b12 returns expected B_12 for individual-level",
                "experimental data that is a subset of cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata, subset = simdata$uoa2 == 1)
  weighted_design <- ate(des, data = simdata[simdata$uoa2 == 1,])
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata[simdata$uoa2 == 1,],
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
  expect_message(b12 <- propertee:::.get_b12(m), paste(nrow(simdata[simdata$uoa2 == 1,]), "rows"))
  expect_equal(b12,
               crossprod(cmod_eqns, m_eqns))
})

test_that(paste(".get_b12 returns expected B_12 for cluster-level",
                "experimental data that is a subset of cov model data"), {
  data(simdata)
  subset_cluster_ids <- unique(simdata[simdata$uoa2 == 1, c("uoa1", "uoa2")])
  Q_cluster_subset <- data.frame(Reduce(
    rbind,
    by(simdata,
       list(simdata$uoa1, simdata$uoa2),
       function(x) {colMeans(x[, c("uoa1", "uoa2", "x", "y", "z")])}),
  ), row.names = NULL)

  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = Q_cluster_subset)
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
  msk <- !is.na(propertee:::.merge_preserve_order(
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
  expect_equal(suppressMessages(propertee:::.get_b12(m)),
               t(cmod_eqns) %*% m_eqns)
})

test_that(paste(".get_b12 returns expected B_12 for experimental",
                "data that has no overlap with cov model data"), {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  C_no_cluster_ids <- simdata
  C_no_cluster_ids[, var_names(des, "u")] <- NA
  cmod <- lm(y ~ x, data = C_no_cluster_ids)

  m <- as.lmitt(
    lm(y ~ assigned() + force, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  expect_message(b12 <- propertee:::.get_b12(m), "0 rows")
  expect_equal(dim(b12), c(2, 3))
  expect_true(all(b12 == 0))
})

test_that(paste(".get_b12 returns B_12 with correct dimensions when only one",
                "cluster overlapps between the covariance and direct adjustment",
                "samples"), {
  data(simdata)
  set.seed(200)

  cmod_data <- data.frame("y" = rnorm(10), "x" = rnorm(10),
                          "uoa1" = rep(1, 10), "uoa2" = rep(1, 10))
  cmod <- lm(y ~ x, cmod_data)
  des <- rct_design(z ~ cluster(uoa1), simdata, subset = simdata$uoa1 %in% c(1, 5))

  msk <- simdata$uoa1 %in% c(1, 5)
  m <- as.lmitt(
    lm(y ~ assigned(), simdata[msk,],
       weights = ate(des, data = simdata[msk,]),
       offset = cov_adj(cmod, newdata = simdata[msk,]))
  )

  expect_warning(
    expect_message(b12 <- propertee:::.get_b12(m), paste(nrow(cmod_data), "rows")),
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
  rvals <- runif(n = length(which(simdata$uoa1 %in% c(2, 4))))
  names(rvals) <- which(simdata$uoa1 %in% c(2, 4))
  rvals <- sort(rvals)
  simdata[as.numeric(names(rvals)[seq_len(round(length(rvals) / 4))]),
          "x"] <- NA_real_

  cmod <- lm(y ~ x, simdata, subset = uoa1 %in% c(2, 4))
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
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
       list(uoa1 = simdata$uoa1[!is.na(simdata$x) & simdata$uoa1 %in% c(2, 4)],
            uoa2 = simdata$uoa2[!is.na(simdata$x) & simdata$uoa1 %in% c(2, 4)]),
       colSums))
  msk <- which(simdata$uoa1[!is.na(simdata$x)] %in% c(2, 4))
  dmod_eqns <- Reduce(
    rbind,
    by(estfun(as(m, "lm"))[msk, , drop = FALSE],
       list(uoa1 = simdata$uoa1[!is.na(simdata$x)][msk],
            uoa2 = simdata$uoa2[!is.na(simdata$x)][msk]),
       colSums))
  b12 <- suppressMessages(propertee:::.get_b12(m))
  expect_equal(b12, crossprod(cmod_eqns, dmod_eqns))
})

test_that(paste(".get_b12 returns expected value for B12 when no intercept is",
                "included in the direct adjustment model"), {
  data(simdata)
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)

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
  msk <- !is.na(propertee:::.merge_preserve_order(
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
  expect_equal(suppressMessages(propertee:::.get_b12(m)),
               t(cmod_eqns) %*% m_eqns)
})

test_that(".get_a22_inverse correct w/o covariance adjustment", {
  data(simdata)

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  m_as.lmitt <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des, weights = ate(des))

  nq <- nrow(stats::model.frame(m_as.lmitt))
  inv_fim <- nq * chol2inv(m_as.lmitt$qr$qr)

  expect_true(all.equal(propertee:::.get_a22_inverse(m_as.lmitt), inv_fim,
                        check.attributes = FALSE))
  expect_true(all.equal(propertee:::.get_a22_inverse(m_lmitt.form), inv_fim,
                        check.attributes = FALSE))
})

test_that(".get_a22_inverse correct w/ covariance adjustment", {
  set.seed(438)
  data(simdata)
  nc <- 30
  nq <- nrow(simdata)
  n <- nc + nq
  cmod_data <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_data)

  simdata$id <- seq(nq)
  des <- rct_design(z ~ unitid(id), simdata)

  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)
  ))
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des,
                        weights = ate(des), offset = cov_adj(cmod))

  inv_fim <- nq * chol2inv(m_as.lmitt$qr$qr)

  expect_true(all.equal(propertee:::.get_a22_inverse(m_as.lmitt), inv_fim,
                        check.attributes = FALSE))
  expect_true(all.equal(propertee:::.get_a22_inverse(m_lmitt.form), inv_fim,
                        check.attributes = FALSE))
})

test_that(".get_tilde_a22_inverse correct w/o covariance adjustment", {
  data(simdata)

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  m_as.lmitt <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des, weights = ate(des))

  expect_equal(.get_tilde_a22_inverse(m_as.lmitt), .get_a22_inverse(m_as.lmitt))
  expect_equal(.get_tilde_a22_inverse(m_lmitt.form), .get_a22_inverse(m_lmitt.form))
})

test_that(".get_a22_inverse correct w/o covariance adjustment", {
  set.seed(438)
  data(simdata)
  nc <- 30
  nq <- nrow(simdata)
  n <- nc + nq
  cmod_data <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_data)

  simdata$id <- seq(nq)
  des <- rct_design(z ~ unitid(id), simdata)

  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)
  ))
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des,
                        weights = ate(des), offset = cov_adj(cmod))

  inv_fim <- n * chol2inv(m_as.lmitt$qr$qr)

  expect_true(all.equal(propertee:::.get_tilde_a22_inverse(m_as.lmitt), inv_fim,
                        check.attributes = FALSE))
  expect_true(all.equal(propertee:::.get_tilde_a22_inverse(m_lmitt.form), inv_fim,
                        check.attributes = FALSE))
})

test_that(".get_b22 returns correct value for lm object w/o offset", {
  data(simdata)

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  expect_equal(propertee:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
  expect_equal(propertee:::.get_b22(m, type = "HC1"),
               vmat * nuoas / (nuoas - 1L) * (nq - 1L) / (nq - 2L))
})

test_that(".get_b22 returns correct value for lm object w offset", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  expect_equal(propertee:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
  expect_equal(propertee:::.get_b22(m, type = "HC1"),
               vmat * nuoas / (nuoas - 1L) * (nq - 1L) / (nq - 2L))
})

test_that(".get_b22 fails with invalid custom cluster argument", {
  data(simdata)
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  expect_error(propertee:::.get_b22(m, cluster = c("cid3")),
               "cid3 are missing from the direct adjustment")
  expect_error(propertee:::.get_b22(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
})

test_that(".get_b22 produces correct estimates with valid custom cluster argument", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  expect_equal(propertee:::.get_b22(m, cluster = uoanames, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
})

test_that(".get_b22 with one clustering column", {
  data(simdata)

  simdata[simdata$uoa1 == 4, "z"] <- 0
  simdata[simdata$uoa1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(uoa1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(simdata$uoa1)
  nuoas <- length(levels(uoas))

  nq <- nrow(simdata)
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoa_matrix <- stats::model.matrix(as.formula("~ -1 + as.factor(uoa1)"),
                                    stats::expand.model.frame(m, "uoa1")[, "uoa1", drop = FALSE])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(propertee:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
})

test_that(".get_b22 returns corrrect value for glm fit with Gaussian family", {
  data(simdata)

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  expect_equal(propertee:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for poisson glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  expect_equal(propertee:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for quasipoisson glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  expect_equal(propertee:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_b22 returns correct value for binomial glm", {
  data(simdata)

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  expect_equal(propertee:::.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1))
})

test_that(".get_a11_inverse returns correct value for lm cmod", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des,
                        weights = ate(des), offset = cov_adj(cmod))

  nc <- nrow(stats::model.frame(cmod))
  fim <- crossprod(stats::model.matrix(cmod))
  expect_equal(propertee:::.get_a11_inverse(m_as.lmitt), nc * solve(fim))
  expect_equal(propertee:::.get_a11_inverse(m_lmitt.form), nc * solve(fim))
})

test_that(".get_a11_inverse returns correct value for Gaussian glm cmod", {
  data(simdata)
  cmod <- glm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des,
                        weights = ate(des), offset = cov_adj(cmod))

  nc <- nrow(stats::model.frame(cmod))
  dispersion <- sum((cmod$weights * cmod$residuals)^2) / sum(cmod$weights)
  fim <- crossprod(stats::model.matrix(cmod), stats::model.matrix(cmod))
  expect_equal(propertee:::.get_a11_inverse(m_as.lmitt), nc * dispersion * solve(fim))
  expect_equal(propertee:::.get_a11_inverse(m_lmitt.form), nc * dispersion * solve(fim))
})

test_that(".get_a11_inverse returns correct value for poisson glm cmod", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::poisson())
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des,
                        weights = ate(des), offset = cov_adj(cmod))

  nc <- nrow(stats::model.frame(cmod))
  fim <- crossprod(stats::model.matrix(cmod) * exp(cmod$linear.predictors),
                   stats::model.matrix(cmod))
  # tol due to diffs from chol2inv vs. solve(crossprod)
  expect_equal(propertee:::.get_a11_inverse(m_as.lmitt), nc * solve(fim), tolerance = 1e-4)
  expect_equal(propertee:::.get_a11_inverse(m_lmitt.form), nc * solve(fim), tolerance = 1e-4)
})

test_that(".get_a11_inverse returns correct value for quasipoisson glm cmod", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::quasipoisson())
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des,
                        weights = ate(des), offset = cov_adj(cmod))

  nc <- nrow(stats::model.frame(cmod))
  dispersion <- sum((cmod$weights * cmod$residuals)^2) / sum(cmod$weights)
  fim <- crossprod(stats::model.matrix(cmod) * exp(cmod$linear.predictors),
                   stats::model.matrix(cmod))
  # tol due to diffs from chol2inv vs. solve(crossprod)
  expect_equal(propertee:::.get_a11_inverse(m_as.lmitt), nc * dispersion * solve(fim),
               tolerance = 1e-4)
  expect_equal(propertee:::.get_a11_inverse(m_lmitt.form), nc * dispersion * solve(fim),
               tolerance = 1e-4)
})

test_that(".get_a11_inverse returns correct value for binomial glm cmod", {
  data(simdata)
  cmod <- suppressWarnings(
    glm(round(exp(y) / (1 + exp(y))) ~ x, data = simdata, family = stats::binomial())
  )
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, design = des,
                        weights = ate(des), offset = cov_adj(cmod))

  nc <- nrow(stats::model.frame(cmod))
  fim <- crossprod(stats::model.matrix(cmod) * cmod$fitted.values * (1 - cmod$fitted.values),
                   stats::model.matrix(cmod))
  # tol due to diffs from chol2inv vs. solve(crossprod)
  expect_equal(propertee:::.get_a11_inverse(m_as.lmitt), nc * solve(fim), tolerance = 1e-6)
  expect_equal(propertee:::.get_a11_inverse(m_lmitt.form), nc * solve(fim), tolerance = 1e-6)
})

test_that(".get_b11 returns correct B_11 for one cluster column", {
  data(simdata)

  simdata[simdata$uoa1 == 4, "z"] <- 0
  simdata[simdata$uoa1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(uoa1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des),
       offset = cov_adj(cmod)))

  uoas <- factor(simdata$uoa1)
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
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(des@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des),
       offset = cov_adj(cmod))
  )

  expect_error(propertee:::.get_b11(m, cluster = c("cid3")),
               "cid3 are missing from the covariance adjustment model dataset")
  expect_error(propertee:::.get_b11(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
})

test_that(".get_b11 produces correct estimates with valid custom cluster argument", {
  data(simdata)

  simdata[simdata$uoa1 == 4, "z"] <- 0
  simdata[simdata$uoa1 == 2, "z"] <- 1
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ uoa(uoa1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(simdata$uoa1)
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )

  expect_equal(propertee:::.get_b11(m, cluster = "uoa1"), expected)

  # test different clustering level
  bids <- factor(simdata[, "bid"])
  nbids <- length(levels(bids))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), simdata[, "bid"], colSums))) *
      nbids / (nbids - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(propertee:::.get_b11(m, cluster = "bid"), expected)
})

test_that(".get_b11 handles NA's correctly in custom clustering columns", {
  data(simdata)
  set.seed(200)

  # check case where all clustering columns only have NA's
  cmod_data <- data.frame("y" = rnorm(100), "x" = rnorm(100),
                          "uoa1" = NA_integer_, "uoa2" = NA_integer_)
  cmod <- lm(y ~ x, cmod_data)
  nc <- sum(summary(cmod)$df[1L:2L])

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_warning(propertee:::.get_b11(dmod, cluster = c("uoa1", "uoa2")),
                 "are found to have NA's")
  expect_equal(suppressWarnings(propertee:::.get_b11(dmod, cluster = c("uoa1", "uoa2"),
                                         type = "HC0", cadjust = FALSE)),
               crossprod(sandwich::estfun(cmod))) # there should be no clustering

  # check case where one clustering column doesn't only have NA's
  cmod_data$uoa1 <- rep(seq(6, 10), each = 20)
  cmod <- lm(y ~ x, cmod_data)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))

  expect_warning(propertee:::.get_b11(dmod, cluster = c("uoa1", "uoa2")),
                 "have NA's for some but not all")
})

test_that(".get_b11 returns correct B_11 for multiple cluster columns", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("uoa1", "uoa2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(propertee:::.get_b11(m), expected)
})

test_that(".get_b11 returns correct B_11 for lm cmod (HC1)", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod)))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("uoa1", "uoa2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(propertee:::.get_b11(m), expected)
})

test_that(".get_b11 returns correct B_11 for glm object (HC0)", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::poisson())
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod),
        family = stats::poisson()))

  uoas <- factor(Reduce(function(x, y) paste(x, y, sep = "_"), simdata[, c("uoa1", "uoa2")]))
  nuoas <- length(levels(uoas))

  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), uoas, colSums))) *
      nuoas / (nuoas - 1L)
  )
  expect_equal(propertee:::.get_b11(m), expected)
})

test_that(paste(".get_b11 returns correct B_11 for experimental data that is a",
                "subset of cov model data (also tests NA's in some cluster",
                "but not all cluster columns)"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  nc <- sum(summary(cmod)$df[1L:2L])

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata, subset = simdata$uoa2 == 1)
  weighted_design <- ate(des, data = simdata[simdata$uoa2 == 1,])
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata[simdata$uoa2 == 1,],
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
  expect_equal(propertee:::.get_b11(m), expected)
})


test_that(paste(".get_b11 returns correct B_11 for experimental data that has",
                "no overlap with cov model data"), {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
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
  expect_equal(propertee:::.get_b11(m, cadjust = FALSE), expected)
})

test_that(".get_b11 returns expected B_11 when cmod fit to one cluster", {
  data(simdata)
  cmod <- lm(y ~ x, simdata, subset = uoa1 == 1)
  des <- rct_design(z ~ cluster(uoa1), simdata, subset = simdata$uoa1 %in% c(1, 5))

  msk <- simdata$uoa1 %in% c(1, 5)
  m <- as.lmitt(
    lm(y ~ assigned(), simdata[msk,],
       weights = ate(des, data = simdata[msk,]),
       offset = cov_adj(cmod, newdata = simdata[msk,]))
  )
  expect_warning(propertee:::.get_b11(m),
                 "meat matrix numerically indistinguishable")
  expect_equal(suppressWarnings(propertee:::.get_b11(m)),
               matrix(0, nrow = 2, ncol = 2))
})

test_that(".get_a21 returns correct matrix for lm cmod and lm damod", {
  data(simdata)

  new_df <- simdata
  new_df$uid <- seq_len(nrow(new_df))
  cmod <- lm(y ~ x, new_df)
  des <- rct_design(z ~ unitid(uid), data = new_df)

  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), new_df, weights = ate(des), offset = cov_adj(cmod)
  ))
  m_lmitt.form <- lmitt(y ~ 1, design = des, data = new_df,
                        weights = ate(des), offset = cov_adj(cmod))

  Qmat <- m_as.lmitt$weights * stats::model.matrix(m_as.lmitt)
  Cmat <- stats::model.matrix(cmod)
  nq <- nrow(stats::model.frame(m_as.lmitt))

  a21_as.lmitt <- propertee:::.get_a21(m_as.lmitt)
  expect_equal(dim(a21_as.lmitt), c(2, 2))
  expect_true(all.equal(a21_as.lmitt, crossprod(Qmat, Cmat) / nq,
                        check.atrributes = FALSE))

  a21_lmitt.form <- propertee:::.get_a21(m_lmitt.form)
  expect_equal(dim(a21_lmitt.form), c(2, 2))
  expect_true(all.equal(a21_lmitt.form, crossprod(Qmat, Cmat) / nq,
                        check.attributes = FALSE))
})

test_that(".get_a21 returns correct matrix for glm cmod and lm damod", {
  data(simdata)
  new_df <- simdata
  new_df$id <- seq(nrow(new_df))
  new_df$bin_y <- rbinom(nrow(new_df), 1, round(1 / (1 + exp(-new_df$x))))

  cmod <- suppressWarnings(glm(bin_y ~ x + force, data = new_df, family = stats::binomial()))
  des <- rct_design(z ~ unitid(id), new_df)
  m_as.lmitt <- as.lmitt(lm(
    bin_y ~ assigned(), new_df, weights = ate(des), offset = cov_adj(cmod, by = "id")
  ))
  m_lmitt.form <- lmitt(bin_y ~ 1, design = des, data = new_df,
                        weights = ate(des), offset = cov_adj(cmod, by = "id"))

  Qmat <- stats::model.matrix(m_as.lmitt)
  Cmat <- cmod$prior.weights * cmod$family$mu.eta(cmod$linear.predictors) * stats::model.matrix(cmod)
  nq <- nrow(stats::model.frame(m_as.lmitt))

  a21_as.lmitt <- propertee:::.get_a21(m_as.lmitt)
  expect_equal(dim(a21_as.lmitt), c(2, 3))
  expect_true(all.equal(a21_as.lmitt, crossprod(Qmat, Cmat) / nq,
                        check.atrributes = FALSE))

  a21_lmitt.form <- propertee:::.get_a21(m_lmitt.form)
  expect_equal(dim(a21_lmitt.form), c(2, 3))
  expect_true(all.equal(a21_lmitt.form, crossprod(Qmat, Cmat) / nq,
                        check.attributes = FALSE))
})

test_that(".get_tilde_a21 correct values", {
  set.seed(438)
  data(simdata)

  new_df <- simdata
  nc <- 30
  nq <- nrow(new_df)
  n <- nc + nq
  cmod_data <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_data)

  new_df$id <- seq(nq)
  des <- rct_design(z ~ unitid(id), new_df)

  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), new_df, weights = ate(des), offset = cov_adj(cmod, by = "id")
  ))
  m_lmitt.form <- lmitt(y ~ 1, design = des, data = new_df,
                        weights = ate(des), offset = cov_adj(cmod, by = "id"))

  expect_equal(.get_tilde_a21(m_as.lmitt), nq / n * .get_a21(m_as.lmitt))
  expect_equal(.get_tilde_a21(m_lmitt.form), nq / n * .get_a21(m_lmitt.form))
})

test_that(paste(".get_a21 returns correct matrix when data input for lmitt",
                "has NA values for some covariate values"), {
  data(simdata)
  simdata[simdata$uoa1 == 1, "x"] <- NA_real_
  cmod_data <- subset(simdata, uoa1 %in% c(2, 4))
  m_data <- subset(simdata, uoa1 %in% c(1, 3, 5))
  cmod <- lm(y ~ x, data = cmod_data)
  des <- rct_design(z ~ cluster(uoa1, uoa2), m_data)

  # warning procs twice
  expect_warning(expect_warning(
    m_as.lmitt <- lmitt(lm(y ~ assigned(), m_data, offset = cov_adj(cmod, design = des))),
    "covariance adjustments are NA"), "covariance adjustments are NA")
  expect_warning(
    m_lmitt.form <- lmitt(y ~ 1, design = des, data = m_data, offset = cov_adj(cmod)),
    "covariance adjustments are NA"
  )

  damod_mm <- stats::model.matrix(m_as.lmitt)
  pg <- stats::model.matrix(formula(cmod), m_data)
  nq <- nrow(damod_mm)

  a21_as.lmitt <- suppressWarnings(propertee:::.get_a21(m_as.lmitt))
  expect_true(all.equal(a21_as.lmitt, crossprod(damod_mm, pg) / nq,
                        check.attributes = FALSE))
  a21_lmitt.form <- suppressWarnings(propertee:::.get_a21(m_lmitt.form))
  expect_true(all.equal(a21_lmitt.form, crossprod(damod_mm, pg) / nq,
                        check.attributes = FALSE))
})

test_that(paste(".get_a21 returns correct matrix when data input for lmitt has
                NA values for some treatment assignments"), {
  data(simdata)
  simdata[simdata$uoa1 %in% c(2, 4), "z"] <- NA_integer_
  cmod <- lm(y ~ x, data = simdata, subset = uoa1 %in% c(2, 4))
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
  m_as.lmitt <- lmitt(lm(y ~ assigned(), simdata, offset = cov_adj(cmod, design = des)))
  m_lmitt.form <- lmitt(y ~ 1, simdata, design = des, offset = cov_adj(cmod))

  damod_mm <- stats::model.matrix(m_as.lmitt)
  pg <- stats::model.matrix(formula(cmod), simdata)[!is.na(simdata$z),]
  nq <- nrow(damod_mm)

  a21_as.lmitt <- suppressWarnings(propertee:::.get_a21(m_as.lmitt))
  expect_true(all.equal(a21_as.lmitt, crossprod(damod_mm, pg) / nq,
                        check.attributes = FALSE))
  a21_lmitt.form <- suppressWarnings(propertee:::.get_a21(m_lmitt.form))
  expect_true(all.equal(a21_lmitt.form, crossprod(damod_mm, pg) / nq,
                        check.attributes = FALSE))
})

test_that(".get_a21 returns only full rank columns for less than full rank model", {
  data(simdata)
  copy_simdata <- simdata
  copy_simdata$o_fac <- as.factor(copy_simdata$o)
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), copy_simdata)

  ### lmitt.formula
  damod <- lmitt(y ~ o_fac, data = copy_simdata, design = des, offset = cov_adj(cmod))
  expect_equal(dim(a21 <- .get_a21(damod)),
               c(damod$rank, damod$model$`(offset)`@fitted_covariance_model$rank))
  keep_ix <- damod$qr$pivot[1L:damod$rank]
  expect_equal(rownames(a21), colnames(model.matrix(damod))[keep_ix])
})

test_that("propertee:::.vcov_CR0 returns px2 matrix", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), simdata, weights = ate(des), offset = cov_adj(cmod)))

  vmat <- propertee:::.vcov_CR0(m)
  expect_equal(dim(vmat), c(2, 2))
})

test_that(".vcov_CR0 doesn't accept `type` argument", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), simdata, weights = ate(des), offset = cov_adj(cmod)))

  expect_error(.vcov_CR0(m, type = "CR0"), "Cannot override the `type`")
})

test_that(".vcov_CR0 returns teeMod model sandwich if it has no SandwichLayer", {
  data(simdata)

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  m <- lmitt(y ~ 1, data = simdata, design = des, weights = ate(des))

  uoas <- apply(simdata[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))
  expect_true(all.equal(vcov_tee(m, type = "CR0"),
                        sandwich::sandwich(m,
                                           meat. = sandwich::meatCL,
                                           cluster = uoas),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR0 lm w/o clustering",
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
  dmod_as.lmitt <- as.lmitt(
    lm(damod_form, data = df, subset = !is.na(uid), weights = ate(des),
       offset = cov_adj(cmod))
  )
  dmod_lmitt.form <- lmitt(y ~ 1, design = des, data = df, subset = !is.na(uid),
                           weights = ate(des), offset = cov_adj(cmod))
  onemod <- lm(y ~ x1 + x2 + z, data = df, subset = !is.na(uid),
               weights = ate(des))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(dmod_as.lmitt)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$uid),])
  nc <- nrow(Xstar)
  nq <- nrow(Z)
  n <- nc + nq

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(dmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(dmod_as.lmitt), crossprod(Z, X * dmod_as.lmitt$weights) / nq)
  expect_equal(bread. <- sandwich::bread(dmod_as.lmitt),
               n * solve(crossprod(Z * dmod_as.lmitt$weights, Z)))

  ef_damod <- utils::getS3method("estfun", "lm")(dmod_as.lmitt)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = nc, ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat. <- crossprod(estfun(dmod_as.lmitt)) / n,
               crossprod(ef_damod - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21)) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meat(dmod_as.lmitt, adjust = FALSE))

  # with perfect group balance, a22inv %*% a21 ((Z'WZ)^(-1)Z'WX) should be 0
  # for the trt row
  expect_equal((bread. %*% a21)[2,],
               setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- propertee:::.vcov_CR0(dmod_as.lmitt, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1/n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be equal to the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(dmod_as.lmitt),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR0(dmod_lmitt.form, cluster = seq_len(n), cadjust = FALSE),
                        check.attributes = FALSE))

  # var_hat(z) should be smaller than var_hat(z) from onemod
  expect_true(all(diag(vmat) <
                  diag(vcov(onemod))[c(1, 4)]))
})

test_that(paste("HC0 .vcov_CR0 lm w/o clustering",
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
  dmod_as.lmitt <- as.lmitt(
    lm(damod_form, data = df, subset = !is.na(uid), weights = ate(des),
       offset = cov_adj(cmod))
  )
  dmod_lmitt.form <- lmitt(y ~ 1, design = des, data = df, subset = !is.na(uid),
                           weights = ate(des), offset = cov_adj(cmod))
  onemod <- lm(y ~ x1 + x2 + z, data = df, subset = !is.na(uid),
               weights = ate(des))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(dmod_as.lmitt)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$uid),])
  nc <- nrow(Xstar)
  nq <- nrow(Z)
  n <- nc + nq

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(dmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(dmod_as.lmitt), crossprod(Z, X * dmod_as.lmitt$weights) / nq)
  expect_equal(bread. <- sandwich::bread(dmod_as.lmitt),
               n * solve(crossprod(Z * dmod_as.lmitt$weights, Z)))

  ef_damod <- utils::getS3method("estfun", "lm")(dmod_as.lmitt)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = nc, ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat. <- crossprod(estfun(dmod_as.lmitt)) / n,
               crossprod(ef_damod - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21)) / n)


  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meat(dmod_as.lmitt, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- propertee:::.vcov_CR0(dmod_as.lmitt, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1/n) * bread. %*% meat. %*% bread.,
                        check.attributes = FALSE))

  # vmat should be equal the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(dmod_as.lmitt),
                       check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR0(dmod_lmitt.form, cluster = seq_len(n), cadjust = FALSE),
                        check.attributes = FALSE))

  # var_hat(z) should be smaller than var_hat(z) from onemod
  expect_true(all(diag(vmat) <
                    diag(vcov(onemod))[c(1, 4)]))
})

test_that(paste("HC0 .vcov_CR0 lm w/ clustering",
                "balanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  set.seed(50)
  nclusts_C <- nclusts_Q <- 4
  MI <- 16

  # trt variable
  z <- c(rep(rep(c(0, 1), each = MI), nclusts_C / 2),
         rep(rep(c(0, 1), each = MI), nclusts_Q / 2))

  # p and mu for each matched pair's categorical and continuous cov
  px1 <- round(stats::runif((nclusts_C + nclusts_Q) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((nclusts_C) / 2, 55, 90))
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
  error_sd <- round(stats::runif(nclusts_C + nclusts_Q, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(nclusts_C + nclusts_Q) - 1)) {
    msk <- (1 + i * MI):((i + 1) * MI)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = MI) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps
  df$cid <- c(rep(NA_integer_, nclusts_C * MI), rep(seq_len(nclusts_Q), each = MI))

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()

  cmod <- lm(cmod_form, df, subset = is.na(cid))
  des <- rct_design(z ~ cluster(cid), df, subset = !is.na(df$cid))
  dmod_as.lmitt <- as.lmitt(
    lm(damod_form, data = df, subset = !is.na(cid), weights = ate(des),
       offset = cov_adj(cmod))
  )
  dmod_lmitt.form <- lmitt(y ~ 1, design = des, data = df, subset = !is.na(cid),
                           weights = ate(des), offset = cov_adj(cmod))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(dmod_as.lmitt)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$cid),])
  nq <- nrow(Z)
  nc <- nrow(Xstar)
  n <- nc + nq

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- propertee:::.get_a11_inverse(dmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- propertee:::.get_a21(dmod_as.lmitt), crossprod(Z, X * dmod_as.lmitt$weights) / nq)
  expect_equal(bread. <- sandwich::bread(dmod_as.lmitt),
               n * solve(crossprod(Z * dmod_as.lmitt$weights, Z)))

  ids <- c(df$cid[!is.na(df$cid)],
           paste0(length(df$cid[!is.na(df$cid)]) + seq_len(length(df$cid[is.na(df$cid)])),
                  "*"))
  ef_damod <- utils::getS3method("estfun", "lm")(dmod_as.lmitt)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = nc, ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat. <- crossprod(Reduce(rbind, by(estfun(dmod_as.lmitt), ids, colSums))) / n,
               crossprod(Reduce(
                 rbind,
                 by(ef_damod - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21), ids, colSums))) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meatCL(dmod_as.lmitt, cluster = ids, cadjust = FALSE))

  # with perfect group balance, a22inv %*% a21 should have 0's in the trt row
  expect_equal((bread. %*% a21)[2,], setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- propertee:::.vcov_CR0(dmod_as.lmitt, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% bread.,
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR0(dmod_lmitt.form, cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be equal the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(dmod_as.lmitt,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR0 lm w/ clustering",
                "imbalanced grps",
                "no cmod/damod data overlap", sep = ", "), {
  set.seed(50)
  nclusts_C <- nclusts_Q <- 4
  MI <- 16

  # trt variable
  z <- c(rep(rep(c(0, 1), each = MI), nclusts_C / 2),
         rep(rep(c(0, 1), each = MI), nclusts_Q / 2))

  # p and mu for each matched pair's categorical and continuous cov
  px1 <- round(stats::runif((nclusts_C + nclusts_Q) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((nclusts_C + nclusts_Q) / 2, 55, 90))
  x1 <- c(mapply(function(x) c(rbinom(MI, 1, x), rbinom(MI, 1, x)), px1))
  x2 <- c(Reduce(cbind, Map(function(x) stats::rgamma(MI * 2, x), mux2)))

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2)

  # generate clustered errors
  error_sd <- round(stats::runif(nclusts_C + nclusts_Q, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(nclusts_C + nclusts_Q) - 1)) {
    msk <- (1 + i * MI):((i + 1) * MI)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = MI) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps
  df$cid <- c(rep(NA_integer_, nclusts_C * MI), rep(seq_len(nclusts_Q), each = MI))

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()

  cmod <- lm(cmod_form, df, subset = is.na(cid))
  des <- rct_design(z ~ cluster(cid), df, subset = !is.na(df$cid))
  dmod_as.lmitt <- as.lmitt(
    lm(damod_form, data = df, subset = !is.na(cid), weights = ate(des),
       offset = cov_adj(cmod))
  )
  dmod_lmitt.form <- lmitt(y ~ 1, data = df, design = des, subset = !is.na(cid),
                           weights = ate(des), offset = cov_adj(cmod))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(dmod_as.lmitt)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$cid),])
  nc <- nrow(Xstar)
  nq <- nrow(Z)
  n <- nc + nq


  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(dmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(dmod_as.lmitt), crossprod(Z, X * dmod_as.lmitt$weights) / nq)
  expect_equal(bread. <- sandwich::bread(dmod_as.lmitt),
               n * solve(crossprod(Z * dmod_as.lmitt$weights, Z)))

  ids <- c(df$cid[!is.na(df$cid)],
           paste0(length(df$cid[!is.na(df$cid)]) + seq_len(length(df$cid[is.na(df$cid)])),
                  "*"))
  ef_damod <- utils::getS3method("estfun", "lm")(dmod_as.lmitt)
  ef_damod <- rbind(ef_damod, matrix(0, nrow = nc, ncol = ncol(ef_damod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nq, ncol = ncol(ef_cmod)), ef_cmod)
  expect_equal(meat. <- crossprod(Reduce(rbind, by(estfun(dmod_as.lmitt), ids, colSums))) / n,
               crossprod(Reduce(
                 rbind,
                 by(ef_damod -  nq / nc * ef_cmod %*% t(a11inv) %*% t(a21), ids, colSums))) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meatCL(dmod_as.lmitt, cluster = ids, cadjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR0(dmod_as.lmitt, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% bread.,
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR0(dmod_lmitt.form, cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(dmod_as.lmitt,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR0 lm w/o clustering",
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
  dmod_as.lmitt <- as.lmitt(
    lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
  )
  dmod_lmitt.form <- lmitt(y ~ 1, data = df, design = des, weights = ate(des),
                           offset = cov_adj(cmod))

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(dmod_as.lmitt)
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df)
  nc <- nrow(Xstar)
  nq <- n <- nrow(Z)

  expect_equal(a11inv <- .get_a11_inverse(dmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(dmod_as.lmitt), crossprod(Z, X * dmod_as.lmitt$weights) / nq)
  expect_equal(bread. <- sandwich::bread(dmod_as.lmitt),
               n * solve(crossprod(Z * dmod_as.lmitt$weights, Z)))

  ef_damod <- utils::getS3method("estfun", "lm")(dmod_as.lmitt)
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_damod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
  expect_equal(meat. <- crossprod(estfun(dmod_as.lmitt)) / n,
               crossprod(ef_damod - nq / nc * ef_cmod %*% a11inv %*% t(a21)) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meat(dmod_as.lmitt, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR0(dmod_as.lmitt, cluster = df$uid, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR0(dmod_lmitt.form, cluster = df$uid, cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(dmod_as.lmitt, adjust = FALSE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR0 lm w/ clustering",
                "imbalanced grps",
                "cmod is a strict subset of damod data", sep = ", "), {
  set.seed(50)
  nclusts_C <- nclusts_Q <- 4
  MI <- 16

  # trt variable
  z <- c(rep(rep(c(0, 1), each = MI), nclusts_C / 2),
         rep(rep(c(0, 1), each = MI), nclusts_Q / 2))

  # p and mu for each matched pair's categorical and continuous cov
  px1 <- round(stats::runif((nclusts_C + nclusts_Q) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((nclusts_C) / 2, 55, 90))
  x1 <- c(mapply(function(x) c(rbinom(MI, 1, x), rbinom(MI, 1, x)), px1))
  x2 <- c(Reduce(cbind, Map(function(x) stats::rgamma(MI * 2, x), mux2)))

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2)

  # generate clustered errors
  error_sd <- round(stats::runif(nclusts_C + nclusts_Q, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(nclusts_C + nclusts_Q) - 1)) {
    msk <- (1 + i * MI):((i + 1) * MI)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = MI) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  df$y <- stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps
  df$cid <- rep(seq_len(nclusts_C + nclusts_Q), each = MI)

  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()
  cmod_idx <- df$z == 0

  cmod <- lm(cmod_form, df[cmod_idx,])
  des <- rct_design(z ~ cluster(cid), df)
  dmod_as.lmitt <- as.lmitt(
    lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
  )
  dmod_lmitt.form <- lmitt(y ~ 1, data = df, design = des, weights = ate(des),
                           offset = cov_adj(cmod))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(dmod_as.lmitt)
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df)
  nc <- nrow(Xstar)
  nq <- n <- nrow(Z)

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(dmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(dmod_as.lmitt), crossprod(Z, X * dmod_as.lmitt$weights) / nq)
  expect_equal(bread. <- sandwich::bread(dmod_as.lmitt),
               n * solve(crossprod(Z * dmod_as.lmitt$weights, Z)))

  ids <- df$cid
  ef_damod <- utils::getS3method("estfun", "lm")(dmod_as.lmitt)
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_damod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
  expect_equal(meat. <- crossprod(Reduce(rbind, by(estfun(dmod_as.lmitt), ids, colSums))) / n,
               crossprod(Reduce(
                 rbind,
                 by(ef_damod - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21), ids, colSums))) / n)

  # meat should be the same as the output of sandwich::meatCL
  expect_equal(meat., sandwich::meatCL(dmod_as.lmitt, cluster = ids, cadjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR0(dmod_as.lmitt, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR0(dmod_lmitt.form, cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(dmod_as.lmitt,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR0 binomial glm cmod",
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
  dmod_as.lmitt <- lmitt(
    lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
  )
  dmod_lmitt.form <- lmitt(y ~ 1, data = df, design = des, weights = ate(des),
                           offset = cov_adj(cmod))

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(cmod)
  fit_wstar <- cmod$weights
  rstar <- cmod$residuals
  mu_etastar <- cmod$family$mu.eta(cmod$linear.predictors)
  Z <- stats::model.matrix(dmod_as.lmitt)
  X <- stats::model.matrix(formula(stats::delete.response(terms(cmod))), df)
  w <- dmod_as.lmitt$weights
  r <- dmod_as.lmitt$residuals
  mu_eta <- cmod$family$mu.eta(drop(X %*% cmod$coefficients))
  nc <- nrow(Xstar)
  nq <- n <- nrow(Z)

  expect_equal(a11inv <- .get_a11_inverse(dmod_as.lmitt), nc * solve(crossprod(Xstar * sqrt(fit_wstar))))
  expect_true(all.equal(a21 <- .get_a21(dmod_lmitt.form),
                        crossprod(Z, X * w * mu_eta) / nq,
                        check.attributes = FALSE))
  expect_equal(bread. <- sandwich::bread(dmod_as.lmitt), n * solve(crossprod(Z * w, Z)))

  ef_damod <- utils::getS3method("estfun", "lm")(dmod_as.lmitt)
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_damod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
  expect_equal(meat. <- crossprod(estfun(dmod_as.lmitt)) / n,
               crossprod(ef_damod - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21)) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meat(dmod_as.lmitt, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR0(dmod_as.lmitt, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR0(dmod_lmitt.form, cluster = seq(n), cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(dmod_as.lmitt, adjust = FALSE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR0 binomial glm cmod",
                "w/ clustering",
                "imbalanced grps",
                "damod is a strict subset of cmod data", sep = ", "), {
  set.seed(50)
  nclusts_C <- nclusts_Q <- 4
  MI <- 16

  # trt variable
  z <- c(rep(rep(0, MI * 2), nclusts_C / 2),
         rep(rep(c(0, 1), each = MI), nclusts_Q / 2))

  # p and mu for each matched pair's categorical and continuous cov
  px1 <- round(stats::runif((nclusts_C + nclusts_Q) / 2, 0.05, 0.94), 1)
  mux2 <- round(stats::runif((nclusts_C + nclusts_Q) / 2, 55, 90))
  x1 <- c(mapply(function(x) c(rbinom(MI, 1, x), rbinom(MI, 1, x)), px1))
  x2 <- c(Reduce(cbind, Map(function(x) stats::rgamma(MI * 2, x), mux2)))

  df <- data.frame("z" = z, "x1" = x1, "x2" = x2,
                   "uid" = paste0("u", seq_len(nclusts_C * MI + nclusts_Q * MI)))

  # generate clustered errors
  error_sd <- round(stats::runif(nclusts_C + nclusts_Q, 1, 3), 1)
  icc <- 0.2
  eps <- stats::rnorm(nrow(df))
  Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for (i in (seq_len(nclusts_C + nclusts_Q) - 1)) {
    msk <- (1 + i * MI):((i + 1) * MI)
    Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = MI) + icc
  }
  A <- chol(Sigma)
  eps <- t(A) %*% eps

  # generate y
  theta <- c(65, 1.5, -0.01, -2) # intercept, x1, x2, and treatment coeffs
  py <- 1 / (1 + exp(-scale(stats::model.matrix(~ x1 + x2 + z, df),
                            scale = FALSE) %*% theta))
  df$y <- rbinom((nclusts_C + nclusts_Q) * MI, 1, py)
  df$cid <- rep(seq_len(nclusts_C + nclusts_Q), each = MI)

  # set trt to NA for cmod rows
  df[1:(nclusts_C * MI), "z"] <- NA_integer_

  # model
  cmod_form <- y ~ x1 + x2
  damod_form <- y ~ assigned()
  cmod <- glm(cmod_form, data = df, family = stats::binomial())
  des <- rct_design(z ~ uoa(cid), df[!is.na(df$z),])
  dmod_as.lmitt <- lmitt(
    lm(damod_form, data = df[!is.na(df$z),], weights = ate(des), offset = cov_adj(cmod, by = "uid"))
  )
  dmod_lmitt.form <- lmitt(y ~ 1, data = df[!is.na(df$z),], design = des,
                           weights = ate(des), offset = cov_adj(cmod, by = "uid"))

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(cmod)
  fit_wstar <- cmod$weights
  rstar <- cmod$residuals
  mu_etastar <- cmod$family$mu.eta(cmod$linear.predictors)
  Z <- stats::model.matrix(dmod_as.lmitt)
  X <- stats::model.matrix(formula(stats::delete.response(terms(cmod))), df[!is.na(df$z),])
  w <- dmod_as.lmitt$weights
  r <- dmod_as.lmitt$residuals
  mu_eta <- cmod$family$mu.eta(drop(X %*% cmod$coefficients))
  nq <- nrow(Z)
  nc <- nrow(Xstar)
  n <- nrow(Xstar)

  expect_equal(a11inv <- .get_a11_inverse(dmod_as.lmitt), nc * solve(crossprod(Xstar * sqrt(fit_wstar))))
  expect_equal(a21 <- .get_a21(dmod_as.lmitt), crossprod(Z, X * w * mu_eta) / nq)
  expect_equal(bread. <- sandwich::bread(dmod_as.lmitt), n * solve(crossprod(Z * w, Z)))

  ef_order <- c(gsub("u", "", sort(df$uid[!is.na(df$z)])), rownames(df[is.na(df$z),]))
  ids <- factor(df$cid, levels = unique(df$cid))[as.numeric(ef_order)]
  ef_cmod <- estfun(cmod)[ef_order,]

  unextended_ef_damod <- estfun(as(dmod_as.lmitt, "lm"))
  ef_damod <- matrix(0, nrow = nrow(ef_cmod), ncol = ncol(unextended_ef_damod),
                     dimnames = list(rownames = seq_len(nrow(ef_cmod)),
                                     colnames = colnames(unextended_ef_damod)))
  Q_order <- match(sort(df$uid[!is.na(df$z)]), df$uid[!is.na(df$z)])
  ef_damod[1L:sum(!is.na(df$z)),] <- estfun(as(dmod_as.lmitt, "lm"))[Q_order,]
  expect_equal(meat. <- crossprod(Reduce(rbind, by(estfun(dmod_as.lmitt), ids, colSums))) / n,
               (crossprod(Reduce(rbind, by(ef_damod, ids, colSums))) -
                  crossprod(Reduce(rbind, by(ef_damod, ids, colSums)),
                            Reduce(rbind, by(nq / sqrt(n * nc) * ef_cmod, ids, colSums))) %*% a11inv %*% t(a21) -
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(nq / sqrt(n * nc) * ef_cmod, ids, colSums)),
                                               Reduce(rbind, by(ef_damod, ids, colSums))) +
                  a21 %*% a11inv %*% crossprod(Reduce(rbind, by(nq / sqrt(n * nc) * ef_cmod, ids, colSums))) %*%
                  a11inv %*% t(a21)
               ) / n)

  # meat should be the same as the output of sandwich::meatCL
  expect_equal(meat., sandwich::meatCL(dmod_as.lmitt, cluster = ids, cadjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR0(dmod_as.lmitt, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR0(dmod_lmitt.form, cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(dmod_as.lmitt,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
})

test_that("HC1/CR1/MB1", {
  set.seed(8431432)
  # no clustering
  n <- 50
  Sigma <- diag(runif(n, 0.01, 2))
  x1 <- rnorm(n)
  z <- rep(c(0, 1), each = n / 2)
  y <- -0.5 * x1 -0.1 * z + chol(Sigma) %*% rnorm(n)
  dat <- data.frame(x1=x1, z=z, y=y, c=sample(seq_len(n)))

  des <- rct_design(z ~ unitid(c), dat)
  cmod <- lm(y ~ x1, dat)
  dmod <- lmitt(y ~ 1, data = dat, design = des,
                weights = ate(des), offset = cov_adj(cmod))

  # CR1
  expect_true(
    all.equal(cr1_vmat <- vcov_tee(dmod, type = "CR1"),
              50 / 48 * .vcov_CR0(dmod, cluster = .make_uoa_ids(dmod, "CR"), cadjust=FALSE),
              check.attributes = FALSE)
  )
  expect_equal(attr(cr1_vmat, "type"), "CR1")
  expect_true(any(grepl("CR1", capture.output(summary(dmod, vcov.type = "CR1")))))

  # HC1
  expect_true(all.equal(hc1_vmat <- vcov_tee(dmod, type = "HC1"), cr1_vmat, check.attributes = FALSE))
  expect_equal(attr(hc1_vmat, "type"), "HC1")
  expect_true(any(grepl("HC1", capture.output(summary(dmod, vcov.type = "HC1")))))

  # MB1
  expect_true(all.equal(mb1_vmat <- vcov_tee(dmod, type = "MB1"), cr1_vmat, check.attributes = FALSE))
  expect_equal(attr(mb1_vmat, "type"), "MB1")
  expect_true(any(grepl("MB1", capture.output(summary(dmod, vcov.type = "MB1")))))

  # clustering
  n <- 50
  g <- 5
  icc <- 0.1
  Sigma <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(g)) {
    Sigma[((i - 1) * n / g + 1):(i * (n / g)), ((i - 1) * n / g + 1):(i * (n / g))] <- (
      (1-icc) * diag(1, nrow=n/g, ncol=n/g) + icc
    )
  }
  x1 <- rnorm(n)
  z <- rep(c(0, 1), c((g-2) * n/g, (g-3) * n/g))
  y <- -0.5 * x1 -0.1 * z + chol(Sigma) %*% rnorm(n)
  dat <- data.frame(x1=x1, z=z, y=y, c=rep(seq_len(g), each = n/g))

  des <- rct_design(z ~ unitid(c), dat)
  cmod <- lm(y ~ x1, dat)
  dmod <- lmitt(y ~ 1, data = dat, design = des,
                weights = ate(des), offset = cov_adj(cmod))

  # CR1
  expect_true(
    all.equal(cr1_vmat <- vcov_tee(dmod, type = "CR1"),
              g / (g-1) * (n-1) / (n-2) * .vcov_CR0(dmod, cluster = .make_uoa_ids(dmod, "CR"),
                                                      cadjust=FALSE),
              check.attributes = FALSE)
  )
  expect_equal(attr(cr1_vmat, "type"), "CR1")
  expect_true(any(grepl("CR1", capture.output(summary(dmod, vcov.type = "CR1")))))
})

test_that("type attribute", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
  damod <- lmitt(y ~ 1, data = simdata, weights = "ate", design = des)
  expect_identical(attr(vcov(damod), "type"), "CR0")
  expect_identical(attr(vcov(damod, type = "CR0"), "type"), "CR0")
  expect_identical(attr(vcov(damod, type = "MB0"), "type"), "MB0")
  expect_identical(attr(vcov(damod, type = "HC0"), "type"), "HC0")
})

test_that("#119 flagging vcov_tee entries as NA", {
  ### factor moderator variable
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
  mod <- lmitt(y ~ as.factor(o), data = simdata, design = des)
  expect_warning(vc <- vcov_tee(mod),
                 "will be returned as NA: as.factor(o)1, as.factor(o)3",
                 fixed = TRUE)

  # Issue is in subgroups w/ moderator=1:z=0, moderator=1:z=1, and
  # moderator=3, so .check_df_moderator_estimates should NA those vcov entries
  na_dim <- c(1, 3, 5)
  expect_true(all(
    abs(diag(sandwich::sandwich(mod, meat. = sandwich::meatCL,
                                cluster = .make_uoa_ids(mod, "CR")))[na_dim])
     < .Machine$double.eps)
  )
  expect_true(all(is.na(vc[na_dim, ])))
  expect_true(all(is.na(vc[, na_dim])))
  expect_true(all(!is.na(vc[-na_dim, -na_dim])))

  ### different factor moderator variable (may not produce negative diagonals
  ### but subgroups with moderator=1:z=0 and moderator=1:z=1 should be
  ### numerically 0)
  set.seed(37)
  ddata <- data.frame(modr = factor(c(rep(1, 8), rep(2, 12))),
                      y = rnorm(20),
                      z = c(rep(rep(c(0, 1), each = 4), 2), rep(1, 4)),
                      ass_id = rep(seq(5), each = 4))
  ddes <- rct_design(z ~ unitid(ass_id), ddata)
  dmod <- lmitt(y ~ modr, design = ddes, data = ddata)
  expect_warning(vc <- vcov_tee(dmod),
                 "will be returned as NA: modr1",
                 fixed = TRUE)

  na_dim <- c(1, 3)
  expect_true(all(
    abs(
      diag(sandwich::sandwich(dmod, meat. = sandwich::meatCL,
                              cluster = .make_uoa_ids(dmod, "CR")))[na_dim]
    ) < .Machine$double.eps)
  )
  expect_true(all(is.na(vc[na_dim, ])))
  expect_true(all(is.na(vc[, na_dim])))
  expect_true(all(!is.na(vc[-na_dim, -na_dim])))

  ### valid factor moderator variable; vcov diagonals shouldn't be numerically 0
  ### when using the sandwich package and shouldn't have NA's when using vcov_tee
  ddata <- data.frame(modr = factor(c(rep(1, 8), rep(2, 12))),
                      y = rnorm(20),
                      z = c(rep(rep(c(0, 1), each = 4), 2), rep(1, 4)),
                      ass_id = rep(seq(10), each = 2))
  ddes <- rct_design(z ~ unitid(ass_id), ddata)
  dmod <- lmitt(y ~ modr, design = ddes, data = ddata)
  vc <- vcov_tee(dmod)
  expect_true(all(
    abs(diag(sandwich::sandwich(dmod, meat. = sandwich::meatCL,
                                cluster = .make_uoa_ids(dmod, "CR")))[na_dim])
    > .Machine$double.eps)
  )
  expect_true(all(!is.na(vc)))

  #### lmitt.lm
  ## we chose to implement special logic for subgroup-level d.f. only
  ## for lmitt objects created with lmitt.formula; the commented-out
  ## tests that follow would have tested similar logic for lmitt.lm
  ##objects. (this distinction is reflected in `lmitt.lm()`'s leaving
  ## the "moderator" slot empty for teeMod objects,
  ## so .check_df_moderator_estimates will return
  ## the initial vcov matrix, and we need not check it here.)
  # mod <- lmitt(lm(y ~ as.factor(o) + as.factor(o):assigned(des), data = simdata), design = des)
  # vc <- vcov_tee(mod)[5:7, 5:7]
  #
  # #****************************************
  # ### Setting these to NA manually only for testing purposes!
  # vc[1, ] <- NA
  # vc[, 1] <- NA
  # ### Remove these once #119 is addressed!!!!!
  # #****************************************
  #
  # # Issue is in subgroup o_fac=1, so the first entry in the vcov matrix
  # expect_true(all(is.na(vc[1, ])))
  # expect_true(all(is.na(vc[, 1])))
  # expect_true(all(!is.na(vc[-1, -1])))
  #
  # ### valid continuous moderator variable
  damod <- lmitt(y ~ o, data = simdata, design = des)
  vc <- vcov(damod)
  expect_true(all(!is.na(vc)))
  #
  # ### invalid continuous moderator variable
  simdata <- simdata[(simdata$uoa1 == 2 & simdata$uoa2 == 2) |
                      (simdata$uoa1 == 2 & simdata$uoa2 == 1),]
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
  damod <- lmitt(y ~ o, data = simdata, design = des)
  expect_warning(vc <- vcov(damod), "will be returned as NA: o")
  na_dim <- which(grepl("(Intercept)", "z._o", colnames(vc)))
  expect_true(all(is.na(vc[na_dim, ])))
  expect_true(all(is.na(vc[, na_dim])))
  expect_true(all(!is.na(vc[-na_dim, -na_dim])))
})

test_that(".check_df_moderator_estimates other warnings", {
  data(simdata)

  # fail without a teeMod model
  nodamod <- lm(y ~ x, simdata)
  vmat <- vcov(nodamod)
  cluster_ids <- apply(simdata[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))
  expect_error(.check_df_moderator_estimates(vmat, nodamod, cluster_ids),
               "must be a teeMod")

  # invalid `data` arg
  des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
  damod <- lmitt(y ~ factor(o), design = des, data = simdata)
  expect_error(.check_df_moderator_estimates(vmat, damod, cluster_ids,
                                             cbind(y = simdata$y, z = simdata$z),
                                             envir = parent.frame()),
               "`data` must be a dataframe")
})

test_that("#123 ensure PreSandwich are converted to Sandwich", {
  data(simdata)
  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
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
  v1 <- vcov_tee(dmod1)

  wts2 <- ett(des, data = copy_simdata)
  offst2 <- cov_adj(cmod, copy_simdata, by = "schoolid")
  lm2 <- lm(y ~ assigned(), copy_simdata, weights = wts2, offset = offst2)
  dmod2 <- lmitt(lm2, design = des)
  expect_true(is(damod$model$`(offset)`, "SandwichLayer"))
  v2 <- vcov_tee(dmod2)

  expect_true(all.equal(v1, v2))

})

test_that("model-based SE's cluster units of assignment in small blocks at block level", {
  desdata <- data.frame(bid = rep(c(1, 2), each = 20),
                        uoa_id = c(rep(c(1, 2), each = 10), rep(seq(3, 6), each = 5)),
                        a = c(rep(c(0, 1), each = 10), rep(rep(c(0, 1), each = 5), 2)))
  desdata$y <- rnorm(40)
  des <- rct_design(a ~ unitid(uoa_id) + block(bid), desdata)
  suppressMessages(mod <- lmitt(y ~ 1, design = des, data = desdata))
  vc_w_small_block_clusters <- vcov_tee(mod)
  vc_w_no_small_block_clusters <- .vcov_CR0(mod,
                                            cluster = .make_uoa_ids(mod, "DB"),
                                            by = "uoa_id")
  expect_true(vc_w_small_block_clusters[2, 2] != vc_w_no_small_block_clusters[2, 2])
})

test_that("#177 vcov with by argument", {
  set.seed(23)
  cmod_data <- data.frame(yr = rep(c("00", "01", "02"), 5),
                          id = rep(letters[1:5], each = 3),
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_data)
  desdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newdes <- rct_design(a ~ unitid(id), desdat)
  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, design = newdes, data = analysis_dat, offset = cov_adj(cmod))
  expect_error(vcov(mod), "not uniquely specified. Provide a `by` argument")
  expect_silent(vcov_tee(mod, by = "by_col"))
})

test_that("#177 vcov with by", {
  set.seed(23)
  cmod_data <- data.frame(yr = rep(c("00", "01", "02"), 5),
                          id = rep(letters[1:5], each = 3),
                          x = rnorm(5 * 3),
                          y = rnorm(5 * 3),
                          by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_data)
  desdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newdes <- rct_design(a ~ unitid(id), desdat)
  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, design = newdes, data = analysis_dat, offset = cov_adj(cmod, by = "by_col"))
  expect_equal(vcov(mod), vcov(mod, by = "by_col"))
})

test_that("vcov_tee does not error when asking for design-based SE for model
          with covariance adjustment but without absorbed block effects",{
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  damod_off <- lmitt(y ~ 1, design = des, data = simdata, weights = ate(des),
                     offset = cov_adj(cmod))
  expect_true(!is.na(vcov_tee(damod_off, type = "DB0")))
})

test_that("vcov_tee errors when asking for design-based SE for a model
          with moderators",{
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  damod_sub <- lmitt(y ~ dose, design = des, data = simdata, weights = ate(des))
  expect_error(
    vcov_tee(damod_sub, type = "DB0"), 
    "cannot be computed for teeMod models with moderators"
  )
})

test_that(".merge_block_id_cols correctly combines multiple block columns", {
  df <- data.frame(
    block1 = c("A", "A", "B", "B", "B", "B"),
    block2 = c(1, 1, 1, 1, 2, 2)
  )
  df_comb <- data.frame(
    block1 = c(1, 1, 2, 2, 3, 3)
  )
  expect_equal(
    propertee:::.merge_block_id_cols(df=df, ids=c("block1", "block2"))$block1,
    df_comb$block1
  )
})

test_that(".get_DB_wo_covadj_se returns correct value for designs with small blocks",{
  # generate data
  nbs <- rep(2, 10)
  n <- sum(nbs) # sample size
  B <- length(nbs) # number of blocks
  ws <- round(rnorm(n=n, mean=50, sd=10))
  yobs <- rnorm(n=n) # observed y's
  zobs <- c() # treatment assignment, 1 or 2
  for (b in 1:B){
    zobs <- c(zobs, sample(c(1,2)))
  }
  
  # calculate variance
  pi_all <- matrix(rep(1/2,2*n), nrow=n) # assignment probabilities, n by 2
  nbk_all <- matrix(rep(1,2*n), nrow=n) # nbk for all units, n by 2
  gammas <- cbind(nbk_all[,1] * ws, nbk_all[,2] * ws) # gamma, n by K
  nu <- c() # nu_b
  
  for (k in 1:2){
    indk <- (zobs == k)
    thetak <- sum(ws[indk] * yobs[indk] / pi_all[indk, k]) /
      sum(ws[indk] / pi_all[indk, k]) # Hajek estimators
    gammas[indk,k] <- gammas[indk,k] / pi_all[indk, k] * (yobs[indk] - thetak)
    for (b in 1:B){
      in_b <- (sum(nbs[1:b])-nbs[b]+1):(sum(nbs[1:b])) # indices of block b
      indbk <- (zobs[in_b] == k)
      if (k > 1){
        indb0 <- (zobs[in_b] == 1)
        nu <- c(nu, (gammas[in_b,k][indbk] - gammas[in_b,1][indb0])^2)
      }
    }
  }
  data <- data.frame(cid = 1:n, bid = rep(1:B, each=2), y = yobs, z = zobs-1, w = ws)
  des <- rct_design(z ~ cluster(cid) + block(bid), data)
  damod <- lmitt(y ~ 1, design = des, data = data, weights = ate(des) * data$w)
  expect_equal(vcov_tee(damod, type="DB0")[1,1], sum(nu) / sum(ws)^2) 
})

test_that(".get_DB_wo_covadj_se returns correct value for designs with a few large blocks",{
  # generate data
  nbs <- rep(10, 2)
  n <- sum(nbs) # sample size
  B <- length(nbs) # number of blocks
  ws <- round(rnorm(n=n, mean=50, sd=10))
  yobs <- rnorm(n=n) # observed y's
  zobs <- rep(1, n) # treatment assignment, 1 or 2
  zobs[c(sample(1:10, 5), sample(11:20, 4))] <- 2
  
  # calculate variance
  pi_all <- matrix(c(rep(1/2,n), rep(c(0.6,0.4),10)), nrow=n, byrow = TRUE)  
  # assignment probabilities, n by 2
  nbk <- matrix(c(5,6,5,4), nrow=B)  # nbk, B by 2
  gammas <- cbind(c(rep(5,10),rep(6,10)) * ws, 
                  c(rep(5,10),rep(4,10)) * ws) # gamma, n by K
  gamsbk <- matrix(nrow=B, ncol=2)  # s^2 b,j
  nu <- c()
  
  for (k in 1:2){
    indk <- (zobs == k)
    thetak <- sum(ws[indk] * yobs[indk] / pi_all[indk, k]) /
      sum(ws[indk] / pi_all[indk, k])  # Hajek estimator
    gammas[indk,k] <- gammas[indk,k] / pi_all[indk, k] * (yobs[indk] - thetak)
    for (b in 1:B){
      in_b <- (sum(nbs[1:b])-nbs[b]+1):(sum(nbs[1:b])) # indices of block b
      gamsbk[b,k] <- var(gammas[in_b,k][zobs[in_b] == k])
      if (k > 1){
        nu <- c(nu, gamsbk[b,1] / nbk[b,1] + gamsbk[b,k] / nbk[b,k])
      }
    }
  }
  data <- data.frame(cid = 1:n, bid = rep(1:B, each=10), y = yobs, z = zobs-1, w = ws)
  des <- rct_design(z ~ cluster(cid) + block(bid), data)
  damod <- lmitt(y ~ 1, design = des, data = data, weights = ate(des) * data$w)
  
  #ainv <- .get_DB_a_inverse(damod)
  #meat <- .get_DB_meat(damod)
  #vmat <- as.matrix((ainv %*% meat %*% t(ainv))[3,3])
  
  expect_equal(vcov_tee(damod, type="DB0")[1,1], sum(nu) / sum(ws)^2) 
  #expect_true(vmat[1,1] != sum(nu) / sum(ws)^2) 
})

test_that(".get_appinv_atp returns correct (A_{pp}^{-1} A_{tau p}^T)
          for tee models with absorbed intercept", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  damod_abs <- lmitt(y ~ 1, design = des, data = simdata, weights = ate(des),
                     absorb = TRUE)
  damod_abs_off <- lmitt(y ~ 1, design = des, data = simdata, weights = ate(des),
                         offset = cov_adj(cmod), absorb = TRUE)
  
  bid <- simdata$bid
  B <- cbind(as.integer(bid == 1), as.integer(bid == 2), as.integer(bid == 3))
  goal <- matrix(0, nrow = 3, ncol = 1)
  for (blk in 1:3){
    goal[blk] <- sum(damod_abs$weights * damod_abs$residuals * B[,blk]) / 
      sum(damod_abs$weights * B[,blk])
  }
  expect_true(all.equal(goal, .get_appinv_atp(damod_abs, db = TRUE)))
  
  for (blk in 1:3){
    goal[blk] <- sum(damod_abs_off$weights * damod_abs_off$residuals * B[,blk]) / 
      sum(damod_abs_off$weights * B[,blk])
  }
  expect_true(all.equal(goal, propertee:::.get_appinv_atp(damod_abs_off, db = TRUE)))
})

test_that(".get_phi_tilde returns correct grave{phi}
          for tee models with absorbed intercept", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  damod_abs <- lmitt(y ~ 1, design = des, data = simdata, weights = ate(des),
                     absorb = TRUE)
  
  bid <- simdata$bid
  B <- cbind(as.integer(bid == 1), as.integer(bid == 2), as.integer(bid == 3))
  Z <- cbind(as.integer(simdata$z == 0), as.integer(simdata$z == 1))
  
  ws <- damod_abs$weights
  p <- matrix(0, nrow = 2, ncol = 3)
  for (k in 1:2)
    for (j in 1:3){
      p[k, j] <- sum(ws * Z[,k] * B[,j]) / sum(ws * B[,j])
    }
  goal <- matrix(0, nrow = 50, ncol = 3)
  for (s in 1:3){
    goal[,s] <- ws * (Z[,2] - p[2,s]) * B[,s]
  }
  expect_true(all.equal(goal, propertee:::.get_phi_tilde(damod_abs, db = TRUE)))
})

