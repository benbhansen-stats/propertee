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
  expect_equal(vmat, suppressMessages(.vcovMB_CR0(damod)))
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

  # # test other arg types
  # expect_equal(vcovDA(dmod), vcovDA(dmod, cluster = simdata[, c("uid"), drop = FALSE]))
  # 
  # uoas <- matrix(simdata$uid)
  # colnames(uoas) <- "uid"
  # expect_equal(vcovDA(dmod), vcovDA(dmod, cluster = uoas))
  # 
  # expect_equal(vcovDA(dmod), vcovDA(dmod, cluster = list(uid = simdata$uid)))
  # 
  # expect_equal(vcovDA(dmod), vcovDA(dmod, cluster = uid))
  # 
  # uid <- simdata$uid
  # expect_equal(vcovDA(dmod), vcovDA(dmod, cluster = uid))
  # 
  # uid <- as.numeric(simdata$uid)
  # expect_equal(vcovDA(dmod), vcovDA(dmod, cluster = uid))
  
  # # can also specify different cluster level
  # expect_equal(vcovDA(dmod, cluster = c("cid1", "cid2")),
  #              vcovDA(dmod, cluster = simdata[, c("cid1", "cid2")]))
  # 
  # expect_equal(vcovDA(dmod, cluster = c("cid1", "cid2")),
  #              vcovDA(dmod, cluster = cbind(cid1 = simdata$cid1,
  #                                           cid2 = simdata$cid2)))
  # 
  # expect_equal(vcovDA(dmod, cluster = c("cid1", "cid2")),
  #              vcovDA(dmod, cluster = list(cid1 = simdata$cid1,
  #                                          cid2 = simdata$cid2)))
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
  
  expect_warning(suppressMessages(vcovDA(dmod, cluster = c("cid1", "cid2"))),
                 "these observations should be treated as IID")
  
  expected <- suppressMessages(vcovDA(dmod))
  expect_equal(
    suppressMessages(suppressWarnings(vcovDA(dmod, cluster = c("cid1", "cid2")))),
    expected
  )
})

test_that("variance helper functions fail without a DirectAdjusted model", {
  data(simdata)
  cmod <- lm(y ~ z, data = simdata)

  expect_error(.vcovMB_CR0(cmod), "must be a DirectAdjusted")
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
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = offset)
  )

  expect_error(.vcovMB_CR0(m), "must have an offset of class")
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
  msk <- !is.na(.merge_preserve_order(
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

  expect_message(b12 <- .get_b12(m), paste(nrow(simdata), "rows"))
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
  zname <- .txt_fn(m)

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
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))

  expect_message(b12 <- .get_b12(m), paste(nrow(simdata), "rows"))
  expect_equal(b12, t(cmod_eqns) %*% m_eqns)
})

test_that(".get_b12 fails with invalid custom cluster argument", {
  data(simdata)
  
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)
  
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  
  expect_error(.get_b12(m, cluster = c("cid3")),
               "cid3 are missing from the covariance adjustment model dataset")
  expect_error(.get_b12(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
})

test_that(".get_b12 produces correct estimates with valid custom cluster argument", {
  data(simdata)
  simdata[simdata$bid == 1, "z"] <- 0
  simdata[simdata$bid == 2, "z"] <- 0
  simdata[simdata$bid == 3, "z"] <- 1
  
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ block(bid) + cluster(cid1, cid2), data = simdata)

  m <- lmitt(y ~ assigned(), data = simdata, design = des, offset = cov_adj(cmod))
  
  cmod_eqns <- Reduce(
    rbind,
    by(estfun(cmod), list(simdata$cid1, simdata$cid2), colSums)
  )
  dmod_eqns <- Reduce(
    rbind,
    by(estfun(m), list(simdata$cid1, simdata$cid2), colSums)
  )
  expected <- crossprod(cmod_eqns, dmod_eqns)
  
  # default (columns specified in `cluster` argument of Design) matches expected
  expect_equal(suppressMessages(.get_b12(m)), expected)
  expect_equal(suppressMessages(.get_b12(m, cluster = c("cid1", "cid2"))), expected)
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
  
  expect_message(b12 <- .get_b12(dmod, cluster = c("cid1", "cid2")), "0 rows")
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
  
  expect_message(b12 <- .get_b12(dmod, cluster = c("cid1", "cid2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))

  # second test is where the non-NA cluster ID's overlap with design
  cmod_data$cid1 <- rep(seq(1, 5), each = 20)
  cmod <- lm(y ~ x, cmod_data)
  
  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))
  
  expect_message(b12 <- .get_b12(dmod, cluster = c("cid1", "cid2")), "0 rows")
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
    by(estfun(cmod), list(cmod_data$cid1, cmod_data$cid2), colSums)
  )
  dmod_eqns <- Reduce(
    rbind,
    by(estfun(dmod), list(simdata$cid1, simdata$cid2), colSums)
  )
  expected <- crossprod(cmod_eqns, dmod_eqns)

  expect_message(b12 <- .get_b12(dmod, cluster = c("cid1", "cid2")),
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
  zname <- .txt_fn(m)

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
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))

  expect_message(b12 <- .get_b12(m), paste(nrow(simdata[simdata$cid2 == 1,]), "rows"))
  expect_equal(b12,
               t(cmod_eqns) %*% m_eqns)
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
  zname <- .txt_fn(m)

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
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))

  expect_equal(suppressMessages(.get_b12(m)),
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

  expect_message(b12 <- .get_b12(m), "0 rows")
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
    expect_message(b12 <- .get_b12(m), paste(nrow(cmod_data), "rows")),
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
    by(sandwich::estfun(m)[msk, , drop = FALSE],
       list(cid1 = simdata$cid1[!is.na(simdata$x)][msk],
            cid2 = simdata$cid2[!is.na(simdata$x)][msk]),
       colSums))

  b12 <- suppressMessages(.get_b12(m))
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
  zname <- .txt_fn(m)

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
    sort = FALSE)[zname])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, , drop = FALSE],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))

  expect_equal(suppressMessages(.get_b12(m)),
               t(cmod_eqns) %*% m_eqns)
})

test_that(".get_a22_inverse returns correct value for lm", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(des)))

  fim <- crossprod(stats::model.matrix(m) * m$weights, stats::model.matrix(m))
  expect_equal(.get_a22_inverse(m), solve(fim))
})

test_that(".get_a22_inverse returns correct value for glm fit with Gaussian family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(glm(y ~ assigned(), data = simdata, weights = ate(des)))

  fim <- crossprod(stats::model.matrix(m) * sqrt(m$weights))
  expect_equal(.get_a22_inverse(m), solve(fim))
})

test_that(".get_a22_inverse returns correct value for glm fit with poisson family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des),
        family = stats::poisson())
  )

  fim <- crossprod(stats::model.matrix(m) * sqrt(m$weights))
  expect_equal(.get_a22_inverse(m), solve(fim))
})

test_that(".get_a22_inverse returns correct value for glm fit with quasipoisson family", {
  data(simdata)

  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(des),
        family = stats::quasipoisson())
  )

  dispersion <- sum((m$weights * m$residuals)^2) / sum(m$weights)
  fim <- crossprod(stats::model.matrix(m) * sqrt(m$weights))
  expect_equal(.get_a22_inverse(m), solve(fim))
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

  expect_equal(.get_b22(m, type = "HC0"),
               vmat * nuoas / (nuoas - 1L))
  expect_equal(.get_b22(m, type = "HC1"),
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
  
  expect_error(.get_b22(m, cluster = c("cid3")),
               "cid3 are missing from the ITT effect model dataset")
  expect_error(.get_b22(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
  
  simdata$cid3 <- NA_integer_
  des <- rct_design(z ~ cluster(cid1, cid2, cid3), data = simdata)
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  expect_error(.get_b22(m, cluster = "cid3"),
               "cannot handle NAs")
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

  expect_equal(.get_b22(m, cluster = uoanames, type = "HC0"),
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

  expect_equal(.get_b22(m, type = "HC0"),
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

  expect_equal(.get_b22(m, type = "HC0"),
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

  expect_equal(.get_b22(m, type = "HC0"),
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

  expect_equal(.get_b22(m, type = "HC0"),
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

  expect_equal(.get_b22(m, type = "HC0"),
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
  expect_equal(.get_a11_inverse(m), solve(fim), tolerance = 1e-4) # tol due to diffs from chol2inv vs. solve(crossprod)
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
  expect_equal(.get_a11_inverse(m), dispersion * solve(fim), tolerance = 1e-4) # tol due to diffs from chol2inv vs. solve(crossprod)
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

test_that(".get_b11 returns correct B_11 for one cluster column", {
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

  expect_equal(.get_b11(m), expected)
})

test_that(".get_b11 fails with invalid custom cluster argument", {
  data(simdata)
  
  cmod <- lm(y ~ x, simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  nuoas <- nrow(des@structure)
  
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )
  
  expect_error(.get_b11(m, cluster = c("cid3")),
               "cid3 are missing from the covariance adjustment model dataset")
  expect_error(.get_b11(m, cluster = c(TRUE, FALSE)),
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
  
  expect_equal(.get_b11(m, cluster = "cid1"), expected)
  
  # test different clustering level
  bids <- factor(simdata[, "bid"])
  nbids <- length(levels(bids))
  
  nc <- sum(summary(cmod)$df[1L:2L])
  expected <- (
    crossprod(Reduce(rbind, by(sandwich::estfun(cmod), simdata[, "bid"], colSums))) *
      nbids / (nbids - 1L) * (nc - 1L) / (nc - 2L)
  )
  expect_equal(.get_b11(m, cluster = "bid"), expected)
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
  
  expect_warning(.get_b11(dmod, cluster = c("cid1", "cid2")),
                 "cid1, cid2 are found to have NA's")
  expect_equal(suppressWarnings(.get_b11(dmod, cluster = c("cid1", "cid2"),
                                         type = "HC0", cadjust = FALSE)),
               crossprod(sandwich::estfun(cmod))) # there should be no clustering
  
  # check case where one clustering column doesn't only have NA's
  cmod_data$cid1 <- rep(seq(6, 10), each = 20)
  cmod <- lm(y ~ x, cmod_data)
  
  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata)
  dmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, design = des)))
  
  expect_warning(.get_b11(dmod, cluster = c("cid1", "cid2")),
                 "Only cid1 will be used to cluster")
  expect_equal(suppressWarnings(.get_b11(dmod, cluster = c("cid1", "cid2"))),
               .get_b11(dmod, cluster = c("cid1")))
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

  expect_equal(.get_b11(m), expected)
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

  expect_equal(.get_b11(m), expected)
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
    lm(y ~ assigned(), data = simdata[simdata$cid2 == 1,],
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

  expect_equal(.get_b11(m, cadjust = FALSE), expected)
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

  expect_warning(.get_b11(m),
                 "meat matrix numerically indistinguishable")
  expect_equal(suppressWarnings(.get_b11(m)),
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

  expect_message(vmat <- .vcovMB_CR0(m), paste(nrow(simdata), "rows"))
  expect_equal(dim(vmat), c(2, 2))
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

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  a11inv <- solve(crossprod(Xstar))
  expect_equal(a11inv, flexida:::.get_a11_inverse(damod))
  b22 <- crossprod(Z, diag((damod$weights * damod$residuals)^2)) %*% Z
  expect_equal(b22, flexida:::.get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(Z * damod$weights, Z))
  expect_equal(a22inv, flexida:::.get_a22_inverse(damod))
  a21 <- crossprod(Z, X * damod$weights)
  expect_equal(a21, flexida:::.get_a21(damod))
  b12 <- matrix(0, nrow = dim(Xstar)[2], ncol = dim(Z)[2])
  expect_equal(b12, suppressMessages(flexida:::.get_b12(damod)))
  b11 <- crossprod(Xstar * cmod$residuals^2, Xstar)
  expect_equal(b11, flexida:::.get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(N * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meat
  expect_equal(b22 / N, sandwich::meat(damod, adjust = FALSE))

  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  ## HEURISTIC CHECKS
  # with perfect group balance, a22inv %*% a21 should be 0 for the trt row
  expect_equal((a22inv %*% a21)[2,],
               setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_message(vmat <- .vcovMB_CR0(damod, cadjust = FALSE),
                 "0 rows")
  expect_equal(vmat,
               a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  # var_hat(z) should be equal to than that given by sandwich
  expect_equal(diag(vmat)[2],
               diag(sandwich::sandwich(damod))[2])

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

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  a11inv <- solve(crossprod(Xstar))
  expect_equal(a11inv, flexida:::.get_a11_inverse(damod))
  b22 <- crossprod(Z, diag((damod$weights * damod$residuals)^2)) %*% Z
  expect_equal(b22, flexida:::.get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(Z * damod$weights, Z))
  expect_equal(a22inv, flexida:::.get_a22_inverse(damod))
  a21 <- crossprod(Z, X * damod$weights)
  expect_equal(a21, flexida:::.get_a21(damod))
  b12 <- matrix(0, nrow = dim(Xstar)[2], ncol = dim(Z)[2])
  expect_equal(b12, suppressMessages(flexida:::.get_b12(damod)))
  b11 <- crossprod(Xstar * cmod$residuals^2, Xstar)
  expect_equal(b11, flexida:::.get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(N / 4 * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meat
  expect_equal(b22 / (N / 4), sandwich::meat(damod, adjust = FALSE))

  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_message(vmat <- .vcovMB_CR0(damod, cadjust = FALSE), "0 rows")
  expect_equal(vmat,
               a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  # var_hat(z) should be greater than that given by sandwich
  expect_true(all(diag(vmat) >
                  diag(sandwich::sandwich(damod))))

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


  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  a11inv <- solve(crossprod(Xstar))
  expect_equal(a11inv, .get_a11_inverse(damod))
  b22 <- crossprod(
    Reduce(rbind,
           by(Z * damod$weights * damod$residuals, df$cid[!is.na(df$cid)],
              colSums))
  )
  expect_equal(b22, .get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(Z * damod$weights, Z))
  expect_equal(a22inv, .get_a22_inverse(damod))
  a21 <- crossprod(Z, X * damod$weights)
  expect_equal(a21, .get_a21(damod))
  b12 <- matrix(0, nrow = dim(Xstar)[2], ncol = dim(Z)[2])
  expect_equal(b12, suppressMessages(.get_b12(damod)))
  b11 <- crossprod(Xstar * cmod$residuals^2, Xstar)
  expect_equal(b11, .get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(NQ * MI * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meatCL
  expect_equal(b22 / (NQ * MI),
               sandwich::meatCL(damod, cluster = factor(df$cid[!is.na(df$cid)]),
                                cadjust = FALSE))
  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  ## HEURISTIC CHECKS
  # with perfect group balance, a22inv %*% a21 should have 0's in the trt row
  expect_equal((a22inv %*% a21)[2,],
               setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_message(vmat <- .vcovMB_CR0(damod, cadjust = FALSE), "0 rows")
  expect_equal(vmat,
               a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  # var_hat(z) should be equal to that given by sandwich
  expect_equal(diag(vmat)[2],
               diag(sandwich::sandwich(damod,
                                       meat. = sandwich::meatCL,
                                       cluster = factor(df$cid[!is.na(df$cid)]),
                                       cadjust = FALSE))[2])
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


  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  a11inv <- solve(crossprod(Xstar))
  expect_equal(a11inv, .get_a11_inverse(damod))
  b22 <- crossprod(
    Reduce(rbind,
           by(Z * damod$weights * damod$residuals, df$cid[!is.na(df$cid)],
              colSums))
  )
  expect_equal(b22, .get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(Z * damod$weights, Z))
  expect_equal(a22inv, .get_a22_inverse(damod))
  a21 <- crossprod(Z, X * damod$weights)
  expect_equal(a21, .get_a21(damod))
  b12 <- matrix(0, nrow = dim(Xstar)[2], ncol = dim(Z)[2])
  expect_equal(b12, suppressMessages(.get_b12(damod)))
  b11 <- crossprod(Xstar * cmod$residuals^2, Xstar)
  expect_equal(b11, .get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(NQ * MI * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meatCL
  expect_equal(b22 / (NQ * MI),
               sandwich::meatCL(damod, cluster = factor(df$cid[!is.na(df$cid)]),
                                cadjust = FALSE))
  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  ## HEURISTIC CHECKS
  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_message(vmat <- .vcovMB_CR0(damod, cadjust = FALSE), "0 rows")
  expect_equal(vmat,
               a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  # var_hat(z) should be greater than that given by sandwich
  expect_true(diag(vmat)[2] >
               diag(sandwich::sandwich(
                 damod,
                 meat. = sandwich::meatCL,
                 cluster = factor(df$cid[!is.na(df$cid)]),
                 cadjust = FALSE))[2])
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

  a11inv <- solve(crossprod(Xstar))
  expect_equal(a11inv, flexida:::.get_a11_inverse(damod))
  b22 <- crossprod(Z, diag((damod$weights * damod$residuals)^2)) %*% Z
  expect_equal(b22, flexida:::.get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(Z * damod$weights, Z))
  expect_equal(a22inv, flexida:::.get_a22_inverse(damod))
  a21 <- crossprod(Z, X * damod$weights)
  expect_equal(a21, flexida:::.get_a21(damod))
  b12 <- crossprod(
    Reduce(rbind, by(Xstar * cmod$residuals, df$uid[cmod_idx], colSums)),
    Reduce(rbind,
           by(Z[cmod_idx,] * (damod$residuals * damod$weights)[cmod_idx],
              df$uid[cmod_idx],
              colSums))
  )
  expect_equal(b12, suppressMessages(flexida:::.get_b12(damod)))
  b11 <- crossprod(Xstar * cmod$residuals^2, Xstar)
  expect_equal(b11, flexida:::.get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(N * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meat
  expect_equal(sandwich::meat(damod, adjust = FALSE),
               sandwich::meatCL(damod, cadjust = FALSE, type = "HC0"))
  expect_equal(b22 / N, sandwich::meat(damod, adjust = FALSE))

  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  ## HEURISTIC CHECKS
  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_message(vmat <- .vcovMB_CR0(damod, cadjust = FALSE),
                 paste(nrow(df[df$z==0,]), "rows"))
  expect_equal(vmat,
               a22inv %*%
                 (b22 - a21 %*% a11inv %*% b12 - t(b12) %*% a11inv %*% t(a21) +
                    (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  # var_hat(z) should be greater than that given by sandwich
  expect_true(diag(vmat)[2] >
                diag(sandwich::sandwich(damod, adjust = FALSE))[2])
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
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df[!is.na(df$cid),])

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  a11inv <- solve(crossprod(Xstar))
  expect_equal(a11inv, .get_a11_inverse(damod))
  b22 <- crossprod(
    Reduce(rbind,
           by(Z * damod$weights * damod$residuals, df$cid[!is.na(df$cid)],
              colSums))
  )
  expect_equal(b22, .get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(Z * damod$weights, Z))
  expect_equal(a22inv, .get_a22_inverse(damod))
  a21 <- crossprod(Z, X * damod$weights)
  expect_equal(a21, .get_a21(damod))
  b12 <- matrix(0, nrow = dim(Xstar)[2], ncol = dim(Z)[2])
  expect_equal(b12, suppressMessages(.get_b12(damod)))
  b11 <- crossprod(Xstar * cmod$residuals^2, Xstar)
  expect_equal(b11, .get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(NQ * MI * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meatCL
  expect_equal(b22 / (NQ * MI),
               sandwich::meatCL(damod, cluster = factor(df$cid[!is.na(df$cid)]),
                                cadjust = FALSE))
  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  ## HEURISTIC CHECKS
  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_message(vmat <- .vcovMB_CR0(damod, cadjust = FALSE), "0 rows")
  expect_equal(vmat,
               a22inv %*%
                 (b22 + (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  # var_hat(z) should be greater than that given by sandwich
  expect_true(diag(vmat)[2] >
                diag(sandwich::sandwich(
                  damod,
                  meat. = sandwich::meatCL,
                  cluster = factor(df$cid[!is.na(df$cid)]),
                  cadjust = FALSE))[2])
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

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  a11inv <- solve(crossprod(Xstar))
  expect_equal(a11inv, flexida:::.get_a11_inverse(damod))
  b22 <- crossprod(
    Reduce(rbind,
           by(Z * damod$weights * damod$residuals, df$cid, colSums))
  )
  expect_equal(b22, flexida:::.get_b22(damod, cadjust = FALSE, type = "HC0"))
  a22inv <- solve(crossprod(Z * damod$weights, Z))
  expect_equal(a22inv, flexida:::.get_a22_inverse(damod))
  a21 <- crossprod(Z, X * damod$weights)
  expect_equal(a21, flexida:::.get_a21(damod))
  b12 <- crossprod(
    Reduce(rbind, by(Xstar * cmod$residuals, df$cid[cmod_idx], colSums)),
    Reduce(rbind,
           by(Z[cmod_idx,] * (damod$residuals * damod$weights)[cmod_idx],
              df$cid[cmod_idx],
              colSums))
  )
  expect_equal(b12, suppressMessages(flexida:::.get_b12(damod)))
  b11 <- crossprod(
    Reduce(rbind, by(Xstar * cmod$residuals, df$cid[cmod_idx], colSums))
  )
  expect_equal(b11, flexida:::.get_b11(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal((NC + NQ) * MI * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meatCL
  expect_equal(b22 / ((NC + NQ) * MI),
               sandwich::meatCL(damod, cluster = factor(df$cid),
                                cadjust = FALSE))
  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv,
               sandwich::sandwich(cmod,
                                  meat. = sandwich::meatCL,
                                  cluster = factor(df$cid[cmod_idx]),
                                  cadjust = FALSE,
                                  type = "HC0"))

  ## HEURISTIC CHECKS
  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_message(vmat <- .vcovMB_CR0(damod, cadjust = FALSE),
                 paste(nrow(df[df$z == 0,]), "rows"))
  expect_equal(vmat,
               a22inv %*%
                 (b22 - a21 %*% a11inv %*% b12 - t(b12) %*% a11inv %*% t(a21) +
                    (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)

  # check that, since (a22inv %*% a21)["z",] != 0 and b21 != 0, flexida's
  # var_hat(z) should be greater than that given by sandwich
  expect_true(diag(vmat)[2] >
                diag(sandwich::sandwich(damod, adjust = FALSE))[2])
})


test_that(paste("HC0 .vcovMB_CR0 binomial glm",
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

  a11inv <- solve(crossprod(Xstar * sqrt(fit_wstar)))
  expect_equal(a11inv, .get_a11_inverse(damod))
  a21 <- crossprod(Z, X * w * mu_eta)
  expect_equal(a21, .get_a21(damod))
  a22inv <- solve(crossprod(Z * w, Z))
  expect_equal(a22inv, .get_a22_inverse(damod))
  b11 <- crossprod(Xstar * fit_wstar * rstar)
  expect_equal(b11, .get_b11(damod, cadjust = FALSE, type = "HC0"))
  b12 <- crossprod(Xstar * rstar * fit_wstar, (Z * w * r)[df$z == 0,])
  expect_equal(b12, suppressMessages(.get_b12(damod)))
  b22 <- crossprod(Z * w * r)
  expect_equal(b22, .get_b22(damod, cadjust = FALSE, type = "HC0"))

  ## COMPARE OUTPUTS TO SANDWICH
  # after scaling, a22inv is equivalent to sandwich::bread
  expect_equal(N * a22inv, sandwich::bread(damod))

  # after scaling, b22 is equivalent to sandwich::meat
  expect_equal(b22 / N, sandwich::meat(damod))

  # check the cmod sandwich
  expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

  ## HEURISTIC CHECKS
  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((a22inv %*% a21)[2,] == 0))

  # check .vcovMB_CR0 matches manual matrix multiplication
  expect_message(vmat <- .vcovMB_CR0(damod, cadjust = FALSE),
                 paste(nrow(df[df$z == 0,]), "rows"))
  expect_equal(vmat,
               a22inv %*%
                 (b22 - a21 %*% a11inv %*% b12 - t(b12) %*% a11inv %*% t(a21) +
                    (a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))) %*%
                 a22inv)
})
