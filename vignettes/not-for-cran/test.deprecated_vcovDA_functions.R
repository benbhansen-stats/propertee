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
                 "are found to have NA's")
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
                 "have NA's for some but not all")
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
                  
                  expect_equal(suppressMessages(.get_b12(m)),
                               t(cmod_eqns) %*% m_eqns)
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
