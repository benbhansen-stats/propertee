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
    by((m$weights * m$residuals * model.matrix(m))[msk, ],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))

  expect_equal(get_overlap_vcov_matrix(m),
               t(cmod_eqns) %*% m_eqns)
}

test_that("get_overlap_vcov_matrix used with DA model without SandwichLayer offset", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  offset <- stats::predict(cmod, simdata)
  m <- as.DirectAdjusted(
    lm(y ~ z, data = simdata, weights = ate(des), offset = offset)
  )

  expect_error(get_overlap_vcov_matrix(m),
               "must have an offset of class")
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for individual-level",
                "experimental data identical to cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)

  m <- as.DirectAdjusted(
   lm(y ~ z, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  test_correct_b12(m)
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for cluster-level",
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

  m <- as.DirectAdjusted(
   lm(y ~ z, data = Q_cluster, weights = ate(des), offset = cov_adj(cmod))
  )

  test_correct_b12(m)
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for individual-level",
                "experimental data that is a subset of cov model data)"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata, subset = simdata$cid2 == 1)
  weighted_design <- ate(des, data = simdata[simdata$cid2 == 1,])
  m <- as.DirectAdjusted(
   lm(y ~ z, data = simdata[simdata$cid2 == 1,],
      weights = weighted_design, offset = cov_adj(cmod))
  )

  test_correct_b12(m)
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for cluster-level",
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
  
  m <- as.DirectAdjusted(
    lm(y ~ z, data = Q_cluster_subset,
       weights = weighted_design, offset = cov_adj(cmod))
  )
  
  test_correct_b12(m)
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for experimental",
                "data that has no overlap with cov model data"), {
  data(simdata)
  des <- rct_design(z ~ cluster(cid1, cid2), data = simdata)
  C_no_cluster_ids <- simdata
  C_no_cluster_ids[, var_names(des, "u")] <- NA
  cmod <- lm(y ~ x, data = C_no_cluster_ids)

  m <- as.DirectAdjusted(
    lm(y ~ z + force, data = simdata, weights = ate(des), offset = cov_adj(cmod))
  )

  expect_equal(dim(get_overlap_vcov_matrix(m)), c(2, 3))
  expect_true(all(get_overlap_vcov_matrix(m) == 0))
})
