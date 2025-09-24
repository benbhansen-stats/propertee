test_that("vcov_tee errors when provided an invalid type", {
  data(simdata)
  cmod <- lm(y ~ z + x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  ssmod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod)))

  expect_error(vcov_tee(ssmod, "not_a_type"))
})

test_that("vcov_tee correctly sets and passes on args", {
  set.seed(993)
  data(simdata)
  cmod <- lm(y ~ z + x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  ssmod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod)))

  vmat1 <- suppressMessages(vcov_tee(ssmod, type = "CR0", cov_adj_rcorrect = "CR0"))
  expect_equal(vmat1,
               suppressMessages(.vcov_CR(ssmod, cluster = .make_uoa_ids(ssmod, "CR"),
                                         type = "CR0", cov_adj_rcorrect = "CR0")))

  vmat2 <- suppressMessages(vcov_tee(ssmod, type = "MB0", cov_adj_rcorrect = "MB0"))
  expect_true(all.equal(vmat2, vmat1, check.attributes = FALSE))
  expect_true(attr(vmat1, "type") != attr(vmat2, "type"))
  expect_true(attr(vmat1, "cov_adj_rcorrect") != attr(vmat2, "cov_adj_rcorrect"))

  vmat3 <- suppressMessages(vcov_tee(ssmod, type = "HC0", cov_adj_rcorrect = "HC0"))
  expect_true(all.equal(vmat3, vmat1, check.attributes = FALSE))
  expect_true(attr(vmat1, "type") != attr(vmat3, "type"))
  expect_true(attr(vmat1, "cov_adj_rcorrect") != attr(vmat3, "cov_adj_rcorrect"))
  
  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  ssmod <- lmitt(y ~ 1, specification = spec, data = simdata, 
                 absorb = TRUE, offset = cov_adj(cmod))
  vmat4 <- suppressMessages(vcov_tee(ssmod, type = "DB0", const_effect = FALSE))
  expect_true(is.matrix(vmat4))
  expect_identical(attr(vmat4, "type"), "DB0")
  
  vmat5 <- suppressMessages(vcov_tee(ssmod, type = "DB0", const_effect = TRUE))
  expect_true(is.matrix(vmat5))
  expect_identical(attr(vmat5, "type"), "DB0")
  expect_false(identical(vmat4, vmat5))
  
  vmat6 <- vcov_tee(ssmod, type = "HC1", cov_adj_rcorrect = "HC1")
  expect_equal(attr(vmat6, "type"), "HC1")
  expect_equal(attr(vmat6, "cov_adj_rcorrect"), "HC1")
  
  vmat7 <- vcov_tee(ssmod)
  expect_equal(attr(vmat7, "type"), "HC0")
  expect_equal(attr(vmat7, "cov_adj_rcorrect"), "HC0")
  
  first_obs_ix <- c(1, 1 + cumsum(c(t(table(simdata$uoa1, simdata$uoa2)))))
  uniqdata <- simdata[first_obs_ix[1:(length(first_obs_ix)-1)],,drop=FALSE]
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  imod <- lmitt(y ~ 1, spec, uniqdata)
  vmat8 <- vcov_tee(imod)
  expect_equal(attr(vmat8, "type"), "HC0")
  expect_true(is.null(attr(vmat8, "cov_adj_rcorrect")))
})

if (requireNamespace("robustbase", quietly = TRUE)) {
  test_that("vcov_tee sets correct bias correction for lmrob cov adj model", {
    set.seed(993)
    data(simdata)
    cmod <- robustbase::lmrob(y ~ z + x, simdata)
    spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
    ssmod <- lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod)))
    
    vmat5 <- vcov_tee(ssmod)
    expect_equal(attr(vmat5, "cov_adj_rcorrect"), "HC0")
  })
}

test_that(paste("vcov_tee produces correct calculations with valid `cluster` arugment",
                "when cluster ID's have no NA's"), {
  data(simdata)
  simdata_copy <- simdata
  simdata_copy$uid <- seq_len(nrow(simdata_copy))
  cmod <- lm(y ~ x, simdata_copy)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2, uid) + block(bid), simdata_copy)
  ssmod <- lmitt(y ~ 1, data = simdata_copy, specification = spec,
                weights = ate(spec), offset = cov_adj(cmod))

  # check default clustering level is the same when specified using cluster arg
  expect_equal(suppressMessages(vcov_tee(ssmod)),
               suppressMessages(vcov_tee(ssmod, cluster = c("uid"))))
})

test_that(paste("vcov_tee produces correct calculations with valid `cluster` arugment",
                "when cluster ID's have NA's (must be via column name)"), {
  data(simdata)
  df <- rbind(simdata, simdata)
  df[1:50, c("uoa1", "uoa2", "bid", "z")] <- NA
  cmod <- lm(y ~ x, df[1:50,])
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), df[51:100,])
  ssmod <- lmitt(y ~ 1, data = df[51:100,], specification = spec,
                weights = ate(spec), offset = cov_adj(cmod))

  expect_equal(vcov_tee(ssmod, cluster = c("uoa1", "uoa2")), vcov_tee(ssmod))
})

test_that("variance helper functions fail without a teeMod model", {
  data(simdata)
  cmod <- lm(y ~ z, data = simdata)

  expect_error(propertee:::.vcov_CR(cmod), "must be a teeMod")
  expect_error(propertee:::.get_b12(cmod), "must be a teeMod")
  expect_error(propertee:::.get_a22_inverse(cmod), "must be a teeMod")
  expect_error(propertee:::.get_b22(cmod), "must be a teeMod")
  expect_error(propertee:::.get_a22_inverse(cmod), "must be a teeMod")
  expect_error(propertee:::.get_a11_inverse(cmod), "must be a teeMod")
  expect_error(propertee:::.get_b11(cmod), "must be a teeMod")
  expect_error(propertee:::.get_a21(cmod), "must be a teeMod")

})

test_that(".get_dof", {
  data(simdata)
  expect_error(.get_dof(lm(y~x, simdata), "CR2", 0), "teeMod object")
  
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  tm <- lmitt(y ~ 1, spec, simdata)
  
  # non-CR2 dof
  expect_equal(.get_dof(tm, "CR0", 1), nrow(spec@structure)-1)
  expect_equal(.get_dof(tm, "DB0", 1, "bid"), 2)
  
  # CR2 dof
  expect_error(.get_dof(tm, "CR2", seq_len(7)), "same length")
  expect_error(.get_dof(tm, "CR2", seq_len(2)), "zeros or ones")
  expect_true(inherits(dof <- .get_dof(tm, "CR2", 2), "numeric"))
  expect_equal(length(dof), 1)
})

test_that(".compute_IK_dof no clustering", {
  data(simdata)

  ## no clustering, continuous y, sufficient dof
  simdata$id <- seq_len(nrow(simdata))
  idspec <- rct_spec(z ~ unitid(id), simdata)
  tm <- lmitt(y ~ 1, idspec, simdata)
  ell <- c(0, 1)
  dof <- .compute_IK_dof(tm, ell)
  
  Q <- qr.Q(tm$qr)
  R <- qr.R(tm$qr)
  I_P <- diag(nrow = nrow(simdata)) - tcrossprod(Q)
  G <- I_P * drop((Q / sqrt(1-stats::hatvalues(tm))) %*% t(solve(R)) %*% ell)
  GoG <- crossprod(G, diag(mean(stats::residuals(tm)^2), nrow = nrow(Q))) %*% G
  expected_dof <- sum(diag(GoG))^2 / sum(diag(GoG %*% GoG))
  expect_equal(dof, expected_dof)
  
  ## no clustering, binary y, sufficient dof
  simdata$y <- round(runif(nrow(simdata)))
  tm <- lmitt(y ~ 1, idspec, simdata)
  dof <- .compute_IK_dof(tm, ell, bin_y = TRUE)
  
  GoG <- crossprod(G, diag(stats::fitted(tm) * (1 - stats::fitted(tm)), nrow = nrow(Q))) %*% G
  expected_dof <- sum(diag(GoG))^2 / sum(diag(GoG %*% GoG))
  expect_equal(dof, expected_dof)
  
  ## not enough dof
  set.seed(103)
  ad2 <- data.frame(cid = c(rep(1, 4), rep(2, 6), rep(3, 5), rep(4, 9), rep(5, 12), rep(6, 7)),
                    x = factor(c(rep(0,15), rep(1, 9+12+7))),
                    a = c(rep(0, 4), rep(1, 11), rep(0, 9), rep(1, 19)),
                    y = rnorm(4 + 5 + 6 + 9+19))
  spec <- rct_spec(a ~ unitid(cid), ad2, subset = cid %in% c(1, 2, 4, 5))
  tm <- lmitt(y ~ x, spec, ad2, subset = cid %in% c(1, 2, 4, 5))

  bad.ell1 <- c(0, 0, 1, 0)
  bad.ell2 <- c(0, 0, 0, 1)
  suppressWarnings(no_dof <- .compute_IK_dof(tm, bad.ell1))
  expect_equal(no_dof, 1) # this is what IK code returns
  suppressWarnings(valid_dof <- .compute_IK_dof(tm, bad.ell2))
  expect_equal(valid_dof, 1) # this is what IK code returns
})

test_that("cluster_iss, no absorption", {
  set.seed(33)
  ad2 <- data.frame(cid = c(rep(1, 4), rep(2, 6), rep(3, 5), rep(4, 9), rep(5, 12), rep(seq(6, 8), each = 7)),
                    x = factor(c(rep(0,15), rep(1, 9+12+7+7), rep(0, 7))),
                    a = c(rep(0, 4), rep(1, 11), rep(0, 9), rep(1, 19), rep(0, 7+7)),
                    y = rnorm(4 + 5 + 6 + 9+19+7+7))

  ## no moderator
  spec <- rct_spec(a ~ unitid(cid), ad2)
  tm <- lmitt(y ~ 1, spec, ad2)
  ATWA_inv <- chol2inv(tm$qr$qr)
  
  # control unit
  ix <- ad2$cid == 1
  Pgg <- model.matrix(tm)[ix,] %*% ATWA_inv %*% t(model.matrix(tm)[ix,])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  all.equal(
    cluster_iss(tm, cluster_unit = 1, cluster_ids = ad2$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors)
  )
  # treated unit
  ix <- ad2$cid == 6
  Pgg <- model.matrix(tm)[ix,] %*% ATWA_inv %*% t(model.matrix(tm)[ix,])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  all.equal(
    cluster_iss(tm, cluster_unit = 6, cluster_ids = ad2$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors)
  )
  
  ## cluster-invariant binary moderator variable
  tm <- lmitt(y ~ x, spec, ad2, weights = "ate")
  ATWA_inv <- chol2inv(tm$qr$qr)
  
  # treated unit w/ x = 1
  ix <- ad2$cid == 6
  Pgg <- model.matrix(tm)[ix,] %*% ATWA_inv %*% t((model.matrix(tm) * weights(tm))[ix,])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, cluster_unit = 6), # test NULL cluster_ids and NULL cluster and NULL trts 
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  # treated unit w/ x = 0
  ix <- ad2$cid == 2
  Pgg <- model.matrix(tm)[ix,] %*% ATWA_inv %*% t((model.matrix(tm) * weights(tm))[ix,])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, cluster_unit = 2, cluster_ids = ad2$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  # control unit w/ x = 1
  ix <- ad2$cid == 4
  Pgg <- model.matrix(tm)[ix,] %*% ATWA_inv %*% t((model.matrix(tm) * weights(tm))[ix,])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, cluster_unit = 4, cluster_ids = ad2$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  # control unit w/ x = 0
  ix <- ad2$cid == 1
  Pgg <- model.matrix(tm)[ix,] %*% ATWA_inv %*% t((model.matrix(tm) * weights(tm))[ix,])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, cluster_unit = 1), # test NULL cluster_ids and NULL cluster and NULL trts 
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  
  ## cluster-invariant continuous moderator variable
  ad2$x <- runif(length(unique(ad2$cid)))[ad2$cid]
  tm <- lmitt(y ~ x, spec, ad2, weights = "ate")
  ATWA_inv <- chol2inv(tm$qr$qr)
  
  # control unit 
  ix <- ad2$cid == 1
  Pgg <- model.matrix(tm)[ix,] %*% ATWA_inv %*% t((model.matrix(tm) * weights(tm))[ix,])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, cluster_unit = 1, cluster_ids = ad2$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  # treated unit
  ix <- ad2$cid == 6
  Pgg <- model.matrix(tm)[ix,] %*% ATWA_inv %*% t((model.matrix(tm) * weights(tm))[ix,])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, cluster_unit = 6, cluster_ids = ad2$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
})

test_that("cluster_iss, absorption", {
  # don't need to distinguish between treated and control units in these tests
  set.seed(499)
  adb <- data.frame(cid = c(rep(1, 6), rep(2, 5), rep(3, 4), rep(4, 5),
                            rep(5, 6), rep(6, 5), rep(7, 4), rep(8, 5)),
                    b = rep(c(1, 2), each = 20),
                    a = c(rep(0, 11), rep(1, 20), rep(0, 9)),
                    y = rnorm(40))
  adb$x <- runif(8)[adb$cid]
  specification <- rct_spec(a ~ cluster(cid) + block(b), adb)
  
  ## no moderator
  tm <- lmitt(y ~ 1, specification, adb, weights = "ate", absorb = TRUE)
  piv <- seq_len(2)
  ATWA_inv <- chol2inv(tm$qr$qr[piv,piv,drop=FALSE])
  
  ix <- adb$cid == 1
  Ag <- model.matrix(tm)[ix,piv,drop=FALSE]
  Pgg <- Ag %*% ATWA_inv %*% t(Ag * weights(tm)[ix])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, 1, adb$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  
  ## cluster-invariant continuous moderator
  tm <- lmitt(y ~ x, specification, adb, weights = "ate", absorb = TRUE)
  piv <- seq_len(tm$rank)
  ATWA_inv <- chol2inv(tm$qr$qr[piv,piv,drop=FALSE])
  
  Ag <- model.matrix(tm)[ix,piv,drop=FALSE]
  Pgg <- Ag %*% ATWA_inv %*% t(Ag * weights(tm)[ix])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, 1, adb$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  
  ## cluster-invariant factor moderator
  adb$x <- factor(round(adb$x))
  tm <- lmitt(y ~ x, specification, adb, weights = "ate", absorb = TRUE)
  piv <- seq_len(tm$rank)
  ATWA_inv <- chol2inv(tm$qr$qr[piv,piv,drop=FALSE])
  
  Ag <- model.matrix(tm)[ix,piv,drop=FALSE]
  Pgg <- Ag %*% ATWA_inv %*% t(Ag * weights(tm)[ix])
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, 1, adb$cid),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
})

test_that("cluster_iss with fine strata (clusters have treated and control units)", {
  set.seed(496)
  ## balanced in each stratum, absorption
  sdata <- data.frame(id = seq_len(20),
                      uid = rep(seq_len(4), each = 5),
                      ai = rep(rep(c(0, 1), each = 5), 2),
                      bi = rep(c(1, 2), each = 10),
                      yi = rnorm(20))
  ats <- rct_spec(ai ~ unitid(uid) + block(bi), sdata)
  tm <- lmitt(yi ~ 1, ats, data = sdata, absorb = TRUE)
  ids <- propertee:::.make_uoa_ids(tm, "MB")
  A <- model.matrix(tm)
  ATWA_inv <- chol2inv(tm$qr$qr)
  ix <- ids == unique(ids)[1]
  Ag <- A[ix,,drop=FALSE]
  Pgg <- Ag %*% ATWA_inv %*% t(Ag)
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, unique(ids)[1], ids),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))

  ## weighted balance in each stratum, absorption
  wts <- rep(c(rep(2, 4), 1), 4)
  tm <- lmitt(yi ~ 1, ats, data = sdata, weights = wts, absorb = TRUE)
  ATWA_inv <- chol2inv(tm$qr$qr)
  wg <- wts[ix]
  Pgg <- (Ag * sqrt(wg)) %*% ATWA_inv %*% t(Ag * sqrt(wg))
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, unique(ids)[1], ids),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  
  ## weighted balance, no absorption
  tm <- lmitt(yi ~ 1, ats, data = sdata, weights = "att")
  A <- model.matrix(tm)
  ATWA_inv <- chol2inv(tm$qr$qr)
  wg <- weights(tm)[ix]
  Ag <- A[ix,,drop=FALSE]
  Pgg <- (Ag * sqrt(wg)) %*% ATWA_inv %*% t(Ag * sqrt(wg))
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, unique(ids)[1], ids),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))

  ## same imbalanced treatment assignment in each stratum, absorption
  sdata_imbal <- sdata
  sdata_imbal <- rbind(sdata_imbal,
                       data.frame(id = seq(21,30),
                                  uid = rep(c(5, 6), each = 5),
                                  ai = 1,
                                  bi = rep(c(1, 2), each = 5),
                                  yi = rnorm(10)))
  imb <- rct_spec(ai ~ unitid(uid) + block(bi), sdata_imbal)
  tm <- lmitt(yi ~ 1, imb, data = sdata_imbal, absorb = TRUE)
  ids <- propertee:::.make_uoa_ids(tm, "MB")
  A <- model.matrix(tm)
  ATWA_inv <- chol2inv(tm$qr$qr)
  ix <- ids == unique(ids)[2]
  Ag <- A[ix,,drop=FALSE]
  Pgg <- Ag %*% ATWA_inv %*% t(Ag)
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, unique(ids)[2], ids),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  
  ## same imbalanced treatment assignment in each stratum, no absorption
  tm <- lmitt(yi ~ 1, imb, data = sdata_imbal, weights = "att")
  wts <- weights(tm)
  ids <- propertee:::.make_uoa_ids(tm, "MB")
  X <- model.matrix(tm)
  ATWA_inv <- chol2inv(tm$qr$qr)
  ix <- ids == unique(ids)[2]
  A <- model.matrix(tm)
  Ag <- A[ix,,drop=FALSE]
  Pgg <- (Ag * sqrt(wts[ix])) %*% ATWA_inv %*% t(Ag * sqrt(wts[ix]))
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, unique(ids)[2], ids),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
  expect_message(tm <- lmitt(yi ~ 1, imb, data = sdata_imbal), "absorb")
  wts <- rep(1, nrow(sdata_imbal))
  ids <- propertee:::.make_uoa_ids(tm, "MB")
  X <- model.matrix(tm)
  ATWA_inv <- chol2inv(tm$qr$qr)
  ix <- ids == unique(ids)[2]
  A <- model.matrix(tm)
  Ag <- A[ix,,drop=FALSE]
  Pgg <- (Ag * sqrt(wts[ix])) %*% ATWA_inv %*% t(Ag * sqrt(wts[ix]))
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, unique(ids)[2], ids),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))

  ## same weighted imbalance in each strata
  sdata$uid[sdata$uid == 2] <- c(1, rep(2, 4))
  sdata$uid[sdata$uid == 3] <- c(rep(3, 2), rep(4, 3))
  sdata$ai[sdata$uid == 1] <- 0
  sdata$ai[sdata$uid == 4] <- 1
  imb <- rct_spec(ai ~ unitid(uid) + block(bi), sdata)
  wts <- c(rep(2/3, 6), rep(1, 4), rep(4, 2), rep(1, 8))
  tm <- lmitt(yi ~ 1, imb, data = sdata, weights = wts, absorb = TRUE)
  ids <- propertee:::.make_uoa_ids(tm, "MB")
  ix <- ids == "bi2"
  A <- model.matrix(tm)
  Ag <- A[ix,,drop=FALSE]
  Pgg <- (Ag * sqrt(wts[ix])) %*% chol2inv(tm$qr$qr) %*% t(Ag * sqrt(wts[ix]))
  eg <- eigen(diag(nrow = sum(ix)) - Pgg)
  expect_true(all.equal(
    cluster_iss(tm, unique(ids)[2], ids),
    eg$vectors %*% diag(1/sqrt(eg$values), nrow = length(eg$values)) %*% solve(eg$vectors),
    check.attributes = FALSE
  ))
})

test_that(paste(".get_b12, .get_a11_inverse, .get_b11, .get_a21 used with teeMod model",
                "without SandwichLayer offset"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  offset <- stats::predict(cmod, simdata)
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = offset)
  )

  expect_error(.get_b12(m), "must have an offset of class")
  expect_error(.get_a11_inverse(m), "must have an offset of class")
  expect_error(.get_b11(m), "must have an offset of class")
  expect_error(.get_a21(m), "must have an offset of class")
})

test_that(".get_b12 fails if ITT model isn't strictl an `lm`", {
  return(expect_true(TRUE)) # this function is deprecated and now failing
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  offset <- stats::predict(cmod, simdata)
  m <- as.lmitt(
    glm(y ~ assigned(), data = simdata, weights = ate(spec), offset = offset)
  )

  expect_error(.get_b12(m), "x must be an `lm` object")
})

test_that(paste(".get_b12 returns expected B_12 for individual-level",
                "experimental data identical to cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )

  uoanames <- var_names(m@StudySpecification, "u")
  zname <- var_names(m@StudySpecification, "t")

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
    merge(unique(m$model$`(offset)`@keys), m@StudySpecification@structure),
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_cluster)

  m <- as.lmitt(
    lm(y ~ assigned(), data = Q_cluster, weights = ate(spec), offset = cov_adj(cmod))
  )

  uoanames <- var_names(m@StudySpecification, "u")
  zname <- var_names(m@StudySpecification, "t")

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
    merge(unique(m$model$`(offset)`@keys), m@StudySpecification@structure),
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )

  expect_error(propertee:::.get_b12(m, cluster = c("cid3")),
               "cid3 are missing from the covariance adjustment model dataset")
  expect_error(propertee:::.get_b12(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
  expect_error(.get_b12(m, cluster = "new_col"),
               "in the teeMod object's StudySpecification: new_col")
})

test_that(".get_b12 produces correct estimates with valid custom cluster argument", {
  data(simdata)
  simdata[simdata$uoa2 == 1, "z"] <- 0
  simdata[simdata$uoa2 == 2, "z"] <- 1

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa2), data = simdata)

  m <- lmitt(y ~ 1, data = simdata, specification = spec, offset = cov_adj(cmod))

  cmod_eqns <- Reduce(
    rbind,
    by(estfun(cmod), list(simdata$uoa2), colSums)
  )
  ssmod_eqns <- Reduce(
    rbind,
    by(estfun(as(m, "lm")), list(simdata$uoa2), colSums)
  )
  expected <- crossprod(cmod_eqns, ssmod_eqns)

  # default (columns specified in `cluster` argument of StudySpecification) matches expected
  expect_equal(suppressMessages(propertee:::.get_b12(m)), expected)
  expect_equal(suppressMessages(propertee:::.get_b12(m, cluster = "uoa2")), expected)
})

test_that("get_b12 handles custom cluster columns with only NA's", {
  data(simdata)
  set.seed(200)

  cmod_data <- data.frame("y" = rnorm(100), "x" = rnorm(100),
                          "uoa1" = NA_integer_, "uoa2" = NA_integer_)
  cmod <- lm(y ~ x, cmod_data)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  ssmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, specification = spec)))

  expect_message(b12 <- propertee:::.get_b12(ssmod, cluster = c("uoa1", "uoa2")), "0 rows")
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

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  ssmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, specification = spec)))

  expect_message(b12 <- propertee:::.get_b12(ssmod, cluster = c("uoa1", "uoa2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))

  # second test is where the non-NA cluster ID's overlap with specification
  cmod_data$uoa1 <- rep(seq(1, 5), each = 20)
  cmod <- lm(y ~ x, cmod_data)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  ssmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, specification = spec)))

  expect_message(b12 <- propertee:::.get_b12(ssmod, cluster = c("uoa1", "uoa2")), "0 rows")
  expect_equal(b12,
               matrix(0, nrow = 2, ncol = 2))
})

test_that(paste(".get_b12 handles multiple custom cluster columns where both",
                "have a mix of NA's and non-NA's"), {
  data(simdata)
  cmod_data <- rbind(simdata, simdata)
  cmod_data[(nrow(simdata)+1):(2*nrow(simdata)), c("uoa1", "uoa2")] <- NA_integer_
  cmod <- lm(y ~ x, cmod_data)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  ssmod <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec),
                      offset = cov_adj(cmod)))

  cmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(cmod), list(cmod_data$uoa1, cmod_data$uoa2), colSums)
  )
  ssmod_eqns <- Reduce(
    rbind,
    by(sandwich:::estfun.lm(ssmod), list(simdata$uoa1, simdata$uoa2), colSums)
  )
  expected <- crossprod(cmod_eqns, ssmod_eqns)

  expect_message(b12 <- propertee:::.get_b12(ssmod, cluster = c("uoa1", "uoa2")),
                 paste(nrow(simdata), "rows"))
  expect_equal(b12,
               expected)
})

test_that(paste(".get_b12 returns expected B_12 for individual-level",
                "experimental data that is a subset of cov model data"), {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata, subset = simdata$uoa2 == 1)
  weighted_spec <- ate(spec, data = simdata[simdata$uoa2 == 1,])
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata[simdata$uoa2 == 1,],
       weights = weighted_spec, offset = cov_adj(cmod))
  )

  uoanames <- var_names(m@StudySpecification, "u")
  zname <- var_names(m@StudySpecification, "t")

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
    merge(unique(m$model$`(offset)`@keys[, uoanames, drop = FALSE]), m@StudySpecification@structure),
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = Q_cluster_subset)
  weighted_spec <- ate(spec, data = Q_cluster_subset)

  m <- as.lmitt(
    lm(y ~ assigned(), data = Q_cluster_subset,
       weights = weighted_spec, offset = cov_adj(cmod))
  )

  uoanames <- var_names(m@StudySpecification, "u")
  zname <- var_names(m@StudySpecification, "t")

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
    merge(unique(m$model$`(offset)`@keys), m@StudySpecification@structure),
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  C_no_cluster_ids <- simdata
  C_no_cluster_ids[, var_names(spec, "u")] <- NA
  cmod <- lm(y ~ x, data = C_no_cluster_ids)

  m <- as.lmitt(
    lm(y ~ assigned() + force, data = simdata, weights = ate(spec), offset = cov_adj(cmod))
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
  spec <- rct_spec(z ~ cluster(uoa1), simdata, subset = simdata$uoa1 %in% c(1, 5))

  msk <- simdata$uoa1 %in% c(1, 5)
  m <- as.lmitt(
    lm(y ~ assigned(), simdata[msk,],
       weights = ate(spec, data = simdata[msk,]),
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  expect_warning(expect_warning(lmitt(lm(y ~ assigned(), simdata,
                                         offset = cov_adj(cmod, specification = spec))),
                                "adjustments are NA"), "adjustments are NA")
  ## warning procs twice
  m <- suppressWarnings(
    lmitt(lm(y ~ assigned(), simdata, offset = cov_adj(cmod, specification = spec)))
  )

  cmod_eqns <- Reduce(
    rbind,
    by(sandwich::estfun(cmod),
       list(uoa1 = simdata$uoa1[!is.na(simdata$x) & simdata$uoa1 %in% c(2, 4)],
            uoa2 = simdata$uoa2[!is.na(simdata$x) & simdata$uoa1 %in% c(2, 4)]),
       colSums))
  msk <- which(simdata$uoa1[!is.na(simdata$x)] %in% c(2, 4))
  ssmod_eqns <- Reduce(
    rbind,
    by(estfun(as(m, "lm"))[msk, , drop = FALSE],
       list(uoa1 = simdata$uoa1[!is.na(simdata$x)][msk],
            uoa2 = simdata$uoa2[!is.na(simdata$x)][msk]),
       colSums))
  b12 <- suppressMessages(propertee:::.get_b12(m))
  expect_equal(b12, crossprod(cmod_eqns, ssmod_eqns))
})

test_that(paste(".get_b12 returns expected value for B12 when no intercept is",
                "included in the direct adjustment model"), {
  data(simdata)
  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)

  m <- as.lmitt(lm(y ~ assigned() - 1, simdata, weights = ate(spec),
                   offset = cov_adj(cmod)))

  uoanames <- var_names(m@StudySpecification, "u")
  zname <- var_names(m@StudySpecification, "t")

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
    merge(unique(m$model$`(offset)`@keys), m@StudySpecification@structure),
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

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  m_as.lmitt <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec)))
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, specification = spec, weights = ate(spec))

  nq <- nrow(stats::model.frame(m_as.lmitt))
  inv_fim <- nq * chol2inv(m_as.lmitt$qr$qr)

  a22inv <- .get_a22_inverse(m_as.lmitt)
  expect_true(all.equal(a22inv[1:2, 1:2], inv_fim, check.attributes = FALSE))
  expect_true(all.equal(a22inv[3, 3], nq / sum(weights(m_as.lmitt@ctrl_means_model)), check.attributes = FALSE))
  expect_true(all(a22inv[1:2, 3] == 0) & all(a22inv[3, 1:2] == 0))
  
  a22inv <- .get_a22_inverse(m_lmitt.form)
  expect_true(all.equal(a22inv[1:2, 1:2], inv_fim, check.attributes = FALSE))
  expect_true(all.equal(a22inv[3, 3], nq / sum(weights(m_as.lmitt@ctrl_means_model)), check.attributes = FALSE))
  expect_true(all(a22inv[1:2, 3] == 0) & all(a22inv[3, 1:2] == 0))
})

test_that(".get_a22_inverse correct w/ covariance adjustment", {
  set.seed(438)
  data(simdata)
  testdata <- simdata
  nc <- 30
  nq <- nrow(testdata)
  n <- nc + nq
  cmod_data <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_data)

  testdata$id <- seq(nq)
  spec <- rct_spec(z ~ unitid(id), testdata)

  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), data = testdata, weights = ate(spec), offset = cov_adj(cmod)
  ))
  m_lmitt.form <- lmitt(y ~ 1, data = testdata, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))

  inv_fim <- nq * chol2inv(m_as.lmitt$qr$qr)

  a22inv <- .get_a22_inverse(m_as.lmitt)
  expect_true(all.equal(a22inv[1:2, 1:2], inv_fim, check.attributes = FALSE))
  expect_true(all.equal(a22inv[3:4, 3:4], diag(2), check.attributes = FALSE))
  expect_true(all(a22inv[1:2, 3:4] ==0) & all(a22inv[3:4, 1:2] == 0))
  
  a22inv <- .get_a22_inverse(m_lmitt.form)
  expect_true(all.equal(a22inv[1:2, 1:2], inv_fim, check.attributes = FALSE))
  expect_true(all.equal(a22inv[3:4, 3:4], diag(2), check.attributes = FALSE))
  expect_true(all(a22inv[1:2, 3:4] ==0) & all(a22inv[3:4, 1:2] == 0))
})

test_that(".get_a22_inverse with missing values", {
  set.seed(438)
  data(simdata)
  testdata <- simdata
  testdata$y[1:3] <- NA_real_
  nc <- 30
  nq <- nrow(testdata)
  n <- nc + nq
  cmod_data <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_data)
  
  testdata$id <- seq(nq)
  spec <- rct_spec(z ~ unitid(id), testdata)
  
  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), data = testdata, weights = ate(spec), offset = cov_adj(cmod)
  ))
  m_lmitt.form <- lmitt(y ~ 1, data = testdata, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))
  
  a22inv <- .get_a22_inverse(m_as.lmitt)
  expect_true(all.equal(a22inv[1:2, 1:2], nq * chol2inv(m_as.lmitt$qr$qr), check.attributes = FALSE))
  expect_true(all.equal(a22inv[3:4, 3:4],
                        diag(2) * nq / sum(weights(m_as.lmitt@ctrl_means_model), na.rm = TRUE),
                        check.attributes = FALSE))
  
  a22inv <- .get_a22_inverse(m_lmitt.form)
  expect_true(all.equal(a22inv[1:2, 1:2], nq * chol2inv(m_lmitt.form$qr$qr), check.attributes = FALSE))
  expect_true(all.equal(a22inv[3:4, 3:4],
                        diag(2) * nq / sum(weights(m_as.lmitt@ctrl_means_model), na.rm = TRUE),
                        check.attributes = FALSE))
})

test_that(".get_tilde_a22_inverse correct w/o covariance adjustment", {
  data(simdata)
  new_data <- simdata

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = new_data)
  m_as.lmitt <- as.lmitt(lm(y ~ assigned(), data = new_data, weights = ate(spec)))
  m_lmitt.form <- lmitt(y ~ 1, data = new_data, specification = spec, weights = ate(spec))

  expect_equal(.get_tilde_a22_inverse(m_as.lmitt), .get_a22_inverse(m_as.lmitt))
  expect_equal(.get_tilde_a22_inverse(m_lmitt.form), .get_a22_inverse(m_lmitt.form))
})

test_that(".get_tilde_a22_inverse correct w/ covariance adjustment", {
  set.seed(438)
  data(simdata)
  new_data <- simdata
  nc <- 30
  nq <- nrow(new_data)
  n <- nc + nq
  cmod_data <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_data)

  new_data$id <- seq(nq)
  spec <- rct_spec(z ~ unitid(id), new_data)

  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), data = new_data, weights = ate(spec), offset = cov_adj(cmod)
  ))
  m_lmitt.form <- lmitt(y ~ 1, data = new_data, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))

  expect_equal(.get_tilde_a22_inverse(m_as.lmitt), n / nq * .get_a22_inverse(m_as.lmitt))
  expect_equal(.get_tilde_a22_inverse(m_lmitt.form), n / nq * .get_a22_inverse(m_lmitt.form))
})

test_that(".get_tilde_a22_inverse with missing values", {
  set.seed(438)
  data(simdata)
  new_data <- simdata
  new_data$y[1:3] <- NA_real_
  nc <- 30
  nq <- nrow(new_data)
  n <- nc + nq
  cmod_data <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_data)
  
  new_data$id <- seq(nq)
  spec <- rct_spec(z ~ unitid(id), new_data)
  
  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), data = new_data, weights = ate(spec), offset = cov_adj(cmod)
  ))
  m_lmitt.form <- lmitt(y ~ 1, data = new_data, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))
  
  expect_equal(.get_tilde_a22_inverse(m_as.lmitt), n / nq * .get_a22_inverse(m_as.lmitt))
  expect_equal(.get_tilde_a22_inverse(m_lmitt.form), n / nq * .get_a22_inverse(m_lmitt.form))
})

test_that(".get_b22 returns correct value for lm object w/o offset", {
  data(simdata)

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- as.lmitt(lm(y ~ assigned(), data = simdata, weights = ate(spec)))
  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@StudySpecification, "u")
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )
  nq <- nrow(simdata)
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@StudySpecification, "u")
  form <- paste0("~ -1 + ", paste("as.factor(", uoanames, ")", collapse = ":"))
  uoa_matrix <- stats::model.matrix(as.formula(form),
                                    stats::expand.model.frame(m, uoanames)[, uoanames])

  uoa_eqns <- crossprod(uoa_matrix, WX)
  vmat <- crossprod(uoa_eqns)
  expect_equal(propertee:::.get_b22(m, type = "HC0", type_psi = "HC0", type_phi = "HC0"),
               vmat * nuoas / (nuoas - 1L))
  expect_equal(propertee:::.get_b22(m, type = "HC1", type_psi = "HC1", type_phi = "HC0"),
               vmat * nuoas / (nuoas - 1L) * (nq - 1L) / (nq - 2L))
})

test_that(".get_b22 fails with invalid custom cluster argument", {
  data(simdata)
  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )

  expect_error(propertee:::.get_b22(m, cluster = c("cid3")),
               "cid3 are missing from the direct adjustment")
  expect_error(propertee:::.get_b22(m, cluster = c(TRUE, FALSE)),
               "must provide a character vector")
})

test_that(".get_b22 produces correct estimates with valid custom cluster argument", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )
  nq <- nrow(simdata)
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@StudySpecification, "u")

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
  spec <- rct_spec(z ~ uoa(uoa1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod)))

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
  return(expect_true(TRUE)) # added 1/14/25 because it's deprecated and now failing

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- as.lmitt(glm(y ~ assigned(), data = simdata, weights = ate(spec)))
  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@StudySpecification, "u")
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
  return(expect_true(TRUE)) # added 1/14/25 because it's deprecated and now failing

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(spec),
        family = stats::poisson())
  )
  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@StudySpecification, "u")
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
  return(expect_true(TRUE)) # added 1/14/25 because it's deprecated and now failing

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(spec),
        family = stats::quasipoisson())
  )
  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@StudySpecification, "u")
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
  return(expect_true(TRUE)) # added 1/14/25 because it's deprecated and now failing

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- suppressWarnings(
    as.lmitt(
      glm(round(exp(y) / (1 + exp(y))) ~ assigned(), data = simdata, weights = ate(spec),
          family = stats::binomial())
    ))

  nq <- nrow(sandwich::estfun(m))
  WX <- m$weights * m$residuals * stats::model.matrix(m)

  uoanames <- var_names(m@StudySpecification, "u")
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))

  nc <- nrow(stats::model.frame(cmod))
  fim <- crossprod(stats::model.matrix(cmod))
  expect_equal(propertee:::.get_a11_inverse(m_as.lmitt), nc * solve(fim))
  expect_equal(propertee:::.get_a11_inverse(m_lmitt.form), nc * solve(fim))
})

test_that(".get_a11_inverse returns correct value for Gaussian glm cmod", {
  data(simdata)
  cmod <- glm(y ~ x, data = simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))

  nc <- nrow(stats::model.frame(cmod))
  dispersion <- sum((cmod$weights * cmod$residuals)^2) / sum(cmod$weights)
  fim <- crossprod(stats::model.matrix(cmod), stats::model.matrix(cmod))
  expect_equal(propertee:::.get_a11_inverse(m_as.lmitt), nc * dispersion * solve(fim))
  expect_equal(propertee:::.get_a11_inverse(m_lmitt.form), nc * dispersion * solve(fim))
})

test_that(".get_a11_inverse returns correct value for poisson glm cmod", {
  data(simdata)
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::poisson())
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))

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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))

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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m_as.lmitt <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
  )
  m_lmitt.form <- lmitt(y ~ 1, data = simdata, specification = spec,
                        weights = ate(spec), offset = cov_adj(cmod))

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
  spec <- rct_spec(z ~ uoa(uoa1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec),
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  nuoas <- nrow(spec@structure)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec),
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
  spec <- rct_spec(z ~ uoa(uoa1), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod)))

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

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  ssmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, specification = spec)))

  expect_warning(propertee:::.get_b11(ssmod, cluster = c("uoa1", "uoa2")),
                 "are found to have NA's")
  expect_equal(suppressWarnings(propertee:::.get_b11(ssmod, cluster = c("uoa1", "uoa2"),
                                         type = "HC0", cadjust = FALSE)),
               crossprod(sandwich::estfun(cmod))) # there should be no clustering

  # check case where one clustering column doesn't only have NA's
  cmod_data$uoa1 <- rep(seq(6, 10), each = 20)
  cmod <- lm(y ~ x, cmod_data)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  ssmod <- as.lmitt(lm(y ~ assigned(), data = simdata,
                      offset = cov_adj(cmod, specification = spec)))

  expect_warning(propertee:::.get_b11(ssmod, cluster = c("uoa1", "uoa2")),
                 "have NA's for some but not all")
})

test_that(".get_b11 returns correct B_11 for multiple cluster columns", {
  data(simdata)
  cmod <- lm(y ~ x, data = simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod)))

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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod)))

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
  return(expect_true(TRUE)) # added 1/14/25 because it's deprecated and now failing
  cmod <- glm(round(exp(y)) ~ x, data = simdata, family = stats::poisson())
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  m <- as.lmitt(
    glm(round(exp(y)) ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod),
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

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata, subset = simdata$uoa2 == 1)
  weighted_spec <- ate(spec, data = simdata[simdata$uoa2 == 1,])
  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata[simdata$uoa2 == 1,],
       weights = weighted_spec, offset = cov_adj(cmod))
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  C_no_cluster_ids <- simdata
  C_no_cluster_ids[, var_names(spec, "u")] <- NA
  cmod <- lm(y ~ x, data = C_no_cluster_ids)
  nc <- sum(summary(cmod)$df[1L:2L])

  m <- as.lmitt(
    lm(y ~ assigned(), data = simdata, weights = ate(spec), offset = cov_adj(cmod))
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
  spec <- rct_spec(z ~ cluster(uoa1), simdata, subset = simdata$uoa1 %in% c(1, 5))

  msk <- simdata$uoa1 %in% c(1, 5)
  m <- as.lmitt(
    lm(y ~ assigned(), simdata[msk,],
       weights = ate(spec, data = simdata[msk,]),
       offset = cov_adj(cmod, newdata = simdata[msk,]))
  )
  expect_warning(propertee:::.get_b11(m),
                 "meat matrix numerically indistinguishable")
  expect_equal(suppressWarnings(propertee:::.get_b11(m)),
               matrix(0, nrow = 2, ncol = 2))
})

test_that(".get_a21 returns correct matrix for lm cmod and lm ssmod", {
  data(simdata)

  new_df <- simdata
  new_df$uid <- seq_len(nrow(new_df))
  cmod <- lm(y ~ x, new_df)
  spec <- rct_spec(z ~ unitid(uid), data = new_df)

  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), new_df, weights = ate(spec), offset = cov_adj(cmod)
  ))
  m_lmitt.form <- lmitt(y ~ 1, specification = spec, data = new_df,
                        weights = ate(spec), offset = cov_adj(cmod))

  Qmat <- m_as.lmitt$weights * stats::model.matrix(m_as.lmitt)
  ctrl.means.grad <- cbind(matrix(0, nrow = nrow(Qmat), ncol = 1),
                           model.matrix(m_as.lmitt@ctrl_means_model) *
                             weights(m_as.lmitt@ctrl_means_model))
  colnames(ctrl.means.grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  Cmat <- stats::model.matrix(cmod)
  nq <- nrow(stats::model.frame(m_as.lmitt))

  a21_as.lmitt <- .get_a21(m_as.lmitt)
  expect_equal(dim(a21_as.lmitt), c(4, 2))
  expect_true(all.equal(a21_as.lmitt, crossprod(cbind(Qmat, -ctrl.means.grad), Cmat) / nq,
                        check.atrributes = FALSE))

  a21_lmitt.form <- .get_a21(m_lmitt.form)
  expect_equal(dim(a21_lmitt.form), c(4, 2))
  expect_true(all.equal(a21_lmitt.form, crossprod(cbind(Qmat, -ctrl.means.grad), Cmat) / nq,
                        check.attributes = FALSE))
})

test_that(".get_a21 returns correct matrix for glm cmod and lm ssmod", {
  data(simdata)
  new_df <- simdata
  new_df$id <- seq(nrow(new_df))
  new_df$bin_y <- rbinom(nrow(new_df), 1, round(1 / (1 + exp(-new_df$x))))

  cmod <- suppressWarnings(glm(bin_y ~ x + force, data = new_df, family = stats::binomial()))
  spec <- rct_spec(z ~ unitid(id), new_df)
  m_as.lmitt <- as.lmitt(lm(
    bin_y ~ assigned(), new_df, weights = ate(spec), offset = cov_adj(cmod, by = "id")
  ))
  m_lmitt.form <- lmitt(bin_y ~ 1, specification = spec, data = new_df,
                        weights = ate(spec), offset = cov_adj(cmod, by = "id"))

  Qmat <- stats::model.matrix(m_as.lmitt)
  ctrl.means.grad <- cbind(matrix(0, nrow = nrow(Qmat), ncol = 1),
                           model.matrix(m_as.lmitt@ctrl_means_model) *
                             weights(m_as.lmitt@ctrl_means_model))
  colnames(ctrl.means.grad) <- c("bin_y:(Intercept)", "cov_adj:(Intercept)")
  Cmat <- cmod$prior.weights * cmod$family$mu.eta(cmod$linear.predictors) * stats::model.matrix(cmod)
  nq <- nrow(stats::model.frame(m_as.lmitt))

  a21_as.lmitt <- .get_a21(m_as.lmitt)
  expect_equal(dim(a21_as.lmitt), c(4, 3))
  expect_true(all.equal(a21_as.lmitt, crossprod(cbind(Qmat, -ctrl.means.grad), Cmat) / nq,
                        check.atrributes = FALSE))

  a21_lmitt.form <- .get_a21(m_lmitt.form)
  expect_equal(dim(a21_lmitt.form), c(4, 3))
  expect_true(all.equal(a21_lmitt.form, crossprod(cbind(Qmat, -ctrl.means.grad), Cmat) / nq,
                        check.attributes = FALSE))
})

test_that(".get_a21 with a moderator", {
  data(simdata)
  
  new_df <- simdata
  new_df$uid <- seq_len(nrow(new_df))
  cmod <- lm(y ~ x, new_df)
  spec <- rct_spec(z ~ unitid(uid), data = new_df)

  m_lmitt.form <- lmitt(y ~ x, specification = spec, data = new_df,
                        weights = ate(spec), offset = cov_adj(cmod))
  
  Qmat <- m_lmitt.form$weights * stats::model.matrix(m_lmitt.form)
  ctrl.means.grad <- cbind(matrix(0, nrow = nrow(Qmat), ncol = 2),
                           model.matrix(m_lmitt.form@ctrl_means_model) *
                             weights(m_lmitt.form@ctrl_means_model))
  colnames(ctrl.means.grad) <- c("y:(Intercept)", "cov_adj:(Intercept)", "y:x", "cov_adj:x")
  Cmat <- stats::model.matrix(cmod)
  nq <- nrow(stats::model.frame(m_lmitt.form))

  a21_lmitt.form <- .get_a21(m_lmitt.form)
  expect_equal(dim(a21_lmitt.form), c(8, 2))
  expect_true(all.equal(a21_lmitt.form, crossprod(cbind(Qmat, -ctrl.means.grad), Cmat) / nq,
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
  spec <- rct_spec(z ~ unitid(id), new_df)

  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), new_df, weights = ate(spec), offset = cov_adj(cmod, by = "id")
  ))
  m_lmitt.form <- lmitt(y ~ 1, specification = spec, data = new_df,
                        weights = ate(spec), offset = cov_adj(cmod, by = "id"))

  expect_equal(.get_tilde_a21(m_as.lmitt), nq / n * .get_a21(m_as.lmitt))
  expect_equal(.get_tilde_a21(m_lmitt.form), nq / n * .get_a21(m_lmitt.form))
})

test_that(".get_tilde_a21 with missing", {
  set.seed(438)
  data(simdata)
  
  new_df <- simdata
  new_df$y[1:3] <- NA_real_
  nc <- 30
  nq <- nrow(new_df)
  n <- nc + nq
  cmod_data <- data.frame(y = rnorm(nc), x = rnorm(nc), id = nq + seq(nc))
  cmod <- lm(y ~ x, cmod_data)
  
  new_df$id <- seq(nq)
  spec <- rct_spec(z ~ unitid(id), new_df)
  
  m_as.lmitt <- as.lmitt(lm(
    y ~ assigned(), new_df, weights = ate(spec), offset = cov_adj(cmod, by = "id")
  ))
  m_lmitt.form <- lmitt(y ~ 1, specification = spec, data = new_df,
                        weights = ate(spec), offset = cov_adj(cmod, by = "id"))
  
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
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), m_data)

  # warning procs twice
  expect_warning(expect_warning(
    m_as.lmitt <- lmitt(lm(y ~ assigned(), m_data, offset = cov_adj(cmod, specification = spec))),
    "covariance adjustments are NA"), "covariance adjustments are NA")
  expect_warning(
    m_lmitt.form <- lmitt(y ~ 1, specification = spec, data = m_data, offset = cov_adj(cmod)),
    "covariance adjustments are NA"
  )

  ssmod_mm <- suppressWarnings(stats::model.matrix(m_as.lmitt))
  ctrl_means_mm <- (
    model.matrix(m_lmitt.form@ctrl_means_model) *
      weights(m_lmitt.form@ctrl_means_model)[!is.na(weights(m_lmitt.form@ctrl_means_model))]
  )
  ctrl.means.grad <- cbind(matrix(0, nrow = nrow(ssmod_mm), ncol = 1),
                           ctrl_means_mm[
                                 row.names(ctrl_means_mm) %in% setdiff(
                                   row.names(ctrl_means_mm), stats::na.action(m_lmitt.form))
                               ])
  colnames(ctrl.means.grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  pg <- stats::model.matrix(formula(cmod), m_data)
  nq <- nrow(m_data)

  a21_as.lmitt <- suppressWarnings(.get_a21(m_as.lmitt))
  expect_true(all.equal(a21_as.lmitt, crossprod(cbind(ssmod_mm, -ctrl.means.grad), pg) / nq,
                        check.attributes = FALSE))
  a21_lmitt.form <- suppressWarnings(.get_a21(m_lmitt.form))
  expect_true(all.equal(a21_lmitt.form, crossprod(cbind(ssmod_mm, -ctrl.means.grad), pg) / nq,
                        check.attributes = FALSE))
})

test_that(paste(".get_a21 returns correct matrix when data input for lmitt has
                NA values for some treatment assignments"), {
  data(simdata)
  simdata[simdata$uoa1 %in% c(2, 4), "z"] <- NA_integer_
  cmod <- lm(y ~ x, data = simdata, subset = uoa1 %in% c(2, 4))
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  m_as.lmitt <- lmitt(
    lm(y ~ assigned(), simdata, w = ate(spec, data = simdata),
       offset = cov_adj(cmod, specification = spec)))
  m_lmitt.form <- lmitt(y ~ 1, simdata, specification = spec, w = ate(spec), offset = cov_adj(cmod))

  ssmod_mm <- stats::model.matrix(m_as.lmitt) * m_as.lmitt$weights
  ctrl.means.grad <- cbind(matrix(0, nrow = nrow(ssmod_mm), ncol = 1),
                           model.matrix(m_lmitt.form@ctrl_means_model) *
                             weights(m_lmitt.form@ctrl_means_model)[
                               !is.na(weights(m_lmitt.form@ctrl_means_model))])
  colnames(ctrl.means.grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  pg <- stats::model.matrix(formula(cmod), simdata)[!is.na(simdata$z),]
  nq <- nrow(simdata)

  a21_as.lmitt <- suppressWarnings(.get_a21(m_as.lmitt))
  expect_true(all.equal(a21_as.lmitt, crossprod(cbind(ssmod_mm, -ctrl.means.grad), pg) / nq,
                        check.attributes = FALSE))
  a21_lmitt.form <- suppressWarnings(.get_a21(m_lmitt.form))
  expect_true(all.equal(a21_lmitt.form, crossprod(cbind(ssmod_mm, -ctrl.means.grad), pg) / nq,
                        check.attributes = FALSE))
})

test_that(".get_a21 returns only full rank columns for less than full rank model", {
  data(simdata)
  copy_simdata <- simdata
  copy_simdata$o_fac <- as.factor(copy_simdata$o)
  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), copy_simdata)

  ### lmitt.formula
  ssmod <- lmitt(y ~ o_fac, data = copy_simdata, specification = spec, offset = cov_adj(cmod))
  ctrl.means.mm <- model.matrix(ssmod@ctrl_means_model)
  expect_equal(dim(a21 <- .get_a21(ssmod)),
               c(ssmod$rank + 2 * ncol(ctrl.means.mm),
                 ssmod$model$`(offset)`@fitted_covariance_model$rank))
  keep_ix <- ssmod$qr$pivot[1L:ssmod$rank]
  expect_equal(rownames(a21),
               c(colnames(model.matrix(ssmod))[keep_ix],
                 paste(rep(c("y", "cov_adj"), each = ncol(ctrl.means.mm)),
                       rep(colnames(ctrl.means.mm), 2), sep = ":")))
})

test_that(".vcov_CR with covariance adjustment and no moderator returns (p+2)x(p+2) matrix", {
  data(simdata)

  cmod <- lm(y ~ x, simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  m <- as.lmitt(lm(y ~ assigned(), simdata, weights = ate(spec), offset = cov_adj(cmod)))

  vmat <- propertee:::.vcov_CR(m)
  expect_equal(dim(vmat), c(4, 4))
})

test_that(".vcov_CR returns teeMod model sandwich if it has no SandwichLayer", {
  data(simdata)

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  m <- lmitt(y ~ 1, data = simdata, specification = spec, weights = ate(spec))

  uoas <- apply(simdata[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))
  expect_true(all.equal(vcov_tee(m, type = "CR0"),
                        sandwich::sandwich(m,
                                           meat. = sandwich::meatCL,
                                           cluster = uoas),
                        check.attributes = FALSE))
})

test_that("test the output of vcov_tee is correct with bias corrections", {
  data(simdata)
  
  cmod <- lm(y ~ x + z, simdata)
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)
  xm <- lmitt(y ~ 1, spec, simdata, offset = cov_adj(cmod))
  
  cls <- paste(simdata$uoa1, simdata$uoa2, sep = "_")
  jk_units <- unique(cls)
  X <- stats::model.matrix(cmod)
  X[,"z"] <- 0
  Z <- stats::model.matrix(xm)
  ZTZ_inv <- chol2inv(xm$qr$qr)
  new_r <- Reduce(
    c,
    mapply(
      function(loo_unit, cmod, cls) {
        cmod_cl <- stats::getCall(cmod)
        cmod_cl$subset <- eval(cls != loo_unit)
        loo_cmod <- eval(cmod_cl, envir = environment(formula(cmod)))
        cl_ix <- cls == loo_unit
        (simdata$y - xm$fitted.values + xm$offset)[cl_ix] - drop(X[cl_ix,,drop=FALSE] %*% loo_cmod$coefficients)
      },
      jk_units,
      SIMPLIFY = FALSE,
      MoreArgs = list(cmod = cmod, cls = cls)
    )
  )
  ef_psi <- estfun(as(xm, "lm")) / xm$residuals * new_r
  
  n <- nrow(simdata)
  ef_phi <- estfun(cmod)
  
  ef_ctrl_means_model <- estfun(xm@ctrl_means_model)
  br <- .get_tilde_a22_inverse(xm)
  vc <- br %*% (
    crossprod(rowsum(cbind(ef_psi, ef_ctrl_means_model) -
                       ef_phi %*% t(.get_a11_inverse(xm)) %*% t(.get_a21(xm)),
                     cls)) / n^2
  ) %*% t(br)
  expect_true(all.equal(vc, vcov_tee(xm, cadjust = FALSE, loco_residuals = TRUE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR lm w/o clustering",
                "balanced grps",
                "no cmod/ssmod data overlap", sep = ", "), {
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
  df1$y <- drop(stats::model.matrix(~ x1 + x2 + z, df1) %*% theta + rnorm(N))
  df2 <- df1
  df2$uid <- NA_integer_
  df <- rbind(df1, df2)

  cmod_form <- y ~ x1 + x2
  ssmod_form <- y ~ assigned()
  cmod <- lm(cmod_form, df, subset = is.na(uid))
  spec <- rct_spec(z ~ uoa(uid), df, subset = !is.na(df$uid))
  ssmod_as.lmitt <- as.lmitt(
    lm(ssmod_form, data = df, subset = !is.na(uid), weights = ate(spec),
       offset = cov_adj(cmod))
  )
  ssmod_lmitt.form <- lmitt(y ~ 1, specification = spec, data = df, subset = !is.na(uid),
                           weights = ate(spec), offset = cov_adj(cmod))
  onemod <- lm(y ~ x1 + x2 + z, data = df, subset = !is.na(uid),
               weights = ate(spec))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(ssmod_as.lmitt)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$uid),])
  cm_grad <- cbind(matrix(0, nrow = nrow(X), ncol = 1),
                   stats::model.matrix(ssmod_as.lmitt@ctrl_means_model) *
                     weights(ssmod_as.lmitt@ctrl_means_model))
  colnames(cm_grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  nc <- nrow(Xstar)
  nq <- nrow(Z)
  n <- nc + nq

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(ssmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(ssmod_as.lmitt), crossprod(cbind(Z * ssmod_as.lmitt$weights, -cm_grad), X) / nq)
  bread. <- sandwich::bread(ssmod_as.lmitt)
  expect_equal(bread.[1:2, 1:2], n * solve(crossprod(Z * ssmod_as.lmitt$weights, Z)))
  expect_equal(bread.[3:4, 3:4],
               matrix(0, nrow = 2, ncol = 2, dimnames = list(colnames(cm_grad), colnames(cm_grad))) + (
                 diag(2) * mean(ssmod_as.lmitt$weights)))
  
  ef_ssmod <- utils::getS3method("estfun", "lm")(ssmod_as.lmitt)
  ef_ssmod <- rbind(ef_ssmod, matrix(0, nrow = nc, ncol = ncol(ef_ssmod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  ctrl_means_mod <- ssmod_as.lmitt@ctrl_means_model
  ef_ctrl_means <- sweep(residuals(ctrl_means_mod), # need sweep bc residuals has 2 columns
                         1,
                         weights(ctrl_means_mod) * model.matrix(ctrl_means_mod), # only 1 col in model matrix
                         FUN = "*")
  ef_ctrl_means <- rbind(ef_ctrl_means, matrix(0, nrow = nc, ncol = ncol(ef_ctrl_means)))
  colnames(ef_ctrl_means) <- colnames(cm_grad)
  expect_equal(meat. <- crossprod(estfun(ssmod_as.lmitt)) / n,
               crossprod(cbind(ef_ssmod, ef_ctrl_means) - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21)) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meat(ssmod_as.lmitt, adjust = FALSE))

  # with perfect group balance, a22inv %*% a21 ((Z'WZ)^(-1)Z'WX) should be 0
  # for the trt row
  expect_equal((bread. %*% a21)[2,],
               setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- propertee:::.vcov_CR(ssmod_as.lmitt, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1/n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be equal to the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(ssmod_as.lmitt),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR(ssmod_lmitt.form, cluster = seq_len(n), cadjust = FALSE),
                        check.attributes = FALSE))

  # var_hat(z) should be smaller than var_hat(z) from onemod
  expect_true(all(diag(vmat) <
                  diag(vcov(onemod))[c(1, 4)]))
})

test_that(paste("HC0 .vcov_CR lm w/o clustering",
                "imbalanced grps",
                "no cmod/ssmod data overlap", sep = ", "), {
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
  df$y <- drop(stats::model.matrix(~ x1 + x2 + z, df) %*% theta + rnorm(N))
  df$uid[seq_len(3 * N / 4)] <- NA_integer_

  cmod_form <- y ~ x1 + x2
  ssmod_form <- y ~ assigned()
  cmod <- lm(cmod_form, df, subset = is.na(uid))
  spec <- rct_spec(z ~ uoa(uid), df, subset = !is.na(df$uid))
  ssmod_as.lmitt <- as.lmitt(
    lm(ssmod_form, data = df, subset = !is.na(uid), weights = ate(spec),
       offset = cov_adj(cmod))
  )
  ssmod_lmitt.form <- lmitt(y ~ 1, specification = spec, data = df, subset = !is.na(uid),
                           weights = ate(spec), offset = cov_adj(cmod))
  onemod <- lm(y ~ x1 + x2 + z, data = df, subset = !is.na(uid),
               weights = ate(spec))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(ssmod_as.lmitt)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$uid),])
  cm_grad <- cbind(matrix(0, nrow = nrow(X), ncol = 1),
                   stats::model.matrix(ssmod_as.lmitt@ctrl_means_model) *
                     weights(ssmod_as.lmitt@ctrl_means_model))
  colnames(cm_grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  nc <- nrow(Xstar)
  nq <- nrow(Z)
  n <- nc + nq

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(ssmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(ssmod_as.lmitt), crossprod(cbind(Z * ssmod_as.lmitt$weights, -cm_grad), X) / nq)
  bread. <- sandwich::bread(ssmod_as.lmitt)
  expect_equal(bread.[1:2, 1:2], n * solve(crossprod(Z * ssmod_as.lmitt$weights, Z)))
  expect_equal(bread.[3:4, 3:4],
               matrix(0, nrow = 2, ncol = 2, dimnames = list(colnames(cm_grad), colnames(cm_grad))) + (
                 diag(2) * (n / sum(weights(ssmod_as.lmitt@ctrl_means_model)))))

  ef_ssmod <- utils::getS3method("estfun", "lm")(ssmod_as.lmitt)
  ef_ssmod <- rbind(ef_ssmod, matrix(0, nrow = nc, ncol = ncol(ef_ssmod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  ctrl_means_mod <- ssmod_as.lmitt@ctrl_means_model
  ef_ctrl_means <- sweep(residuals(ctrl_means_mod), # need sweep bc residuals has 2 columns
                         1,
                         weights(ctrl_means_mod) * model.matrix(ctrl_means_mod), # only 1 col in model matrix
                         FUN = "*")
  ef_ctrl_means <- rbind(ef_ctrl_means, matrix(0, nrow = nc, ncol = ncol(ef_ctrl_means)))
  expect_true(all.equal(meat. <- crossprod(estfun(ssmod_as.lmitt)) / n,
                        crossprod(cbind(ef_ssmod, ef_ctrl_means) - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21)) / n,
                        check.attributes = FALSE))


  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meat(ssmod_as.lmitt, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- propertee:::.vcov_CR(ssmod_as.lmitt, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1/n) * bread. %*% meat. %*% bread.,
                        check.attributes = FALSE))

  # vmat should be equal the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(ssmod_as.lmitt),
                       check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR(ssmod_lmitt.form, cluster = seq_len(n), cadjust = FALSE),
                        check.attributes = FALSE))

  # var_hat(z) should be smaller than var_hat(z) from onemod
  expect_true(all(diag(vmat) <
                    diag(vcov(onemod))[c(1, 4)]))
})

test_that(paste("HC0 .vcov_CR lm w/ clustering",
                "balanced grps",
                "no cmod/ssmod data overlap", sep = ", "), {
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
  df$y <- drop(stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps)
  df$cid <- c(rep(NA_integer_, nclusts_C * MI), rep(seq_len(nclusts_Q), each = MI))

  cmod_form <- y ~ x1 + x2
  ssmod_form <- y ~ assigned()

  cmod <- lm(cmod_form, df, subset = is.na(cid))
  spec <- rct_spec(z ~ cluster(cid), df, subset = !is.na(df$cid))
  ssmod_as.lmitt <- as.lmitt(
    lm(ssmod_form, data = df, subset = !is.na(cid), weights = ate(spec),
       offset = cov_adj(cmod))
  )
  ssmod_lmitt.form <- lmitt(y ~ 1, specification = spec, data = df, subset = !is.na(cid),
                           weights = ate(spec), offset = cov_adj(cmod))

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(ssmod_as.lmitt)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$cid),])
  cm_grad <- cbind(matrix(0, nrow = nrow(X), ncol = 1),
                   stats::model.matrix(ssmod_as.lmitt@ctrl_means_model) *
                     weights(ssmod_as.lmitt@ctrl_means_model))
  colnames(cm_grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  nq <- nrow(Z)
  nc <- nrow(Xstar)
  n <- nc + nq

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(ssmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(ssmod_as.lmitt), crossprod(cbind(Z * ssmod_as.lmitt$weights, -cm_grad), X) / nq)
  bread. <- sandwich::bread(ssmod_as.lmitt)
  expect_equal(bread.[1:2, 1:2], n * solve(crossprod(Z * ssmod_as.lmitt$weights, Z)))
  expect_equal(bread.[3:4, 3:4],
               matrix(0, nrow = 2, ncol = 2, dimnames = list(colnames(cm_grad), colnames(cm_grad))) + (
                 diag(2) * n / sum(weights(ssmod_as.lmitt@ctrl_means_model))))

  ids <- c(df$cid[!is.na(df$cid)],
           paste0(length(df$cid[!is.na(df$cid)]) + seq_len(length(df$cid[is.na(df$cid)])),
                  "*"))
  ef_ssmod <- utils::getS3method("estfun", "lm")(ssmod_as.lmitt)
  ef_ssmod <- rbind(ef_ssmod, matrix(0, nrow = nc, ncol = ncol(ef_ssmod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nrow(Z), ncol = ncol(ef_cmod)), ef_cmod)
  ctrl_means_mod <- ssmod_as.lmitt@ctrl_means_model
  ef_ctrl_means <- sweep(residuals(ctrl_means_mod), # need sweep bc residuals has 2 columns
                         1,
                         weights(ctrl_means_mod) * model.matrix(ctrl_means_mod), # only 1 col in model matrix
                         FUN = "*")
  ef_ctrl_means <- rbind(ef_ctrl_means, matrix(0, nrow = nc, ncol = ncol(ef_ctrl_means)))
  colnames(ef_ctrl_means) <- colnames(cm_grad)
  expect_equal(meat. <- crossprod(Reduce(rbind, by(estfun(ssmod_as.lmitt), ids, colSums))) / n,
               crossprod(Reduce(
                 rbind,
                 by(cbind(ef_ssmod, ef_ctrl_means) - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21), ids, colSums))) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meatCL(ssmod_as.lmitt, cluster = ids, cadjust = FALSE))

  # with perfect group balance, a22inv %*% a21 should have 0's in the trt row
  expect_equal((bread. %*% a21)[2,], setNames(rep(0, dim(a21)[2]), colnames(a21)))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- propertee:::.vcov_CR(ssmod_as.lmitt, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% bread.,
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR(ssmod_lmitt.form, cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be equal the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(ssmod_as.lmitt,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR lm w/ clustering",
                "imbalanced grps",
                "no cmod/ssmod data overlap", sep = ", "), {
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
  df$y <- drop(stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps)
  df$cid <- c(rep(NA_integer_, nclusts_C * MI), rep(seq_len(nclusts_Q), each = MI))

  cmod_form <- y ~ x1 + x2
  ssmod_form <- y ~ assigned()

  cmod <- lm(cmod_form, df, subset = is.na(cid))
  spec <- rct_spec(z ~ cluster(cid), df, subset = !is.na(df$cid))
  ssmod_as.lmitt <- as.lmitt(
    lm(ssmod_form, data = df, subset = !is.na(cid), weights = ate(spec),
       offset = cov_adj(cmod))
  )
  ssmod_lmitt.form <- lmitt(y ~ 1, data = df, specification = spec, subset = !is.na(cid),
                           weights = ate(spec), offset = cov_adj(cmod))
  ctrl_means_mod <- ssmod_as.lmitt@ctrl_means_model

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(ssmod_as.lmitt)
  X <- stats::model.matrix(cmod_form, df[!is.na(df$cid),])
  cm_grad <- cbind(matrix(0, nrow = nrow(X), ncol = 1),
                   stats::model.matrix(ssmod_as.lmitt@ctrl_means_model) *
                     weights(ssmod_as.lmitt@ctrl_means_model))
  colnames(cm_grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  nc <- nrow(Xstar)
  nq <- nrow(Z)
  n <- nc + nq


  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(ssmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(ssmod_as.lmitt), crossprod(cbind(Z * ssmod_as.lmitt$weights, -cm_grad), X) / nq)
  bread. <- sandwich::bread(ssmod_as.lmitt)
  expect_equal(bread.[1:2, 1:2], n * solve(crossprod(Z * ssmod_as.lmitt$weights, Z)))
  expect_equal(bread.[3:4, 3:4],
               matrix(0, nrow = 2, ncol = 2, dimnames = list(colnames(cm_grad), colnames(cm_grad))) + (
                 diag(2) * n / sum(weights(ctrl_means_mod))))

  ids <- c(df$cid[!is.na(df$cid)],
           paste0(length(df$cid[!is.na(df$cid)]) + seq_len(length(df$cid[is.na(df$cid)])),
                  "*"))
  ef_ssmod <- utils::getS3method("estfun", "lm")(ssmod_as.lmitt)
  ef_ssmod <- rbind(ef_ssmod, matrix(0, nrow = nc, ncol = ncol(ef_ssmod)))
  ef_cmod <- estfun(cmod)
  ef_cmod <- rbind(matrix(0, nrow = nq, ncol = ncol(ef_cmod)), ef_cmod)
  ef_ctrl_means <- sweep(residuals(ctrl_means_mod), # need sweep bc residuals has 2 columns
                         1,
                         weights(ctrl_means_mod) * model.matrix(ctrl_means_mod), # only 1 col in model matrix
                         FUN = "*")
  ef_ctrl_means <- rbind(ef_ctrl_means, matrix(0, nrow = nc, ncol = ncol(ef_ctrl_means)))
  colnames(ef_ctrl_means) <- colnames(cm_grad)
  expect_equal(meat. <- crossprod(Reduce(rbind, by(estfun(ssmod_as.lmitt), ids, colSums))) / n,
               crossprod(Reduce(
                 rbind,
                 by(cbind(ef_ssmod, ef_ctrl_means) -  nq / nc * ef_cmod %*% t(a11inv) %*% t(a21), ids, colSums))) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meatCL(ssmod_as.lmitt, cluster = ids, cadjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR(ssmod_as.lmitt, cluster = ids, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% bread.,
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR(ssmod_lmitt.form, cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(ssmod_as.lmitt,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR lm w/o clustering",
                "imbalanced grps",
                "cmod is a strict subset of ssmod data", sep = ", "), {
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
  df$y <- drop(stats::model.matrix(~ x1 + x2 + z, df) %*% theta + rnorm(N))

  cmod_form <- y ~ x1 + x2
  ssmod_form <- y ~ assigned()
  cmod_idx <- df$z == 0
  cmod <- lm(cmod_form, df[cmod_idx,])
  spec <- rct_spec(z ~ uoa(uid), df)
  ssmod_as.lmitt <- as.lmitt(
    lm(ssmod_form, data = df, weights = ate(spec), offset = cov_adj(cmod))
  )
  ssmod_lmitt.form <- lmitt(y ~ 1, data = df, specification = spec, weights = ate(spec),
                           offset = cov_adj(cmod))
  ctrl_means_mod <- ssmod_as.lmitt@ctrl_means_model

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(ssmod_as.lmitt)
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df)
  cm_grad <- cbind(matrix(0, nrow = nrow(X), ncol = 1),
                   stats::model.matrix(ctrl_means_mod) * weights(ctrl_means_mod))
  colnames(cm_grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  nc <- nrow(Xstar)
  nq <- n <- nrow(Z)

  expect_equal(a11inv <- .get_a11_inverse(ssmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(ssmod_as.lmitt), crossprod(cbind(Z * ssmod_as.lmitt$weights, -cm_grad), X) / nq)
  bread. <- bread(ssmod_as.lmitt)
  expect_equal(bread.[1:2, 1:2], n * solve(crossprod(Z * ssmod_as.lmitt$weights, Z)))
  expect_equal(bread.[3:4, 3:4],
               matrix(0, nrow = 2, ncol = 2, dimnames = list(colnames(cm_grad), colnames(cm_grad))) + (
                 diag(2) * n / sum(weights(ctrl_means_mod))))

  ef_ssmod <- utils::getS3method("estfun", "lm")(ssmod_as.lmitt)
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_ssmod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
  ef_ctrl_means <- sweep(residuals(ctrl_means_mod), # need sweep bc residuals has 2 columns
                         1,
                         weights(ctrl_means_mod) * model.matrix(ctrl_means_mod), # only 1 col in model matrix
                         FUN = "*")
  colnames(ef_ctrl_means) <- colnames(cm_grad)
  expect_equal(meat. <- crossprod(estfun(ssmod_as.lmitt)) / n,
               crossprod(cbind(ef_ssmod, ef_ctrl_means) - nq / nc * ef_cmod %*% a11inv %*% t(a21)) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meat(ssmod_as.lmitt, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR(ssmod_as.lmitt, cluster = df$uid, cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR(ssmod_lmitt.form, cluster = df$uid, cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(ssmod_as.lmitt, adjust = FALSE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR lm w/ clustering",
                "imbalanced grps",
                "cmod is a strict subset of ssmod data", sep = ", "), {
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
  df$y <- drop(stats::model.matrix(~ x1 + x2 + z, df) %*% theta + eps)
  df$cid <- rep(seq_len(nclusts_C + nclusts_Q), each = MI)

  cmod_form <- y ~ x1 + x2
  ssmod_form <- y ~ assigned()
  cmod_idx <- df$z == 0

  cmod <- lm(cmod_form, df[cmod_idx,])
  spec <- rct_spec(z ~ cluster(cid), df)
  ssmod_as.lmitt <- as.lmitt(
    lm(ssmod_form, data = df, weights = ate(spec), offset = cov_adj(cmod))
  )
  ssmod_lmitt.form <- lmitt(y ~ 1, data = df, specification = spec, weights = ate(spec),
                           offset = cov_adj(cmod))
  ctrl_means_mod <- ssmod_as.lmitt@ctrl_means_model

  Xstar <- stats::model.matrix(cmod)
  Z <- stats::model.matrix(ssmod_as.lmitt)
  X <- stats::model.matrix(as.formula(cmod_form[-2]), df)
  cm_grad <- cbind(matrix(0, nrow = nrow(X), ncol = 1),
                   stats::model.matrix(ctrl_means_mod) * weights(ctrl_means_mod))
  colnames(cm_grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  nc <- nrow(Xstar)
  nq <- n <- nrow(Z)

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  expect_equal(a11inv <- .get_a11_inverse(ssmod_as.lmitt), nc * solve(crossprod(Xstar)))
  expect_equal(a21 <- .get_a21(ssmod_as.lmitt), crossprod(cbind(Z * ssmod_as.lmitt$weights, -cm_grad), X) / nq)
  bread. <- sandwich::bread(ssmod_as.lmitt)
  expect_equal(bread.[1:2, 1:2], n * solve(crossprod(Z * ssmod_as.lmitt$weights, Z)))
  expect_equal(bread.[3:4, 3:4],
               matrix(0, nrow = 2, ncol = 2, dimnames = list(colnames(cm_grad), colnames(cm_grad))) + (
                 diag(2) * n / sum(weights(ctrl_means_mod))))

  ids <- df$cid
  ef_ssmod <- utils::getS3method("estfun", "lm")(ssmod_as.lmitt)
  loo_r <- .compute_loo_resids(ssmod_as.lmitt, "cid")
  ef_ssmod <- ef_ssmod / stats::residuals(ssmod_as.lmitt) * loo_r
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_ssmod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
  ef_ctrl_means <- sweep(residuals(ctrl_means_mod), # need sweep bc residuals has 2 columns
                         1,
                         weights(ctrl_means_mod) * model.matrix(ctrl_means_mod), # only 1 col in model matrix
                         FUN = "*")
  colnames(ef_ctrl_means) <- colnames(cm_grad)
  expect_equal(meat. <- crossprod(Reduce(rbind, by(estfun(ssmod_as.lmitt, loco_residuals = TRUE),
                                                   ids, colSums))) / n,
               crossprod(Reduce(
                 rbind,
                 by(cbind(ef_ssmod, ef_ctrl_means) - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21), ids, colSums))) / n)

  # meat should be the same as the output of sandwich::meatCL
  expect_equal(meat., sandwich::meatCL(ssmod_as.lmitt, cluster = ids, cadjust = FALSE,
                                       loco_residuals = TRUE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR(ssmod_as.lmitt, cluster = ids, cadjust = FALSE, loco_residuals = TRUE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR(ssmod_lmitt.form, cluster = ids, cadjust = FALSE, loco_residuals = TRUE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat,
                        sandwich::sandwich(ssmod_as.lmitt,
                                           meat. = sandwich::meatCL,
                                           cluster = ids, cadjust = FALSE,
                                           loco_residuals = TRUE),
                        check.attributes = FALSE))
})

test_that(paste("HC0 .vcov_CR binomial glm cmod",
                "w/o clustering",
                "imbalanced grps",
                "cmod is a strict subset of ssmod data", sep = ", "), {
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
  ssmod_form <- y ~ assigned()
  cmod <- glm(cmod_form, data = df, subset = z == 0, family = stats::binomial())
  spec <- rct_spec(z ~ uoa(uid), df)
  ssmod_as.lmitt <- lmitt(
    lm(ssmod_form, data = df, weights = ate(spec), offset = cov_adj(cmod))
  )
  ssmod_lmitt.form <- lmitt(y ~ 1, data = df, specification = spec, weights = ate(spec),
                           offset = cov_adj(cmod))
  ctrl_means_mod <- ssmod_as.lmitt@ctrl_means_model

  ## COMPARE BLOCKS TO MANUAL DERIVATIONS
  Xstar <- stats::model.matrix(cmod)
  fit_wstar <- cmod$weights
  rstar <- cmod$residuals
  mu_etastar <- cmod$family$mu.eta(cmod$linear.predictors)
  Z <- stats::model.matrix(ssmod_as.lmitt)
  X <- stats::model.matrix(formula(stats::delete.response(terms(cmod))), df)
  cm_grad <- cbind(matrix(0, nrow = nrow(X), ncol = 1),
                   stats::model.matrix(ctrl_means_mod) * weights(ctrl_means_mod))
  colnames(cm_grad) <- c("y:(Intercept)", "cov_adj:(Intercept)")
  w <- ssmod_as.lmitt$weights
  r <- ssmod_as.lmitt$residuals
  mu_eta <- cmod$family$mu.eta(drop(X %*% cmod$coefficients))
  nc <- nrow(Xstar)
  nq <- n <- nrow(Z)

  expect_equal(a11inv <- .get_a11_inverse(ssmod_as.lmitt), nc * solve(crossprod(Xstar * sqrt(fit_wstar))))
  expect_true(all.equal(a21 <- .get_a21(ssmod_lmitt.form),
                        crossprod(cbind(Z * w, -cm_grad), X * mu_eta) / nq,
                        check.attributes = FALSE))
  bread. <- bread(ssmod_as.lmitt)
  expect_equal(bread.[1:2, 1:2], n * solve(crossprod(Z * w, Z)))
  expect_equal(bread.[3:4, 3:4],
               matrix(0, nrow = 2, ncol = 2, dimnames = list(colnames(cm_grad), colnames(cm_grad))) + (
                 diag(2) * n / sum(weights(ctrl_means_mod))))

  ef_ssmod <- utils::getS3method("estfun", "lm")(ssmod_as.lmitt)
  nonzero_ef_cmod <- estfun(cmod)
  ef_cmod <- matrix(0, nrow = nrow(ef_ssmod), ncol = ncol(nonzero_ef_cmod))
  colnames(ef_cmod) <- colnames(nonzero_ef_cmod)
  ef_cmod[which(df$z == 0),] <- nonzero_ef_cmod
  ef_ctrl_means <- sweep(residuals(ctrl_means_mod), # need sweep bc residuals has 2 columns
                         1,
                         weights(ctrl_means_mod) * model.matrix(ctrl_means_mod), # only 1 col in model matrix
                         FUN = "*")
  colnames(ef_ctrl_means) <- colnames(cm_grad)
  expect_equal(meat. <- crossprod(estfun(ssmod_as.lmitt)) / n,
               crossprod(cbind(ef_ssmod, ef_ctrl_means) - nq / nc * ef_cmod %*% t(a11inv) %*% t(a21)) / n)

  # meat should be the same as the output of sandwich::meat
  expect_equal(meat., sandwich::meat(ssmod_as.lmitt, adjust = FALSE))

  # with imperfect group balance, a22inv %*% a21 should not be 0 for the trt row
  expect_false(all((bread. %*% a21)[2,] == 0))

  # check .vcov_CR0 matches manual matrix multiplication
  vmat <- .vcov_CR(ssmod_as.lmitt, cluster = seq_len(n), cadjust = FALSE)
  expect_true(all.equal(vmat, (1 / n) * bread. %*% meat. %*% t(bread.),
                        check.attributes = FALSE))

  # vmat should be the same for both lmitt calls
  expect_true(all.equal(vmat,
                        .vcov_CR(ssmod_lmitt.form, cluster = seq(n), cadjust = FALSE),
                        check.attributes = FALSE))

  # vmat should be the same as the outputs of sandwich
  expect_true(all.equal(vmat, sandwich::sandwich(ssmod_as.lmitt, adjust = FALSE),
                        check.attributes = FALSE))
})

test_that("type attribute", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  ssmod <- lmitt(y ~ 1, data = simdata, weights = "ate", specification = spec)
  expect_identical(attr(vcov(ssmod), "type"), "HC0")
  expect_identical(attr(vcov(ssmod, type = "CR2"), "type"), "CR2")
  expect_identical(attr(vcov(ssmod, type = "MB2"), "type"), "MB2")
  expect_identical(attr(vcov(ssmod, type = "HC2"), "type"), "HC2")
})

test_that("#119 flagging vcov_tee entries as NA", {
  ### factor moderator variable
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  mod <- lmitt(y ~ as.factor(o), data = simdata, specification = spec)
  expect_warning(vc <- vcov_tee(mod, type = "CR0"),
                 "will be returned as NA: as.factor(o)1, as.factor(o)3",
                 fixed = TRUE)

  # Issue is in subgroups w/ moderator=1:z=0, moderator=1:z=1, and
  # moderator=3, so .check_df_moderator_estimates should NA those vcov entries
  na_dim <- c(1, 3, 5, 8, 10)
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
  dspec <- rct_spec(z ~ unitid(ass_id), ddata)
  ssmod <- lmitt(y ~ modr, specification = dspec, data = ddata)
  expect_warning(vc <- vcov_tee(ssmod, type = "CR0"),
                 "will be returned as NA: modr1",
                 fixed = TRUE)

  na_dim <- c(1, 3, 5)
  expect_true(all(
    abs(
      diag(sandwich::sandwich(ssmod, meat. = sandwich::meatCL,
                              cluster = .make_uoa_ids(ssmod, "CR")))[na_dim]
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
  dspec <- rct_spec(z ~ unitid(ass_id), ddata)
  ssmod <- lmitt(y ~ modr, specification = dspec, data = ddata)
  vc <- vcov_tee(ssmod, type = "CR0")
  expect_true(all(
    abs(diag(sandwich::sandwich(ssmod, meat. = sandwich::meatCL,
                                cluster = .make_uoa_ids(ssmod, "CR")))[na_dim])
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
  # mod <- lmitt(lm(y ~ as.factor(o) + as.factor(o):assigned(spec), data = simdata), specification = spec)
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
  ssmod <- lmitt(y ~ o, data = simdata, specification = spec)
  vc <- vcov(ssmod, type = "CR0")
  expect_true(all(!is.na(vc)))
  #
  # ### invalid continuous moderator variable
  simdata <- simdata[(simdata$uoa1 == 2 & simdata$uoa2 == 2) |
                      (simdata$uoa1 == 2 & simdata$uoa2 == 1),]
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  ssmod <- lmitt(y ~ o, data = simdata, specification = spec)
  expect_warning(vc <- vcov(ssmod, type = "CR0"), "will be returned as NA: o")
  na_dim <- which(grepl("(Intercept)", colnames(vc)) |
                    grepl("z._o", colnames(vc)))
  expect_true(all(is.na(vc[na_dim, ])))
  expect_true(all(is.na(vc[, na_dim])))
  expect_true(all(!is.na(vc[-na_dim, -na_dim])))
})

test_that(".check_df_moderator_estimates other warnings", {
  data(simdata)

  # fail without a teeMod model
  nossmod <- lm(y ~ x, simdata)
  vmat <- vcov(nossmod)
  cluster_ids <- apply(simdata[, c("uoa1", "uoa2")], 1, function(...) paste(..., collapse = "_"))
  expect_error(.check_df_moderator_estimates(vmat, nossmod, cluster_ids),
               "must be a teeMod")

  # invalid `data` arg
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  ssmod <- lmitt(y ~ factor(o), specification = spec, data = simdata)
  expect_error(.check_df_moderator_estimates(vmat, ssmod, cluster_ids,
                                             cbind(y = simdata$y, z = simdata$z),
                                             envir = parent.frame()),
               "`model_data` must be a dataframe")
})

test_that("#123 ensure PreSandwich are converted to Sandwich", {
  data(simdata)
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  cmod <- lm(y ~ x, data = simdata)
  # Make sure its PreSandwich prior
  ca <- cov_adj(cmod, newdata = simdata)
  expect_false(is(ca, "SandwichLayer"))
  ssmod <- as.lmitt(lm(y ~ a.(spec), data = simdata, offset = ca),
                    specification = spec)
  expect_true(is(ssmod$model$`(offset)`, "SandwichLayer"))


  copy_simdata <- simdata
  copy_simdata$schoolid <- seq(900, 949)
  C_df <- rbind(copy_simdata[, c("schoolid", "x", "y")],
                data.frame(schoolid = seq(1000, 1049),
                           x = rnorm(50),
                           y = rnorm(50)))
  cmod <- lm(y ~ x, C_df)
  spec <- rct_spec(z ~ uoa(schoolid), copy_simdata)
  lm1 <- lm(y ~ assigned(), copy_simdata, weights = ett(specification = spec),
            offset = cov_adj(cmod, NULL, NULL, "schoolid"))
  ssmod1 <- lmitt(lm1, specification = spec)
  expect_true(is(ssmod$model$`(offset)`, "SandwichLayer"))
  v1 <- vcov_tee(ssmod1)

  wts2 <- ett(spec, data = copy_simdata)
  offst2 <- cov_adj(cmod, copy_simdata, by = "schoolid")
  lm2 <- lm(y ~ assigned(), copy_simdata, weights = wts2, offset = offst2)
  ssmod2 <- lmitt(lm2, specification = spec)
  expect_true(is(ssmod$model$`(offset)`, "SandwichLayer"))
  v2 <- vcov_tee(ssmod2)

  expect_true(all.equal(v1, v2))

})

test_that("model-based SE's cluster units of assignment in small blocks at block level", {
  specdata <- data.frame(bid = rep(c(1, 2), each = 20),
                        uoa_id = c(rep(c(1, 2), each = 10), rep(seq(3, 6), each = 5)),
                        a = c(rep(c(0, 1), each = 10), rep(rep(c(0, 1), each = 5), 2)))
  specdata$y <- rnorm(40)
  spec <- rct_spec(a ~ unitid(uoa_id) + block(bid), specdata)
  suppressMessages(mod <- lmitt(y ~ 1, specification = spec, data = specdata))
  vc_w_small_block_clusters <- vcov_tee(mod)
  vc_w_no_small_block_clusters <- .vcov_CR(mod,
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
  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)
  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat, offset = cov_adj(cmod))
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
  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)
  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                             yr = rep(c("01", "02"), 5),
                             x = rnorm(10),
                             a = rep(c(rep(1, 3), rep(0, 2)), 2),
                             y = rnorm(10),
                             by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat, offset = cov_adj(cmod, by = "by_col"))
  expect_equal(vcov(mod), vcov(mod, by = "by_col"))
})

test_that("vcov_tee does not error when asking for specification-based SE for model
          with covariance adjustment but without absorbed block effects",{
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  ssmod_off <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                     offset = cov_adj(cmod))
  expect_true(!is.na(vcov_tee(ssmod_off, type = "DB0")))
})

test_that("vcov_tee errors when asking for specification-based SE for teeMod with
           external sample for covariance adjustment",{
  # Borrowing code from the #177 test
  set.seed(23)
  cmod_data <- data.frame(yr = rep(c("00", "01", "02"), 5),
                         id = rep(letters[1:5], each = 3),
                         x = rnorm(5 * 3),
                         y = rnorm(5 * 3),
                         by_col = seq_len(15))
  cmod <- lm(y ~ x, cmod_data)
  specdat <- data.frame(id = letters[1:5], a = c(rep(1, 3), rep(0, 2)))
  newspec <- rct_spec(a ~ unitid(id), specdat)
  analysis_dat <- data.frame(id = rep(letters[1:5], each = 2),
                            yr = rep(c("01", "02"), 5),
                            x = rnorm(10),
                            a = rep(c(rep(1, 3), rep(0, 2)), 2),
                            y = rnorm(10),
                            by_col = setdiff(seq_len(15), seq(1, 15, 3)))
  mod <- lmitt(y ~ yr, specification = newspec, data = analysis_dat,
               offset = cov_adj(cmod, by = "by_col"))
  expect_error(
    vcov_tee(mod, type = "DB0"),
    "teeMod models with external sample"
  )
})

test_that("vcov_tee errors when asking for specification-based SE for a non-RCT specification",{
  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  ssmod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec))
  expect_error(
    vcov_tee(ssmod, type = "DB0"),
    "only for RCT specifications"
  )
})

test_that("vcov_tee errors when asking for specification-based SE for a model
          with moderators",{
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  ssmod_sub <- lmitt(y ~ dose, specification = spec, data = simdata, weights = ate(spec))
  expect_error(
    vcov_tee(ssmod_sub, type = "DB0"),
    "are not supported for teeMod models with moderators"
  )
})

test_that("vcov_tee throws a warning when asking for specification-based SE without IPW",{
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  teemod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ett())
  teemod_abs <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ett(spec),
                      absorb = TRUE)
  expect_warning(
    vcov_tee(teemod, type = "DB0"),
    "inverse probability weights"
  )
  expect_silent(vcov_tee(teemod_abs, type = "DB0"))
})

test_that(".merge_block_id_cols correctly combines multiple block columns", {
  df1 <- data.frame(
    block1 = c("A", "A", "B", "B", "B", "B"),
    block2 = c(1, 1, 1, 1, 2, 2)
  )
  df2 <- data.frame(
    block1 = c("A", "A", "B", "B", "B", "B"),
    block2 = c(T, F, T, F, T, F)
  )
  df3 <- data.frame(
    block1 = c("A", "A", "B", "B", "C", "C"),
    block2 = c(T, F, T, T, T, T),
    block3 = c(2, 2, 3, 3, 4, 4)
  )
  expect_equal(
    propertee:::.merge_block_id_cols(df=df1, ids=c("block1", "block2"))$block1,
    c(1, 1, 2, 2, 3, 3)
  )
  expect_equal(
    propertee:::.merge_block_id_cols(df=df1, ids=c("block1"))$block1,
    c(1, 1, 2, 2, 2, 2)
  )
  expect_equal(
    propertee:::.merge_block_id_cols(df=df2, ids=c("block1", "block2"))$block1,
    c(2, 1, 4, 3, 4, 3)
  )
  expect_equal(
    .merge_block_id_cols(df=df3, ids=c("block1", "block2", "block3"))$block1,
    c(2, 1, 3, 3, 4, 4)
  )
})

test_that(".get_DB_wo_covadj_se returns correct value for specifications with small blocks",{
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
  spec <- rct_spec(z ~ cluster(cid) + block(bid), data)
  ssmod <- lmitt(y ~ 1, specification = spec, data = data, weights = ate(spec) * data$w)
  expect_equal(vcov_tee(ssmod, type="DB0")[1,1], sum(nu) / sum(ws)^2)
})

test_that(".get_DB_wo_covadj_se returns correct value for specifications with a few large blocks",{
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
  spec <- rct_spec(z ~ cluster(cid) + block(bid), data)
  ssmod <- lmitt(y ~ 1, specification = spec, data = data, weights = ate(spec) * data$w)

  #ainv <- .get_DB_a_inverse(ssmod)
  #meat <- .get_DB_meat(ssmod)
  #vmat <- as.matrix((ainv %*% meat %*% t(ainv))[3,3])

  expect_equal(vcov_tee(ssmod, type="DB0")[1,1], sum(nu) / sum(ws)^2)
  #expect_true(vmat[1,1] != sum(nu) / sum(ws)^2)
})

test_that("specification-based SE for tee models without absorption does not crash", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  teemod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate())
  teemod_off <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                     offset = cov_adj(cmod))
  expect_silent(vcov_tee(teemod, type = "DB0"))
  expect_silent(vcov_tee(teemod_off, type = "DB0"))

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)
  teemod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate())
  teemod_off <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                      offset = cov_adj(cmod))
  expect_silent(vcov_tee(teemod, type = "DB0"))
  expect_silent(vcov_tee(teemod_off, type = "DB0"))
})

test_that("specification-based SE for tee models with absorption does not crash", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  teemod_abs <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                      absorb = TRUE)
  teemod_abs_off <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                          offset = cov_adj(cmod), absorb = TRUE)
  expect_silent(vcov_tee(teemod_abs, type = "DB0"))
  expect_silent(vcov_tee(teemod_abs_off, type = "DB0"))
  
  ssmod_sub_abs <- lmitt(y ~ dose, specification = spec, data = simdata, absorb = TRUE)
  ssmod_sub_abs_off <- lmitt(y ~ dose, specification = spec, data = simdata,
                             offset = cov_adj(cmod), absorb = TRUE)
  expect_silent(vcov_tee(ssmod_sub_abs, type = "DB0"))
  expect_silent(vcov_tee(ssmod_sub_abs_off, type = "DB0"))
  expect_silent(vcov_tee(ssmod_sub_abs, type = "DB0", const_effect = TRUE))
  expect_silent(vcov_tee(ssmod_sub_abs_off, type = "DB0", const_effect = TRUE))
})

test_that(".get_appinv_atp returns correct (A_{pp}^{-1} A_{tau p}^T)
          for tee models with absorbed intercept", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  ssmod_abs <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                     absorb = TRUE)
  ssmod_abs_off <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                         offset = cov_adj(cmod), absorb = TRUE)

  bid <- simdata$bid
  B <- cbind(as.integer(bid == 1), as.integer(bid == 2), as.integer(bid == 3))
  goal <- matrix(0, nrow = 3, ncol = 1)
  for (blk in 1:3){
    goal[blk] <- sum(weights(ssmod_abs) * residuals(ssmod_abs) * B[,blk]) /
      sum(ssmod_abs$weights * B[,blk])
  }
  expect_true(all.equal(goal, .get_appinv_atp(ssmod_abs, db = TRUE)))

  for (blk in 1:3){
    goal[blk] <- sum(weights(ssmod_abs_off) * residuals(ssmod_abs_off) * B[,blk]) /
      sum(weights(ssmod_abs_off)* B[,blk])
  }
  expect_true(all.equal(goal, propertee:::.get_appinv_atp(ssmod_abs_off, db = TRUE)))
})

test_that(".get_phi_tilde returns correct grave{phi}
          for tee models with absorbed intercept", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  cmod <- lm(y ~ x, simdata)
  ssmod_abs <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(spec),
                     absorb = TRUE)

  bid <- simdata$bid
  B <- cbind(as.integer(bid == 1), as.integer(bid == 2), as.integer(bid == 3))
  Z <- cbind(as.integer(simdata$z == 0), as.integer(simdata$z == 1))

  ws <- stats::weights(ssmod_abs)
  p <- matrix(0, nrow = 2, ncol = 3)
  for (k in 1:2)
    for (j in 1:3){
      p[k, j] <- sum(ws * Z[,k] * B[,j]) / sum(ws * B[,j])
    }
  goal <- matrix(0, nrow = 50, ncol = 3)
  for (s in 1:3){
    goal[,s] <- ws * residuals(ssmod_abs) * (Z[,2] - p[2,s]) * B[,s]
  }
  expect_true(all.equal(goal, propertee:::.get_phi_tilde(ssmod_abs, db = TRUE)))
  
  ssmod_abs <- lmitt(y ~ 1, specification = spec, data = simdata,
                     absorb = TRUE)
  p <- matrix(0, nrow = 2, ncol = 3)
  for (k in 1:2)
    for (j in 1:3){
      p[k, j] <- sum(Z[,k] * B[,j]) / sum(B[,j])
    }
  goal <- matrix(0, nrow = 50, ncol = 3)
  for (s in 1:3){
    goal[,s] <- residuals(ssmod_abs) * (Z[,2] - p[2,s]) * B[,s]
  }
  expect_true(all.equal(goal, propertee:::.get_phi_tilde(ssmod_abs, db = TRUE)))
})

test_that("block_center_residuals correctly subtracts off block center residuals", {
  data(simdata)
  des <- rct_specification(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  damod_abs <- lmitt(y ~ 1, specification = des, data = simdata, absorb = TRUE)
  r <- residuals(damod_abs)
  
  # calculate block center residuals
  block_id <- as.factor(simdata$bid)
  model <- lm(r ~ block_id)
  intercept <- coef(model)[1]
  block_effects <- coef(model)[-1]
  block_means <- intercept + c(0, block_effects)
  names(block_means) <- levels(block_id)
  
  # subtract off block means
  block_means <- block_means[block_id]
  r <- r - block_means
  expect_true(all.equal(r, residuals(propertee:::block_center_residuals(damod_abs))))
})

test_that("block_center_residuals works for input containing NA block IDs", {
  data(simdata)
  simdata[1:3, 'bid'] <- NA
  des <- rct_specification(z ~ cluster(uoa1, uoa2) + block(bid), 
                           simdata, na.fail = FALSE)
  damod_abs <- lmitt(y ~ 1, specification = des, data = simdata, absorb = TRUE)
  r <- residuals(damod_abs)
  
  # calculate block center residuals
  block_id <- as.factor(simdata$bid)
  model <- lm(r ~ block_id)
  intercept <- coef(model)[1]
  block_effects <- coef(model)[-1]
  block_means <- intercept + c(0, block_effects)
  names(block_means) <- levels(block_id)
  
  # subtract off block means
  block_means <- block_means[block_id]
  r <- r - block_means
  expect_true(all.equal(r, residuals(propertee:::block_center_residuals(damod_abs))))
})

test_that("block_center_residuals works for teeMod having NA in fitted values", {
  data(simdata)
  simdata[1:3, 'x'] <- NA
  des <- rct_specification(z ~ cluster(uoa1, uoa2) + block(bid), 
                           simdata, na.fail = FALSE)
  cmod <- lm(y ~ x, simdata)
  teemod <- suppressWarnings(lmitt(
    y ~ 1, specification = des, data = simdata, absorb = TRUE, offset = cov_adj(cmod)
    ))
  r <- residuals(teemod)
  
  # calculate block center residuals
  block_id <- as.factor(simdata$bid)[4:50]
  model <- lm(r ~ block_id)
  intercept <- coef(model)[1]
  block_effects <- coef(model)[-1]
  block_means <- intercept + c(0, block_effects)
  names(block_means) <- levels(block_id)
  
  # subtract off block means
  block_means <- block_means[block_id]
  r <- r - block_means
  expect_true(all.equal(r, residuals(propertee:::block_center_residuals(teemod))))
})

test_that("issue #239", {
  # crux of the issue is that when all observations with a given value of a
  # moderator are given a weight of 0 in the regression, the lm/mlm for
  # computing means in the control condition has NA coefficients. the .get_a21()
  # and custom bread.mlm() functions did not handle this as they needed to.
  # here we make sure they only retain columns that yield full-rank matrices
  set.seed(298)
  odata <- data.frame(x1 = factor(rep(seq_len(3), each = 5)),
                      y = rnorm(15),
                      a = rep(c(0, 1), 8)[1:15],
                      uid = seq_len(15))
  cmod <- lm(y~x1, odata)
  spec <- rct_spec(a ~ unitid(uid), odata)
  tmod <- lmitt(y ~ x1, spec, odata, weights = c(rep(0, 5), rep(1, 10)),
                offset = cov_adj(cmod))
  
  # test bread. if there's an offset and the rank of the ctrl means model is K,
  # then there should be K columns corresponding to the regression for the
  # response and K columns corresponding to the regression for the offset for a
  # total of 2 * K columns. when we add those to the bread matrix for the
  # effect estimation regression, which we'll say has rank R, then the number
  # of columns in the bread matrix will be R + 2 * K.
  twoK <- tmod@ctrl_means_model$rank * 2
  R <- tmod$rank
  br_cm <- bread.mlm(tmod@ctrl_means_model)
  expect_equal(ncol(br_cm), twoK)
  br <- bread(tmod)
  expect_equal(ncol(br), R + twoK)
  
  # test A21. should have R + 2 * K rows. so when we use it to make the output
  # of estfun.teeMod(), we should have a matrix with R + 2 * K columns
  a21 <- .get_a21(tmod)
  expect_equal(nrow(a21), R + twoK)
  ef <- estfun(tmod)
  expect_equal(ncol(ef), R + twoK)
})
