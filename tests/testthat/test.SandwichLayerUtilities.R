set.seed(20)
cont_x <- rnorm(100)
cat_x <- rbinom(100, 2, 0.2)
y <- cont_x - cat_x + rnorm(100)
C_individ <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                        "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                        "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                        "cont_x" = cont_x,
                        "cat_x" = cat_x,
                        "y" = y)
C_cluster_ids <- data.frame(Reduce(
  rbind,
  lapply(
    split(C_individ, list(C_individ$uoa1, C_individ$uoa2)),
    function(x) unique(x[, c("uoa1", "uoa2")]))),
  row.names = NULL)
C_cluster <- cbind(
  C_cluster_ids,
  data.frame(Reduce(
    rbind,
    by(C_individ,
       list(C_individ$uoa1, C_individ$uoa2),
       function(x) {colMeans(x[, c("cont_x", "y")])}),
  ),
  row.names = NULL)) # cat variable doesn't make sense to include for cluster-level cmod
Q <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                "cont_x" = rnorm(100),
                "cat_x" = rbinom(100, 2, 0.2))
des <- rct_design(t ~ uoa(uoa1, uoa2), data = Q)
uoanames <- var_names(des, "u")
zname <- var_names(des, "t")

test_correct_b12 <- function(m,
                             uoanames,
                             zname) {
  # est eqns for lm are wts * resids * dmatrix
  cmod <- m$model$`(offset)`@fitted_covariance_model
  cmod_eqns <- Reduce(
    rbind,
    by(cmod$residuals * model.matrix(cmod),
       lapply(uoanames, function(col) m$model$`(offset)`@keys[,col]),
       colSums))

  Q <- eval(m$call$data)
  Q_ <- .merge_preserve_order(
    Q,
    merge(unique(m$model$`(offset)`@keys), m@Design@structure),
    by = uoanames,
    all.x = TRUE,
    sort = FALSE)

  msk <- !is.na(Q_[, paste0(zname, ".y")])
  m_eqns <- Reduce(
    rbind,
    by((m$weights * m$residuals * model.matrix(m))[msk, ],
       lapply(uoanames, function(col) Q[msk, col]),
       colSums))

  expect_equal(get_overlap_vcov_matrix(m),
               t(cmod_eqns) %*% m_eqns)
}

test_that("get_overlap_vcov_matrix used with DA model without SandwichLayer offset", {
  cmod <- lm(y ~ cont_x + as.factor(cat_x), data = C_individ)
  offset <- stats::predict(cmod, Q)
  m <- as.DirectAdjusted(
    lm(y ~ t, data = Q, weights = ate(des), offset = offset)
  )

  expect_error(get_overlap_vcov_matrix(m),
               "must have an offset of class")
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for indivividual-level",
                "cov model data with full overlap"), {
  cmod <- lm(y ~ cont_x + as.factor(cat_x), data = C_individ)
  pred_gradient <- model.matrix(~ cont_x + as.factor(cat_x), Q)
  sl <- as.SandwichLayer(
   new("PreSandwichLayer",
       stats::predict(cmod, Q),
       fitted_covariance_model = cmod,
       prediction_gradient = pred_gradient),
   des)
  m <- as.DirectAdjusted(
   lm(y ~ t, data = Q, weights = ate(des), offset = sl)
  )

  test_correct_b12(m, uoanames, zname)
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for cluster-level",
                "cov model data with full overlap"), {
  cmod <- lm(y ~ cont_x, data = C_cluster)
  keys <- C_cluster[, c("uoa1", "uoa2")]
  pred_gradient <- model.matrix(~ cont_x, Q)
  slayer <- as.SandwichLayer(
   new("PreSandwichLayer",
       stats::predict(cmod, Q),
       fitted_covariance_model = cmod,
       prediction_gradient = pred_gradient),
   des)
  m <- as.DirectAdjusted(
   lm(y ~ t, data = Q, weights = ate(des), offset = slayer)
  )

  test_correct_b12(m, uoanames, zname)
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for indivividual-level",
                "cov model data with partial overlap (cov data missing experiment",
                "rows)"), {
  on.exit(C_individ <- data.frame("uoa1" = c(rep(1, 50), rep(2, 50)),
                                  "uoa2" = rep(c(rep(1, 25), rep(2, 25)), 2),
                                  "t" = c(rep(1, 25), rep(0, 50), rep(1, 25)),
                                  "cont_x" = cont_x,
                                  "cat_x" = cat_x,
                                  "y" = y))
  C_individ[C_individ$uoa1 == 1 & C_individ$uoa2 == 1, c("uoa1", "uoa2")] <- NA
  cmod <- lm(y ~ cont_x + as.factor(cat_x), data = C_individ)
  pred_gradient <- model.matrix(~ cont_x + as.factor(cat_x), Q)
  slayer <- as.SandwichLayer(
   new("PreSandwichLayer",
       stats::predict(cmod, Q),
       fitted_covariance_model = cmod,
       prediction_gradient = pred_gradient),
   des)
  m <- as.DirectAdjusted(
   lm(y ~ t, data = Q, weights = ate(des), offset = slayer)
  )

  test_correct_b12(m, uoanames, zname)
})

test_that(paste("get_overlap_vcov_matrix returns expected B_12 for indivividual-level",
                "cov model data with partial overlap (experimental data is subset of",
                "cov model data)"), {
  on.exit(des <- rct_design(t ~ uoa(uoa1, uoa2), data = Q), add = TRUE)
  
  cmod <- lm(y ~ cont_x + as.factor(cat_x), data = C_individ)
  des <- rct_design(t ~ uoa(uoa1, uoa2), data = Q, subset = Q$uoa1 == 1)
  pred_gradient <- model.matrix(~ cont_x + as.factor(cat_x), Q[Q$uoa1 == 1,])
  slayer <- as.SandwichLayer(
   new("PreSandwichLayer",
       stats::predict(cmod, Q[Q$uoa1 == 1,]),
       fitted_covariance_model = cmod,
       prediction_gradient = pred_gradient),
   des)
  weighted_design <- ate(des, data = Q[Q$uoa1 == 1,])
  m <- as.DirectAdjusted(
   lm(y ~ t, data = Q[Q$uoa1 == 1,], weights = weighted_design, offset = slayer)
  )
  
  test_correct_b12(m, uoanames, zname)
})