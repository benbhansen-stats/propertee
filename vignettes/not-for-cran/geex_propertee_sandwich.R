library(propertee)

## NOT CLUSTERED DATA
set.seed(980)
nc <- 2
mi <- 20
df <- data.frame(
  "uid" = seq_len(mi * nc),
  "z" = rep(c(0, 1), each = mi),
  "x1" = rnorm(sd = 2, n = mi * nc),
  "x2" = rep(c(rep(0, mi * 0.5), rep(1, mi * 0.5)), nc / 2)
)
eps <- rnorm(nrow(df))
theta <- c(75, 0.1, -0.5, 2.5)
df$y <- model.matrix(~ x1 + x2 + z, df) %*% theta + eps

cmod_form <- y ~ x1 + x2
damod_form <- y ~ z - 1
msk <- df$z == 0

cmod <- lm(cmod_form, df[msk,])
des <- rct_design(z ~ unitid(uid), df)
damod <- as.lmitt(
  lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
)

a22inv <- .get_a22_inverse(damod)
b22 <- .get_b22(damod, type = "HC0", cadjust = FALSE)
b12 <- .get_b12(damod)
a21 <- .get_a21(damod)
a11inv <- .get_a11_inverse(damod)
b11 <- .get_b11(damod, type = "HC0", cadjust = FALSE)
propertee_meat <- (b22 - a21 %*% a11inv %*% b12 - t(b12) %*% a11inv %*% t(a21) +
                   a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))

## JOSH GEEX
estFun <- function(data){
  function(theta) {
    # covariance model eqns
    if (data$z == 0) {
      Xstar <- model.matrix(y ~ x1 + x2, data)
      cmod_eqns <- drop(data$y - Xstar %*% theta[1:3]) * Xstar
    } else {
      cmod_eqns <- rep(0, 3)
    }
    
    # itt model eqns (no intercept)
    X <- model.matrix(y ~ x1 + x2, data)
    Z <- model.matrix(y ~ z - 1, data)
    damod_eqns <- drop(data$weight * (data$y - X %*% theta[1:3] - Z * theta[4])) * Z
    
    out <- c(cmod_eqns, damod_eqns)
    return(out)
  }
}

geexRes <- geex::m_estimate(estFun,
                            data = cbind(df, "weight" = ate(des, data = df)),
                            root_control = geex::setup_root_control(start = rep(0.1,4)))
print("geex returns cov_hat(tau_hat):")
geexRes@vcov[4,4]
print("multiplying propertee internal matrices manually returns cov_hat(tau_hat):")
a22inv %*% propertee_meat %*% a22inv
print("propertee vcovDA function (multiplied by n_Q) returns cov_hat(tau_hat):")
vcovDA(damod, type = "CR0", cadjust = FALSE)
print("sandwich::sandwich (no clustering) returns cov_hat(tau_hat):")
sandwich::sandwich(damod, type = "HC0", cadjust = FALSE)

## CLUSTERED DATA
set.seed(50)
nc <- 4
mi <- 20
TT <- 2
df <- rbind(
  data.frame(
    "cid" = rep(seq_len(2), each = mi * TT),
    "uid" = rep(seq_len(mi), 2 * TT),
    "z" = rep(c(0, 1), each = mi * TT),
    "t" = rep(rep(seq_len(TT), each = mi), 2),
    "x1" = rep(c(rep(0, mi * 0.8), rep(1, mi * 0.2)), nc / 2 * TT),
    "x2" = rep(c(rep(0, mi * 0.5), rep(1, mi * 0.5)), nc / 2 * TT)),
  data.frame(
    "cid" = rep(seq_len(2) + 2, each = mi * (TT - 1)),
    "uid" = rep(seq_len(mi), 2 * (TT - 1)),
    "z" = rep(c(0, 1), each = mi * (TT - 1)),
    "t" = rep(rep(seq_len(TT - 1) + 1, each = mi), 2),
    "x1" = rep(c(rep(0, mi * 0.4), rep(1, mi * 0.6)), nc / 2 * (TT - 1)),
    "x2" = rep(c(rep(0, mi * 0.6), rep(1, mi * 0.4)), nc / 2 * (TT - 1)))
)

theta <- c(75, 7.5, -1, 2.5)
error_sd <- round(runif(nc, 1, 3), 1)
icc <- 0.2
eps <- rnorm(nrow(df))
Sigma <- matrix(0, nrow = nrow(df), ncol = nrow(df))
for (i in (seq_len(nc) - 1)) {
  msk <- df$cid == (i + 1)
  Sigma[msk, msk] <- diag(error_sd[i + 1] - icc, nrow = sum(msk)) + icc
}
A <- chol(Sigma)
eps <- t(A) %*% eps
df$y <- model.matrix(~ x1 + x2 + z, df) %*% theta + eps

cmod_form <- y ~ x1 + x2
damod_form <- y ~ z - 1
msk <- df$z == 0

cmod <- lm(cmod_form, df[msk,])
des <- rct_design(z ~ cluster(cid), df)
damod <- as.lmitt(
  lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
)

a22inv <- .get_a22_inverse(damod)
b22 <- .get_b22(damod, type = "HC0", cadjust = FALSE)
b12 <- .get_b12(damod)
a21 <- .get_a21(damod)
a11inv <- .get_a11_inverse(damod)
b11 <- .get_b11(damod, type = "HC0", cadjust = FALSE)
propertee_meat <- (b22 - a21 %*% a11inv %*% b12 - t(b12) %*% a11inv %*% t(a21) +
                   a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))

clusterEstFunc <- function(data){
  function(theta) {
    # covariance model eqns
    if (all(data$z == 0)) {
      Xstar <- model.matrix(y ~ x1 + x2, data)
      cmod_agg_func <- ifelse(dim(Xstar)[2] > 1, colSums, sum)
      cmod_eqns <- cmod_agg_func(drop(data$y - Xstar %*% theta[1:3]) * Xstar)
    } else {
      cmod_eqns <- rep(0, 3)
    }
    
    # itt model eqns (no intercept)
    X <- model.matrix(y ~ x1 + x2, data)
    Z <- model.matrix(y ~ z - 1, data)
    damod_agg_func <- ifelse(dim(Z)[2] > 1, colSums, sum)
    damod_eqns <- damod_agg_func(
      drop(data$weight * (data$y - X %*% theta[1:3] - Z * theta[4])) * Z)
    
    out <- c(cmod_eqns, damod_eqns)
    return(out)
  }
}

geexRes <- geex::m_estimate(clusterEstFunc,
                            data = cbind(df, "weight" = ate(des, data = df)),
                            units = "cid",
                            root_control = geex::setup_root_control(start = rep(0.1,4)))

print("geex returns cov_hat(tau_hat):")
geexRes@vcov[4,4]
print("multiplying propertee internal matrices manually returns cov_hat(tau_hat):")
a22inv %*% propertee_meat %*% a22inv
print("propertee vcovDA function (multiplied by n_Q) returns cov_hat(tau_hat):")
vcovDA(damod, type = "CR0", cadjust = FALSE)
print("sandwich::sandwich (no clustering) returns cov_hat(tau_hat):")
sandwich::sandwich(damod, meat. = sandwich::meatCL, cluster = df$cid,
                   type = "HC0", cadjust = FALSE)
