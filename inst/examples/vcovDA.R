### EQUAL SIZES OF CMOD AND DAMOD DATA, NO CLUSTERING, BALANCED RCT
set.seed(50)
n <- 1000
x1_beta <- 0.5
x2_beta <- -0.25
tau <- -2
error_sd <- 0.1

cmod_df <- data.frame(
  "uid" = NA_integer_,
  "x1" = rep(c(0, 1, 0, 1), each = n/4),
  "x2" = rep(c(0, 1), each = n/2)
)
cmod_df$y <- (cmod_df$x1 * x1_beta +
                cmod_df$x2 * x2_beta +
                error_sd * rep(c(-1, 1), n/2))
cmod_form <- ~ x1 + x2
X <- stats::model.matrix(cmod_form, stats::model.frame(cmod_form, cmod_df))
cmod <- lm(y ~ x1 + x2, cmod_df)

damod_df <- data.frame(
  "uid" = seq_len(n),
  "x1" = rep(c(0, 1, 0, 1), each = n/4),
  "x2" = rep(rep(c(0, 1), 4), each = n/8),
  "z" = rep(c(0, 1), each = n/2)
)
damod_df$y <- (damod_df$x1 * x1_beta +
                 damod_df$x2 * x2_beta +
                 damod_df$z * tau + error_sd * rep(c(-1, 1), n/2))
damod_form <- ~ z
C <- stats::model.matrix(damod_form, stats::model.frame(damod_form, damod_df))
Xstar <- stats::model.matrix(cmod_form, stats::model.frame(cmod_form, damod_df))

des <- rct_design(z ~ unitid(uid), damod_df)
m <- as.lmitt(
  lm(y ~ z, data = damod_df, weights = ate(des), offset = cov_adj(cmod))
)

## check blocks are what we expect
a11inv <- solve(crossprod(X))
expect_equal(a11inv, flexida:::.get_a11_inverse(m))
b22 <- 4 * error_sd^2 * crossprod(C) # weights = 2, residuals = error_sd for all
expect_equal(b22, flexida:::.get_b22(m, cadjust = FALSE, type = "HC0"))
a22inv <- solve(crossprod(C * m$weights, C))
expect_equal(a22inv, flexida:::.get_a22_inverse(m))
a21 <- -2 * crossprod(Xstar, C)
expect_equal(a21, flexida:::.get_a21(m))
b12 <- matrix(0, nrow = dim(X)[2], ncol = dim(C)[2])
expect_equal(b12, flexida:::.get_b12(m))
b11 <- error_sd^2 * crossprod(X) # residuals = error_sd for all units
expect_equal(b11, flexida:::.get_b11(m, cadjust = FALSE, type = "HC0"))

## breads are equal
expect_equal(a22inv, sandwich::bread(m) / n)

## meats are not! since there's no overlap of cmod and damod data, cov-adjusted
## meat terms simplifies to b22 + a21 %*% a11inv %*% b11 %*% a11inv %*% a21,
## which is a robust estimate of the cmod vcov matrix sandwiched by covariances
## of the terms in each model
# as shown below, we can compare to sandwich::meat when there's no clustering
expect_equal(sandwich::meat(m, adjust = FALSE),
             sandwich::meatCL(m, cadjust = FALSE, type = "HC0"))

expect_equal(b22, sandwich::meat(m, adjust = FALSE) * n)

# robust cmod vcov matrix is error_sd^2 * (X'X)^(-1)
expect_equal(a11inv %*% b11 %*% a11inv, error_sd^2 * a11inv)
expect_equal(a11inv %*% b11 %*% a11inv, sandwich::sandwich(cmod))

# accounting for cmod vcov should increase our vcov
expect_true(all((b22 + crossprod(a21, a11inv %*% b11 %*% a11inv) %*% a21) >
                  sandwich::meat(m, adjust = FALSE) * n))

## sandwiches are not equal! since meats are not
expect_equal(vcovDA(m, type = "HC0", cadjust = FALSE),
             a22inv %*% (b22 + crossprod(a21, a11inv %*% b11 %*% a11inv) %*%
                           a21) %*% a22inv)
expect_true(all(vcovDA(m, type = "HC0", cadjust = FALSE) >
                  sandwich::sandwich(m)))
