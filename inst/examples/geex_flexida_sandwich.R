## NOT CLUSTERED DATA
set.seed(980)
nc <- 2
mi <- 20
df <- data.frame(
  # "cid" = rep(seq_len(2), each = mi * TT),
  "uid" = seq_len(mi * nc),#rep(seq_len(mi), 2 * TT),
  "z" = rep(c(0, 1), each = mi),
  "x1" = rnorm(sd = 2, n = mi * nc),#rep(c(rep(0, mi * 0.8), rep(1, mi * 0.2)), nc / 2 * TT),
  "x2" = rep(c(rep(0, mi * 0.5), rep(1, mi * 0.5)), nc / 2)
)
eps <- rnorm(nrow(df))
theta <- c(75, 0.1, -0.5, 2.5)
df$y <- model.matrix(~ x1 + x2 + z, df) %*% theta + eps

cmod_form <- y ~ x1 + x2
damod_form <- y ~ z - 1
msk <- df$z == 0

cmod <- lm(cmod_form,rbind(
  df[msk,]
  # , matrix(rep(c(NA_integer_, NA_integer_, 0.4, 82), 100),
  #        byrow = TRUE, nrow = 100,
  #        dimnames = list(seq_len(100), colnames(df))),
  # , c(NA_integer_, NA_integer_, 0.4, 82)
  # , c(NA_integer_, NA_integer_, 0.9, 89)
  # , c(NA_integer_, NA_integer_, -0.5, 72)
  # , c(NA_integer_, NA_integer_, -0.5, 72)
))
# cmod <- lm(cmod_form, df[msk,])
des <- rct_design(z ~ unitid(uid), df)
damod <- as.lmitt(
  lm(damod_form, data = df, weights = ate(des), offset = cov_adj(cmod))
)

sandwich::sandwich(damod, type = "HC0", cadjust = FALSE)
a22inv <- .get_a22_inverse(damod)
b22 <- .get_b22(damod, type = "HC0", cadjust = FALSE)
b12 <- .get_b12(damod)
a21 <- .get_a21(damod)
a11inv <- .get_a11_inverse(damod)
b11 <- .get_b11(damod, type = "HC0", cadjust = FALSE)
flexida_meat <- (b22 - a21 %*% a11inv %*% b12 - t(b12) %*% a11inv %*% t(a21) +
                   a21 %*% a11inv %*% b11 %*% a11inv %*% t(a21))

## JOSH GEEX
estFun <- function(data){
  function(theta) {
    # covariance model eqns (no intercept)
    if (data$z == 0) {
      Xstar <- model.matrix(y ~ x1 + x2, data)
      cmod_eqns <- drop(data$y - Xstar %*% theta[1:3]) * Xstar
    } else {
      cmod_eqns <- rep(0, 3)
    }
    
    # itt model eqns (include overall intercept)
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
print("multiplying flexida internal matrices manually returns cov_hat(tau_hat):")
a22inv %*% flexida_meat %*% a22inv
print("flexida vcovDA function (multiplied by n_Q) returns cov_hat(tau_hat):")
vcovDA(damod, type = "HC0", cadjust = FALSE) * dim(model.matrix(damod))[1]
print("sandwich::sandwich (no clustering) returns cov_hat(tau_hat):")
sandwich::sandwich(damod, type = "HC0", cadjust = FALSE)

