if (requireNamespace("multcomp", quietly = TRUE)) {
  test_that(".make_uoa_ids for mmm object", {
    set.seed(39)
    fct_df <- data.frame(id = seq_len(13),
                         a = factor(c(rep(0, 3), rep(1, 6), rep(2, 4))),
                         x = rnorm(13),
                         y = rnorm(13))
    ca_df <- data.frame(id = seq(14, 26), x = rnorm(13), y = rnorm(13))
    
    ca <- lm(y ~ x, ca_df)
    spec <- rct_spec(a ~ unitid(id), fct_df)
    
    a1 <- lmitt(y ~ 1, spec, fct_df, dichotomy = a == 1 ~ a == 0)
    a2 <- lmitt(y ~ 1, spec, fct_df, dichotomy = a == 2 ~ a == 0)
    a1lm <- lmitt(lm(y ~ assigned(spec, fct_df, a == 1 ~ a == 0), fct_df), spec)
    a2lm <- lmitt(lm(y ~ assigned(spec, fct_df, a == 2 ~ a == 0), fct_df), spec)
    
    mods <- multcomp::mmm(a1, a2)
    expect_true(inherits(ids <- .make_uoa_ids(mods, "CR0"), "factor"))
    expect_equal(length(ids), 13)
    
    mods <- multcomp::mmm(a1lm, a2lm)
    expect_true(inherits(ids <- .make_uoa_ids(mods, "CR0"), "factor"))
    expect_equal(length(ids), 13)
    
    a1 <- lmitt(y ~ 1, spec, fct_df, offset = cov_adj(ca, newdata = fct_df), dichotomy = a == 1 ~ a == 0)
    a2 <- lmitt(y ~ 1, spec, fct_df, offset = cov_adj(ca, newdata = fct_df), dichotomy = a == 2 ~ a == 0)
    a1lm <- lmitt(
      lm(y ~ assigned(spec, fct_df, a == 1 ~ a == 0), fct_df,
         off = cov_adj(ca, newdata = fct_df, sp = spec)),
      spec)
    a2lm <- lmitt(
      lm(y ~ assigned(spec, fct_df, a == 2 ~ a == 0), fct_df,
         off = cov_adj(ca, newdata = fct_df, sp = spec)),
      spec)
    
    mods <- multcomp::mmm(a1, a2)
    expect_true(inherits(ids <- .make_uoa_ids(mods, "CR0"), "factor"))
    expect_equal(length(ids), 26)
    
    mods <- multcomp::mmm(a1lm, a2lm)
    expect_true(inherits(ids <- .make_uoa_ids(mods, "CR0"), "factor"))
    expect_equal(length(ids), 26)
  })
  
  test_that("factor treatment variables (#78) no weights, no cov adj", {
    set.seed(39)
    fct_df <- data.frame(id = seq_len(13),
                         a = factor(c(rep(0, 3), rep(1, 6), rep(2, 4))),
                         y = rnorm(13))
    spec <- rct_spec(a ~ unitid(id), fct_df)
    
    a1 <- lmitt(y ~ 1, spec, fct_df, dichotomy = a == 1 ~ a == 0)
    a2 <- lmitt(y ~ 1, spec, fct_df, dichotomy = a == 2 ~ a == 0)
    a1lm <- lmitt(lm(y ~ assigned(spec, fct_df, a == 1 ~ a == 0), fct_df), spec)
    a2lm <- lmitt(lm(y ~ assigned(spec, fct_df, a == 2 ~ a == 0), fct_df), spec)
    
    # check pt estimates of dichotomized mods
    grp_means <- rowsum(fct_df$y, fct_df$a) / rowsum(numeric(13) + 1, fct_df$a)
    expect_true(all.equal(a1$coefficients[2], grp_means[2,] - grp_means[1,], check.attributes = FALSE))
    expect_true(all.equal(a1$coefficients[2], a1lm$coefficients[2], check.attributes = FALSE))
    expect_true(all.equal(a2$coefficients[2], grp_means[3,] - grp_means[1,], check.attributes = FALSE))
    expect_true(all.equal(a2$coefficients[2], a2lm$coefficients[2], check.attributes = FALSE))
    
    # check vcov matrix for lmitt.formula
    mods <- multcomp::mmm(a1, a2)
    n <- nrow(fct_df)
    tsts <- multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1, 0), nrow = 1)),
                           vcov = function(x) vcov_tee(x, type = "HC0"))
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      n / (n-1) * (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
    
    # check vcov matrix for lmitt.lm
    mods <- multcomp::mmm(a1lm, a2lm)
    tsts <- multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1, 0), nrow = 1)),
                           vcov = function(x) vcov_tee(x, type = "HC0"))
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      n / (n-1) * (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
  })
  
  test_that("factor treatment variables (#78) weights, no cov adj", {
    set.seed(39)
    fct_df <- data.frame(id = seq_len(13),
                         a = factor(c(rep(0, 3), rep(1, 6), rep(2, 4))),
                         y = rnorm(13))
    spec <- rct_spec(a ~ unitid(id), fct_df)
    
    a1 <- lmitt(y ~ 1, spec, fct_df, w = ate(spec), dichotomy = a == 1 ~ a == 0)
    a2 <- lmitt(y ~ 1, spec, fct_df, w = ate(spec), dichotomy = a == 2 ~ a == 0)
    a1lm <- lmitt(
      lm(y ~ assigned(spec, fct_df, a == 1 ~ a == 0), fct_df,
         w = ate(spec, dichotomy = a == 1 ~ a == 0, data = fct_df)),
      spec)
    a2lm <- lmitt(
      lm(y ~ assigned(spec, fct_df, a == 2 ~ a == 0), fct_df,
         w = ate(spec, dichotomy = a == 2 ~ a == 0, data = fct_df)),
      spec)
    
    # check pt estimates of dichotomized mods
    grp_means_d1 <- rowsum(a1$model$y * a1$weights, a1$model$a.) / rowsum(a1$weights, a1$model$a.)
    grp_means_d2 <- rowsum(a2$model$y * a2$weights, a2$model$a.) / rowsum(a2$weights, a2$model$a.)
    expect_true(all.equal(a1$coefficients[2], drop(c(-1, 1) %*% grp_means_d1), check.attributes = FALSE))
    expect_true(all.equal(a1$coefficients[2], a1lm$coefficients[2], check.attributes = FALSE))
    expect_true(all.equal(a2$coefficients[2], drop(c(-1, 1) %*% grp_means_d2), check.attributes = FALSE))
    expect_true(all.equal(a2$coefficients[2], a2lm$coefficients[2], check.attributes = FALSE))
    
    # check vcov matrix for lmitt.formula
    mods <- multcomp::mmm(a1, a2)
    n <- nrow(fct_df)
    tsts <- multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1, 0), nrow = 1)),
                           vcov = function(x) vcov_tee(x, type = "HC0"))
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      n / (n-1) * (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
    
    # check vcov matrix for lmitt.lm
    mods <- multcomp::mmm(a1lm, a2lm)
    tsts <- multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1, 0), nrow = 1)),
                           vcov = function(x) vcov_tee(x, type = "HC0"))
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      n / (n-1) * (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
  })
  
  test_that("factor treatment variables (#78) no weights, cov adj", {
    set.seed(39)
    fct_df <- data.frame(id = seq_len(13),
                         a = factor(c(rep(0, 3), rep(1, 6), rep(2, 4))),
                         x = rnorm(13),
                         y = rnorm(13))
    ca_df <- data.frame(id = seq(14, 26), x = rnorm(13), y = rnorm(13))
    spec <- rct_spec(a ~ unitid(id), fct_df)
    
    ca <- lm(y ~ x, ca_df)
    a1 <- lmitt(y ~ 1, spec, fct_df, offset = cov_adj(ca, newdata = fct_df), dichotomy = a == 1 ~ a == 0)
    a2 <- lmitt(y ~ 1, spec, fct_df, offset = cov_adj(ca, newdata = fct_df), dichotomy = a == 2 ~ a == 0)
    a1lm <- lmitt(
      lm(y ~ assigned(spec, fct_df, a == 1 ~ a == 0), fct_df,
         off = cov_adj(ca, newdata = fct_df, sp = spec)),
      spec)
    a2lm <- lmitt(
      lm(y ~ assigned(spec, fct_df, a == 2 ~ a == 0), fct_df,
         off = cov_adj(ca, newdata = fct_df, sp = spec)),
      spec)
    
    # check pt estimates of dichotomized mods
    grp_means <- rowsum(fct_df$y - cov_adj(ca, newdata = fct_df, sp = spec), fct_df$a) / rowsum(numeric(13) + 1, fct_df$a)
    expect_true(all.equal(a1$coefficients[2], grp_means[2,] - grp_means[1,], check.attributes = FALSE))
    expect_true(all.equal(a1$coefficients[2], a1lm$coefficients[2], check.attributes = FALSE))
    expect_true(all.equal(a2$coefficients[2], grp_means[3,] - grp_means[1,], check.attributes = FALSE))
    expect_true(all.equal(a2$coefficients[2], a2lm$coefficients[2], check.attributes = FALSE))
    
    # check vcov matrix for lmitt.formula
    mods <- multcomp::mmm(a1, a2)
    n <- nrow(fct_df) + nrow(ca_df)
    tsts <- multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1, 0, 0), nrow = 1)),
                           vcov = function(x) vcov_tee(x, type = "HC0", cov_adj_rcorrect = "HC0"))
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      n / (n-1) * (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
    
    # check vcov matrix for lmitt.lm
    mods <- multcomp::mmm(a1lm, a2lm)
    tsts <- multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1, 0, 0), nrow = 1)),
                           vcov = function(x) vcov_tee(x, type = "HC0", cov_adj_rcorrect = "HC0"))
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      n / (n-1) * (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
  })
  
  test_that("combine lmitts (#79) for different DV's", {
    set.seed(39)
    Q_df <- data.frame(id = seq_len(13),
                       a = rbinom(13, 1, 0.5),
                       x = rnorm(13),
                       y1 = rnorm(13),
                       y2 = rnorm(13))
    ca_df <- data.frame(id = seq(14, 26), x = rnorm(13), y1 = rnorm(13), y2 = rnorm(13))
    spec <- rct_spec(a ~ unitid(id), Q_df)
    
    ca1 <- lm(y1 ~ x, ca_df)
    ca2 <- lm(y2 ~ x, ca_df)
    a1 <- lmitt(y1 ~ 1, spec, Q_df, offset = cov_adj(ca1, newdata = Q_df))
    a2 <- lmitt(y2 ~ 1, spec, Q_df, offset = cov_adj(ca2, newdata = Q_df))
    a1lm <- lmitt(
      lm(y1 ~ assigned(spec, Q_df), Q_df,
         off = cov_adj(ca1, newdata = Q_df, sp = spec)),
      spec)
    a2lm <- lmitt(
      lm(y2 ~ assigned(spec, Q_df), Q_df,
         off = cov_adj(ca2, newdata = Q_df, sp = spec)),
      spec)
    
    # check vcov matrix for lmitt.formula
    mods <- multcomp::mmm(a1, a2)
    n <- nrow(Q_df) + nrow(ca_df)
    tsts <- multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1, 0, 0), nrow = 1)),
                           vcov = function(x) vcov_tee(x, type = "HC0", cov_adj_rcorrect = "HC0"))
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      n / (n-1) * (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
    
    # check vcov matrix for lmitt.lm
    mods <- multcomp::mmm(a1lm, a2lm)
    tsts <- multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1, 0, 0), nrow = 1)),
                           vcov = function(x) vcov_tee(x, type = "HC0", cov_adj_rcorrect = "HC0"))
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      n / (n-1) * (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
  })
  
  test_that("combine lmitts (#79) for different moderator variables", {
    set.seed(39)
    Q_df <- data.frame(id = seq_len(45),
                       b = rep(LETTERS[1:3], each = 15),
                       a = rep(c(rep(0, 8), rep(1, 7)), 3),
                       x1 = factor(sample(seq_len(3), 45, re = TRUE)),
                       x2 = factor(sample(letters[1:2], 45, re = TRUE)),
                       y = rnorm(45))
    ca_df <- data.frame(id = seq(46, 58),
                        x1 = factor(sample(seq_len(3), 13, re = TRUE)),
                        x2 = factor(sample(letters[1:2], 13, re = TRUE)),
                        y = rnorm(13))
    spec <- rct_spec(a ~ unitid(id) + block(b), Q_df)
    
    ca <- lm(y ~ x1 + x2, ca_df)
    a1 <- lmitt(y ~ x1, spec, Q_df, w = ate(spec), offset = cov_adj(ca, newdata = Q_df))
    a2 <- lmitt(y ~ x2, spec, Q_df, w = ate(spec), offset = cov_adj(ca, newdata = Q_df))
    
    # check vcov matrix for lmitt.formula
    mods <- multcomp::mmm(a1, a2)
    n <- nrow(Q_df) + nrow(ca_df)
    tsts <- multcomp::glht(
      mods, rbind(c(rep(0, 3), 1, rep(0, 16)), c(rep(0, 14), 1, rep(0, 5))),
      vcov. = function(x, ...) vcov_tee(x, type = "HC0", cov_adj_rcorrect = "HC0", ...), cadjust = FALSE)
    vmat <- summary(tsts)$vcov
    expect_true(all.equal(
      vmat,
      (1 / n) * (bread(mods) %*% (crossprod(estfun(mods)) / n) %*% bread(mods)),
      check.attributes = FALSE
    ))
  })
  
  test_that("error if Q or C df's don't match across models in mmm object", {
    on.exit(options("show.error.messages" = TRUE))
    options("show.error.messages" = FALSE)
    set.seed(39)
    Q_df1 <- data.frame(id = seq_len(45),
                        b = rep(LETTERS[1:3], each = 15),
                        a = rep(c(rep(0, 8), rep(1, 7)), 3),
                        x1 = factor(sample(seq_len(3), 45, re = TRUE)),
                        y1 = rnorm(45),
                        y2 = rnorm(45))
    Q_df2 <- data.frame(id = seq_len(60),
                        b = rep(LETTERS[1:3], each = 20),
                        a = rep(c(rep(0, 8), rep(1, 12)), 3),
                        x1 = factor(sample(seq_len(3), 60, re = TRUE)),
                        y1 = rnorm(60),
                        y2 = rnorm(60))
    ca_df1 <- data.frame(id = seq(46, 58),
                         x1 = factor(sample(seq_len(3), 13, re = TRUE)),
                         y1 = rnorm(13),
                         y2 = rnorm(13))
    ca_df2 <- data.frame(id = seq(46, 59),
                         x1 = factor(sample(seq_len(3), 14, re = TRUE)),
                         y1 = rnorm(14),
                         y2 = rnorm(14))
    spec1 <- rct_spec(a ~ unitid(id) + block(b), Q_df1)
    spec2 <- rct_spec(a ~ unitid(id) + block(b), Q_df2)
    
    a1 <- lmitt(y1 ~ 1, spec1, Q_df1, w = ate(spec1))
    a2 <- lmitt(y2 ~ 1, spec2, Q_df2, w = ate(spec2))
    mods <- multcomp::mmm(a1, a2)
    expect_error(
      multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1), nrow = 1)), vcov = vcov_tee, cadjust = FALSE)
    )
    
    ca1 <- lm(y1 ~ x1, ca_df1)
    ca2 <- lm(y2 ~ x1, ca_df2)
    a1 <- lmitt(y1 ~ 1, spec1, Q_df1, w = ate(spec1), offset = cov_adj(ca1))
    a2 <- lmitt(y2 ~ 1, spec1, Q_df1, w = ate(spec1), offset = cov_adj(ca2))
    mods <- multcomp::mmm(a1, a2)
    expect_error(
      multcomp::glht(mods, multcomp::mlf(matrix(c(0, 1), nrow = 1)), vcov = vcov_tee, cadjust = FALSE)
    )
  })
}
