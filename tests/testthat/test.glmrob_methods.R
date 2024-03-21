if (requireNamespace("robustbase", quietly = TRUE)) {

  test_that("basic tests just to touch glmrobMethods.R file", {
    data(simdata)

    cmod <- robustbase::glmrob(z ~ x, simdata, family = binomial)
    expect_true(is(estfun(cmod), "matrix"))
    expect_true(is(bread(cmod), "matrix"))
    cmod <- robustbase::glmrob(bid ~ x, simdata, family = poisson)
    expect_true(is(estfun(cmod), "matrix"))
    expect_true(is(bread(cmod), "matrix"))
    cmod <- robustbase::glmrob(y ~ x, simdata, family = gaussian)
    expect_true(is(estfun(cmod), "matrix"))
    expect_true(is(bread(cmod), "matrix"))
    cmod <- robustbase::glmrob(force ~ x, simdata, family = Gamma)
    expect_true(is(estfun(cmod), "matrix"))
    expect_true(is(bread(cmod), "matrix"))
  })
}
