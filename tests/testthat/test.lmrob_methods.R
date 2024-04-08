if (requireNamespace("robustbase", quietly = TRUE)) {

  test_that("invalid bread. and estfun.lmrob calls produce expected errors", {
    data(simdata)

    cmod <- robustbase::lmrob(y ~ x, simdata)
    cmod$control$method <- "invalid_method"
    expect_error(estfun(cmod), "supports only SM or MM")
    expect_error(bread(cmod), "supports only SM or MM")

    cmod$control$method <- "SM"
    cmod$control$psi <- NULL
    expect_error(estfun(cmod), "psi is not defined")
    expect_error(bread(cmod), "psi is not defined")

    cmod$control <- NULL
    expect_error(estfun(cmod), "object must have a `control` element")
    expect_error(bread(cmod), "object must have a `control` element")
  })

  test_that("estfun.lmrob produces expected output", {
    data(simdata)

    cmod <- robustbase::lmrob(y ~ x, simdata)
    ef <- estfun(cmod)

    expect_equal(dim(ef), c(nrow(simdata), 3))
    expect_equal(colnames(ef), c("sigma", colnames(stats::model.matrix(cmod))))
    expect_true(is.null(attr(ef, "assign")))
    expect_true(is.null(attr(ef, "contrasts")))
  })

  test_that("bread.lmrob produces expected output", {
    data(simdata)

    cmod <- robustbase::lmrob(y ~ x, simdata)
    ainv <- bread(cmod)

    expect_equal(dim(ainv), c(2, 3))
    expect_equal(dimnames(ainv),
                 list(colnames(stats::model.matrix(cmod)),
                      c("sigma", colnames(stats::model.matrix(cmod)))))
  })

}
