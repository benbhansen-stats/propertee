if (requireNamespace("robustbase", quietly = TRUE)) {
  test_that("the columns of estfun.glmrob have approximately 0 mean",{
    data(lsoSynth)

    mod <- robustbase::glmrob(probation_year1~age_at_entry+lhsgrade_pct+male,
                  family=binomial,data=lsoSynth)


    ee <- estfun(mod)

    expect_true(all(abs(colMeans(ee))<0.001))
  })

  test_that("the sandwich vcov matrix for a glmrob object is approximately the same as the cov returned with the object",{
    data(lsoSynth)

    mod <- robustbase::glmrob(probation_year1~age_at_entry+lhsgrade_pct+male,
                  family=binomial,data=lsoSynth)

    expect_true(max(abs(sandwich::sandwich(mod)-mod$cov))<0.01)
  })
} else {
  expect_true(TRUE)
  # Dummy test to avoid empty test warnings
}
