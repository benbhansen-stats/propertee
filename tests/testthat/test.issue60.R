test_that("issue 60", {
  data(simdata)

  # First, subset that does NOT exclude any clusters
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata,
                    subset = 1:50 <= 47)
  expect_equal(nrow(get_structure(spec)), 10)
  # By passing the subset
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec)


  # Now, subset that excludes clusters
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata,
                    subset = 1:50 <= 33)
  # This removes 3 complete clusters, and 1 row from a clsuter that remains,
  # leaving 33 rows of data with information about 7 clusters
  expect_equal(nrow(unique(simdata[1:50 <= 33, c("uoa1", "uoa2")])),
               nrow(get_structure(spec)))

  # The model should use all rows related to those 7 clusters, which is 34.
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec)

  # This counts the number of rows in simdata whose clusters are found in the
  # specification. (Instead of hardcoding the number, in case simdata changes)
  num_rows_with_clusters_in_spec <-
    sum(apply(simdata[, c("uoa1", "uoa2")], 1, paste, collapse = "") %in%
          apply(units_of_assignment(spec), 1, paste, collapse = ""))

  expect_equal(num_rows_with_clusters_in_spec,
               nrow(mod1$model))

})


test_that("issue 60, with NAs", {
  data(simdata)
  simdata$uoa1[1:10] <- simdata$uoa2[1:10] <- NA_integer_
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata,
                    subset = !is.na(simdata$uoa1) & !is.na(simdata$uoa2))
  expect_silent(damod <- lmitt(y ~ x, data = simdata,
                               subset = !is.na(uoa1) & !is.na(uoa2),
                               specification = spec, weights = ate()))
  # Just testing that this doesn't error/warn
})
