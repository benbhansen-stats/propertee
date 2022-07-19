test_that("issue 60", {
  data(simdata)

  # First, subset that does NOT exclude any clusters
  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata,
                    subset = 1:50 <= 47)
  expect_equal(nrow(structure(des)), 10)
  # By passing the subset
  mod1 <- lmitt(y ~ adopters(), data = simdata, design = des)
  expect_true("adopters()" %in% names(coef(mod1)))


  # Now, subset that excludes clusters
  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata,
                    subset = 1:50 <= 33)
  # This removes 3 complete clusters, and 1 row from a clsuter that remains,
  # leaving 33 rows of data with information about 7 clusters
  expect_equal(nrow(unique(simdata[1:50 <= 33, c("cid1", "cid2")])),
               nrow(structure(des)))

  # The model should use all rows related to those 7 clusters, which is 34.
  mod1 <- lmitt(y ~ adopters(), data = simdata, design = des)

  # This counts the number of rows in simdata whose clusters are found in the
  # design. (Instead of hardcoding the number, in case simdata changes)
  num_rows_with_clusters_in_design <-
    sum(apply(simdata[, c("cid1", "cid2")], 1, paste, collapse = "") %in%
          apply(units_of_assignment(des), 1, paste, collapse = ""))

  expect_equal(num_rows_with_clusters_in_design,
               nrow(mod1$model))

})


test_that("issue 60, with NAs", {
  data(simdata)
  simdata$cid1[1:10] <- simdata$cid2[1:10] <- NA_integer_
  des <- rct_design(z ~ uoa(cid1, cid2), data = simdata,
                    subset = !is.na(simdata$cid1) & !is.na(simdata$cid2))
  expect_silent(damod <- lmitt(y ~ x, data = simdata,
                               subset = !is.na(cid1) & !is.na(cid2),
                               design = des, weights = ate()))
  # Just testing that this doesn't error/warn
})
