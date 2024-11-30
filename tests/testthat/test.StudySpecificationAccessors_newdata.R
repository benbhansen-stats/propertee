test_that("treatments with `newdata`", {
  data(simdata)

  # Case 1: binary = FALSE
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid),
                    data = simdata)

  expect_true(all.equal(treatment(spec, newdata = simdata),
                        simdata[, "z", drop = FALSE]))

  # Case 2b: binary = TRUE, stored treatment is non-binary but has dichotomy
  spec <- rct_spec(o ~ cluster(uoa1, uoa2) + block(bid),
                    data = simdata)

  expect_equal(nrow(treatment(spec, dichotomy = o >= 3 ~ ., newdata = simdata)),
               50)
  expect_true(all(treatment(spec, dichotomy = o >= 3 ~ ., newdata = simdata)[,1] %in% 0:1))

})



test_that("units of assignment and clusters with `newdata`", {
  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all.equal(clusters(spec, newdata = simdata), simdata[, c("uoa1", "uoa2")]))

  clusters(spec) <- 1:10
  expect_length(clusters(spec, newdata = simdata)[, 1], 50)

  # This is only 5 since simdata doesn't have the same modification as two lines
  # above.
  expect_equal(length(unique(clusters(spec, newdata = simdata)[, 1])), 5)


  spec <- rd_spec(z ~ uoa(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all.equal(units_of_assignment(spec, newdata = simdata), simdata[, c("uoa1", "uoa2")]))

  units_of_assignment(spec) <- 1:10
  expect_length(units_of_assignment(spec, newdata = simdata)[, 1], 50)

  # This is only 5 since simdata doesn't have the same modification as two lines
  # above.
  expect_equal(length(unique(units_of_assignment(spec, newdata = simdata)[, 1])), 5)

})

test_that("blocks with `newdata`", {
  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all.equal(blocks(spec, newdata = simdata), simdata[, "bid", drop = FALSE]))

  blocks(spec) <- 1:10
  expect_length(blocks(spec, newdata = simdata)[, 1], 50)
  expect_equal(length(unique(blocks(spec, newdata = simdata)[, 1])), 10)

  # More blocks
  df <- data.frame(bid1 = 1:10,
                   bid2 = 10:1)
  blocks(spec) <- df
  expect_identical(dim(blocks(spec, newdata = simdata)), c(50L, 2L))

  ### multi-dimensional blocks

  simdata2 <- simdata
  simdata2$uoa1 <- 1:50
  simdata2 <- rbind(simdata2, simdata2) # Duplicate rows to make data different
  names(simdata2)[1:2] <- c("cid", "bida")
  spec <- rd_spec(z ~ cluster(cid) + block(bid, bida) + forcing(force),
                   data = simdata2)

  expect_identical(dim(blocks(spec, newdata = simdata2)), c(100L, 2L))
})

test_that("forcing with `newdata`", {
  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all.equal(forcings(spec, newdata = simdata), simdata[, "force", drop = FALSE]))

  forcings(spec) <- 1:10
  expect_length(forcings(spec, newdata = simdata)[, 1], 50)
  expect_equal(length(unique(forcings(spec, newdata = simdata)[, 1])), 10)

  # More forcings
  df <- data.frame(f1 = 1:10,
                   f2 = 10:1)
  forcings(spec) <- df
  expect_identical(dim(forcings(spec, newdata = simdata)), c(50L, 2L))

  ### multi-dimensional forcings

  simdata2 <- simdata
  simdata2$force2 <- rep(rnorm(10, mean = 5), times = rep(c(4, 6), times = 5))
  simdata2 <- rbind(simdata2, simdata2) # Duplicate rows to make data different
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) +
                     forcing(force, force2), data = simdata2)


  expect_identical(dim(forcings(spec, newdata = simdata2)), c(100L, 2L))
})

test_that(".specification_accessors_newdata_validate", {
  expect_warning(.specification_accessors_newdata_validate(1, NULL))

  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  names(simdata)[names(simdata) == "uoa1"] <- "ccc"
  names(simdata)[names(simdata) == "force"] <- "fff"

  expect_true(.specification_accessors_newdata_validate(simdata,
                                                 by = list(force = "fff",
                                                           uoa1 = "ccc")))

})

test_that(".get_col_from_new_data", {
  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  simdata2 <- simdata
  names(simdata2)[names(simdata2) == "uoa1"] <- "ccc"
  names(simdata2)[names(simdata2) == "force"] <- "fff"

  expect_identical(.get_col_from_new_data(spec, simdata, type = "b"),
                   .get_col_from_new_data(spec, simdata2, type = "b",
                                          by = list(force = "fff",
                                                    uoa1 = "ccc")))


})
