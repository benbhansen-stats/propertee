test_that("treatments with `newdata`", {
  data(simdata)

  # Case 1: binary = FALSE
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid),
                    data = simdata)

  expect_true(all.equal(treatment(des, newdata = simdata),
                        simdata[, "z", drop = FALSE]))

  # Case 2a: binary = TRUE, treatment is stored as binary
  expect_true(all.equal(treatment(des, newdata = simdata, binary = TRUE),
                        simdata[, "z", drop = FALSE]))


  # Case 2b: binary = TRUE, stored treatment is non-binary but has dichotomy
  des <- rct_design(o ~ cluster(uoa1, uoa2) + block(bid),
                    data = simdata)

  expect_equal(nrow(treatment(des, dichotomy = o >= 3 ~ ., newdata = simdata, binary = TRUE)),
               50)
  expect_true(all(treatment(des, dichotomy = o >= 3 ~ ., newdata = simdata, binary = TRUE)[,1] %in% 0:1))

  # Case 3a: binary = "ifany", no dichotomy
  des <- rct_design(o ~ cluster(uoa1, uoa2) + block(bid),
                    data = simdata)

  expect_true(all.equal(treatment(des, newdata = simdata, binary = "ifany"),
                        simdata[, "o", drop = FALSE]))

  # Case 3b: binary = "ifany", dichotomy
  des <- rct_design(o ~ cluster(uoa1, uoa2) + block(bid),
                    data = simdata)

  expect_equal(nrow(treatment(des, dichotomy = o >= 3 ~ ., newdata = simdata, binary = "ifany")),
               50)
  expect_true(all(treatment(des, dichotomy = o >= 3 ~ ., newdata = simdata, binary = "ifany")[,1] %in% 0:1))


})



test_that("units of assignment and clusters with `newdata`", {
  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all.equal(clusters(des, newdata = simdata), simdata[, c("uoa1", "uoa2")]))

  clusters(des) <- 1:10
  expect_length(clusters(des, newdata = simdata)[, 1], 50)

  # This is only 5 since simdata doesn't have the same modification as two lines
  # above.
  expect_equal(length(unique(clusters(des, newdata = simdata)[, 1])), 5)


  des <- rd_design(z ~ uoa(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all.equal(units_of_assignment(des, newdata = simdata), simdata[, c("uoa1", "uoa2")]))

  units_of_assignment(des) <- 1:10
  expect_length(units_of_assignment(des, newdata = simdata)[, 1], 50)

  # This is only 5 since simdata doesn't have the same modification as two lines
  # above.
  expect_equal(length(unique(units_of_assignment(des, newdata = simdata)[, 1])), 5)

})

test_that("blocks with `newdata`", {
  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all.equal(blocks(des, newdata = simdata), simdata[, "bid", drop = FALSE]))

  blocks(des) <- 1:10
  expect_length(blocks(des, newdata = simdata)[, 1], 50)
  expect_equal(length(unique(blocks(des, newdata = simdata)[, 1])), 10)

  # More blocks
  df <- data.frame(bid1 = 1:10,
                   bid2 = 10:1)
  blocks(des) <- df
  expect_identical(dim(blocks(des, newdata = simdata)), c(50L, 2L))

  ### multi-dimensional blocks

  simdata2 <- simdata
  simdata2$uoa1 <- 1:50
  simdata2 <- rbind(simdata2, simdata2) # Duplicate rows to make data different
  names(simdata2)[1:2] <- c("cid", "bida")
  des <- rd_design(z ~ cluster(cid) + block(bid, bida) + forcing(force),
                   data = simdata2)

  expect_identical(dim(blocks(des, newdata = simdata2)), c(100L, 2L))
})

test_that("forcing with `newdata`", {
  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all.equal(forcings(des, newdata = simdata), simdata[, "force", drop = FALSE]))

  forcings(des) <- 1:10
  expect_length(forcings(des, newdata = simdata)[, 1], 50)
  expect_equal(length(unique(forcings(des, newdata = simdata)[, 1])), 10)

  # More forcings
  df <- data.frame(f1 = 1:10,
                   f2 = 10:1)
  forcings(des) <- df
  expect_identical(dim(forcings(des, newdata = simdata)), c(50L, 2L))

  ### multi-dimensional forcings

  simdata2 <- simdata
  simdata2$force2 <- rep(rnorm(10, mean = 5), times = rep(c(4, 6), times = 5))
  simdata2 <- rbind(simdata2, simdata2) # Duplicate rows to make data different
  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) +
                     forcing(force, force2), data = simdata2)


  expect_identical(dim(forcings(des, newdata = simdata2)), c(100L, 2L))
})

test_that(".design_accessors_newdata_validate", {
  expect_warning(.design_accessors_newdata_validate(1, NULL))

  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  names(simdata)[names(simdata) == "uoa1"] <- "ccc"
  names(simdata)[names(simdata) == "force"] <- "fff"

  expect_true(.design_accessors_newdata_validate(simdata,
                                                 by = list(force = "fff",
                                                           uoa1 = "ccc")))

})

test_that(".get_col_from_new_data", {
  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  simdata2 <- simdata
  names(simdata2)[names(simdata2) == "uoa1"] <- "ccc"
  names(simdata2)[names(simdata2) == "force"] <- "fff"

  expect_identical(.get_col_from_new_data(des, simdata, type = "b"),
                   .get_col_from_new_data(des, simdata2, type = "b",
                                          by = list(force = "fff",
                                                    uoa1 = "ccc")))


})
