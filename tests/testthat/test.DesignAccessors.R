test_that("Accessing and replacing treatment", {
  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_equal(dim(treatment(des)), c(10, 1))
  expect_true(is.numeric(treatment(des)[, 1]))
  expect_identical(treatment(des), des@structure[, 1, drop = FALSE])

  treatment(des) <- rep(0:1, each = 5)
  expect_equal(treatment(des), data.frame(z = rep(0:1, each = 5)))

  df <- data.frame(z = rep(0:1, each = 5))
  treatment(des) <- df
  expect_equal(treatment(des), data.frame(z = rep(0:1, each = 5)))

  treatment(des)[1:2, 1] <- 1
  expect_equal(treatment(des),
               data.frame(z = rep(c(1, 0, 1),
                                  times = c(2, 3, 5))))

  expect_error(treatment(des) <- 1:5,
               "same number")

  expect_error(treatment(des) <- data.frame(a = c(1, 0, 1, 0, 1)),
               "same number")

  # Continuous treatment, no dichotomy
  des <- rd_design(dose ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(treatment(des), des@structure[, 1, drop = FALSE])

  # Continuous treatment, dichotomy
  des <- rd_design(dose ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(treatment(des), des@structure[, 1, drop = FALSE])
  tt <- treatment(des, dichotomy = dose <= 100 ~ dose == 200)

  expect_equal(dim(tt), c(10, 1))
  expect_true(is.numeric(tt[, 1]))
  expect_true(all(tt[, 1] %in% c(0, 1, NA)))
  expect_equal(colnames(tt), "z__")

  ############ .bin_txt

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(treatment(des)[, 1], .bin_txt(des))

  des <- rd_design(o ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all(.bin_txt(des, dichotomy = o >= 3 ~ .) %in% c(0, 1, NA)))
  expect_identical(.bin_txt(des, dichotomy = o >= 3 ~ .),
                   treatment(des, dichotomy = o >= 3 ~ .)[, 1])

  des <- rd_design(o ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_error(.bin_txt(des), "Must provide a dichotomy")

})

test_that("New .bin_txt, .expand_txt, and .apply_dichotomy errors", {
  des <- rct_design(z ~ unitid(uoa1, uoa2), simdata)
  expect_error(.bin_txt(des, data.frame("uoa1" = seq(100, 110))),
               "Not all unit of assignment variables")

  expect_error(.apply_dichotomy(simdata$dose, dose > 50 ~ dose == 50),
               "expected to be a named `data.frame`")
  expect_error(.apply_dichotomy(treatment(des), dichotomy = dose > 50 ~ dose == 50),
               "Could not find variables specified in dichotomy")
  expect_error(.apply_dichotomy(simdata, . ~ .),
               "At least one side")
  expect_error(.apply_dichotomy(simdata, dichotomy = dose >= 50 ~ dose == 50),
               "dichotomy overlaps")
})

test_that("Accessing and replacing unit of assignment", {
  data(simdata)

  des <- rd_design(z ~ unit_of_assignment(uoa1, uoa2) + block(bid) +
                     forcing(force),
                   data = simdata)

  expect_identical(units_of_assignment(des), des@structure[, 2:3])

  units_of_assignment(des) <- data.frame(uoa1 = 10:1, uoa2 = 1:10)
  expect_identical(var_names(des, "u"), c("uoa1", "uoa2"))
  expect_true(all(data.frame(uoa1 = 10:1, uoa2 = 1:10) ==
                    units_of_assignment(des)))

  units_of_assignment(des) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(var_names(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) ==
                    units_of_assignment(des)))

  m <- matrix(1:20, ncol = 2)
  units_of_assignment(des) <- m
  expect_identical(var_names(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) ==
                    units_of_assignment(des)))

  colnames(m) <- c("qwe", "asd")
  units_of_assignment(des) <- m
  expect_identical(var_names(des, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) ==
                    units_of_assignment(des)))

  units_of_assignment(des)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) ==
                    units_of_assignment(des)))

  # less units_of_assignment
  des2 <- des
  units_of_assignment(des2) <- m[, 1]
  expect_identical(var_names(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == units_of_assignment(des2)))

  des2 <- des
  units_of_assignment(des2) <- 10:1
  expect_identical(var_names(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == units_of_assignment(des2)))

  df <- data.frame(abc = 3:12)
  des2 <- des
  units_of_assignment(des2) <- df
  expect_identical(var_names(des2, "u"), colnames(df))
  expect_true(all(df == units_of_assignment(des2)))

  ########## UnitsOfAssignment, reduce duplicates

  treatment(des) <- rep(0:1, each = 5)

  expect_warning(units_of_assignment(des) <- data.frame(c1 = 1, c2 = c(1, 1:9)),
                 "collapsing")
  expect_equal(nrow(des@structure), 9)

  expect_error(units_of_assignment(des) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more units_of_assignment
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  des <- rd_design(z ~ unit_of_assignment(uoa1, uoa2) + block(bid) +
                     forcing(force), data = simdata)
  des2 <- des
  units_of_assignment(des2) <- df
  expect_identical(var_names(des2, "u"), colnames(df))
  expect_true(all(df == units_of_assignment(des2)))

  m <- matrix(1:40, ncol = 4)
  des2 <- des
  expect_error(units_of_assignment(des2) <- m,
               "must be named")
  colnames(m) <- letters[1:4]
  units_of_assignment(des2) <- m
  expect_identical(var_names(des2, "u"), colnames(m))
  expect_true(all(m == units_of_assignment(des2)))

  # uoa
  des <- rd_design(z ~ uoa(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  expect_identical(units_of_assignment(des), des@structure[, 2:3])

})

test_that("Accessing and replacing clusters", {
  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(clusters(des), des@structure[, 2:3])

  clusters(des) <- data.frame(uoa1 = 10:1, uoa2 = 1:10)
  expect_identical(var_names(des, "u"), c("uoa1", "uoa2"))
  expect_true(all(data.frame(uoa1 = 10:1, uoa2 = 1:10) == clusters(des)))

  clusters(des) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(var_names(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == clusters(des)))

  m <- matrix(1:20, ncol = 2)
  clusters(des) <- m
  expect_identical(var_names(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == clusters(des)))

  colnames(m) <- c("qwe", "asd")
  clusters(des) <- m
  expect_identical(var_names(des, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == clusters(des)))

  clusters(des)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == clusters(des)))

  # less clusters
  des2 <- des
  clusters(des2) <- m[, 1]
  expect_identical(var_names(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == clusters(des2)))

  des2 <- des
  clusters(des2) <- 10:1
  expect_identical(var_names(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == clusters(des2)))

  df <- data.frame(abc = 3:12)
  des2 <- des
  clusters(des2) <- df
  expect_identical(var_names(des2, "u"), colnames(df))
  expect_true(all(df == clusters(des2)))

  ########## Clusters, reduce duplicates

  treatment(des) <- rep(0:1, each = 5)

  expect_warning(clusters(des) <- data.frame(c1 = 1, c2 = c(1, 1:9)),
                 "collapsing")
  expect_equal(nrow(des@structure), 9)

  expect_error(clusters(des) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more clusters
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  des2 <- des
  clusters(des2) <- df
  expect_identical(var_names(des2, "u"), colnames(df))
  expect_true(all(df == clusters(des2)))

  m <- matrix(1:40, ncol = 4)
  des2 <- des
  expect_error(clusters(des2) <- m,
               "must be named")
  colnames(m) <- letters[1:4]
  clusters(des2) <- m
  expect_identical(var_names(des2, "u"), colnames(m))
  expect_true(all(m == clusters(des2)))
})

test_that("Accessing and replacing unitids", {
  data(simdata)

  des <- rd_design(z ~ unitid(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(unitids(des), des@structure[, 2:3])

  unitids(des) <- data.frame(uoa1 = 10:1, uoa2 = 1:10)
  expect_identical(var_names(des, "u"), c("uoa1", "uoa2"))
  expect_true(all(data.frame(uoa1 = 10:1, uoa2 = 1:10) == unitids(des)))

  unitids(des) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(var_names(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == unitids(des)))

  m <- matrix(1:20, ncol = 2)
  unitids(des) <- m
  expect_identical(var_names(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == unitids(des)))

  colnames(m) <- c("qwe", "asd")
  unitids(des) <- m
  expect_identical(var_names(des, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == unitids(des)))

  unitids(des)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == unitids(des)))

  # less unitids
  des2 <- des
  unitids(des2) <- m[, 1]
  expect_identical(var_names(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == unitids(des2)))

  des2 <- des
  unitids(des2) <- 10:1
  expect_identical(var_names(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == unitids(des2)))

  df <- data.frame(abc = 3:12)
  des2 <- des
  unitids(des2) <- df
  expect_identical(var_names(des2, "u"), colnames(df))
  expect_true(all(df == unitids(des2)))

  ########## Unitids, reduce duplicates

  treatment(des) <- rep(0:1, each = 5)

  expect_warning(unitids(des) <- data.frame(c1 = 1, c2 = c(1, 1:9)),
                 "collapsing")
  expect_equal(nrow(des@structure), 9)

  expect_error(unitids(des) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more unitids
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  des <- rd_design(z ~ unitid(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  des2 <- des
  unitids(des2) <- df
  expect_identical(var_names(des2, "u"), colnames(df))
  expect_true(all(df == unitids(des2)))

  m <- matrix(1:40, ncol = 4)
  des2 <- des
  expect_error(unitids(des2) <- m,
               "must be named")
  colnames(m) <- letters[1:4]
  unitids(des2) <- m
  expect_identical(var_names(des2, "u"), colnames(m))
  expect_true(all(m == unitids(des2)))

})

test_that("Accessing and replacing uoa vs unitid vs cluster", {
  data(simdata)

  des <- rd_design(z ~ unitid(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  # created with unitid
  expect_error(clusters(des), "not `cluster")
  expect_error(units_of_assignment(des), "not `unit_of_assignment")

  # created with cluster
  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  expect_error(unitids(des), "not `unitid")
  expect_error(units_of_assignment(des), "not `unit_of_assignment")

  # created with unit_of_assignment
  des <- rd_design(z ~ unit_of_assignment(uoa1, uoa2) + block(bid) +
                     forcing(force), data = simdata)
  expect_error(unitids(des), "not `unitid")
  expect_error(clusters(des), "not `cluster")

  # created with uoa
  des <- rd_design(z ~ uoa(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  expect_error(unitids(des), "not `unitid")
  expect_error(clusters(des), "not `cluster")
})

test_that("Accessing and replacing blocks", {
  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(blocks(des), des@structure[, 4, drop = FALSE])

  blocks(des) <- 1:10
  expect_equal(blocks(des), data.frame(bid = 1:10))

  mat <- matrix(10:1)
  blocks(des) <- mat
  expect_equal(blocks(des), data.frame(bid = 10:1))

  mat <- matrix(1:10)
  colnames(mat) <- "bid"
  blocks(des) <- mat
  expect_equal(blocks(des), data.frame(bid = 1:10))

  mat <- matrix(10:1)
  colnames(mat) <- "abc"
  blocks(des) <- mat
  expect_equal(blocks(des), data.frame(abc = 10:1))

  df <- data.frame(bid = 1:10)
  blocks(des) <- df
  expect_equal(blocks(des), df)

  blocks(des)[1:2, 1] <- 1
  expect_equal(blocks(des), data.frame(bid = c(1, 1, 3:10)))

  # More blocks
  df <- data.frame(bid1 = 1:10,
                   bid2 = 10:1)
  blocks(des) <- df
  expect_equal(blocks(des), df)

  expect_error(blocks(des) <- data.frame(a = 1:5),
               "same number")

  expect_error(blocks(des) <- 1:5,
               "same number")


  ### multi-dimensional blocks

  simdata2 <- simdata
  simdata2$uoa1 <- 1:50
  names(simdata2)[1:2] <- c("cid", "bida")
  des <- rd_design(z ~ cluster(cid) + block(bid, bida) + forcing(force),
                   data = simdata2)

  expect_identical(blocks(des), des@structure[, 3:4])

  blocks(des) <- data.frame(bid = 50:1, bida = 1:50)
  expect_identical(var_names(des, "b"), c("bid", "bida"))
  expect_true(all(data.frame(bid = 50:1, bida = 1:50) == blocks(des)))

  blocks(des) <- data.frame(abc1 = 50:1, abc2 = 1:50)
  expect_identical(var_names(des, "b"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 50:1, abc2 = 1:50) == blocks(des)))

  m <- matrix(1:100, ncol = 2)
  blocks(des) <- m
  expect_identical(var_names(des, "b"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:50, abc2 = 51:100) == blocks(des)))

  colnames(m) <- c("qwe", "asd")
  blocks(des) <- m
  expect_identical(var_names(des, "b"), c("qwe", "asd"))
  expect_true(all(data.frame(qwe = 1:50, asd = 51:100) == blocks(des)))

  blocks(des)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:50),
                             asd = c(100, 52:100)) == blocks(des)))

  # less blocks
  des2 <- des
  blocks(des2) <- 1:50
  expect_identical(var_names(des2, "b"), "qwe")
  expect_true(all(blocks(des2) == 1:50))

  des2 <- des
  blocks(des2) <- matrix(50:1, ncol = 1)
  expect_identical(var_names(des2, "b"), "qwe")
  expect_true(all(blocks(des2) == 50:1))

  des2 <- des
  df <- data.frame(bbb = 1:50)
  blocks(des2) <- df
  expect_identical(var_names(des2, "b"), "bbb")
  expect_true(all(blocks(des2) == df))

  # more blocks

  des2 <- des
  df <- data.frame(bbb = 1:50, ccc = 50:1, ddd = 1:50)
  blocks(des2) <- df
  expect_identical(var_names(des2, "b"), colnames(df))
  expect_true(all(blocks(des2) == df))
})


test_that("Accessing and replacing forcing", {
  data(simdata)

  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(forcings(des), des@structure[, 5, drop = FALSE])

  forcings(des) <- 1:10
  expect_equal(forcings(des), data.frame(force = 1:10))

  mat <- matrix(10:1)
  forcings(des) <- mat
  expect_equal(forcings(des), data.frame(force = 10:1))

  mat <- matrix(1:10)
  colnames(mat) <- "fvar"
  forcings(des) <- mat
  expect_equal(forcings(des), data.frame(fvar = 1:10))

  mat <- matrix(10:1)
  colnames(mat) <- "abc"
  forcings(des) <- mat
  expect_equal(forcings(des), data.frame(abc = 10:1))

  df <- data.frame(fvar = 1:10)
  forcings(des) <- df
  expect_equal(forcings(des), data.frame(fvar = 1:10))

  forcings(des)[1:2, 1] <- 1
  expect_equal(forcings(des), data.frame(fvar = c(1, 1, 3:10)))

  expect_error(forcings(des) <- data.frame(a = 1:5),
               "same number")

  # more forcing

  df <- data.frame(fvar = 1:10, fvar2 = 10:1)
  forcings(des) <- df
  expect_equal(var_names(des, "f"), colnames(df))
  expect_equal(forcings(des), df)

  ### multi-dimensional forcings

  simdata2 <- simdata
  simdata2$force2 <- rep(rnorm(10, mean = 5), times = rep(c(4, 6), times = 5))
  des <- rd_design(z ~ cluster(uoa1, uoa2) + block(bid) +
                     forcing(force, force2), data = simdata2)

  expect_identical(forcings(des), des@structure[, 5:6])

  forcings(des) <- data.frame(force = 10:1, force2 = 1:10)
  expect_identical(var_names(des, "f"), c("force", "force2"))
  expect_true(all(data.frame(force = 10:1, force2 = 1:10) == forcings(des)))

  forcings(des) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(var_names(des, "f"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == forcings(des)))

  m <- matrix(1:20, ncol = 2)
  forcings(des) <- m
  expect_identical(var_names(des, "f"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == forcings(des)))

  colnames(m) <- c("qwe", "asd")
  forcings(des) <- m
  expect_identical(var_names(des, "f"), c("qwe", "asd"))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == forcings(des)))

  forcings(des)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == forcings(des)))

  # less forcings
  des2 <- des
  forcings(des2) <- 1:10
  expect_identical(var_names(des2, "f"), "qwe")
  expect_true(all(forcings(des2) == 1:10))

  des2 <- des
  forcings(des2) <- matrix(10:1, ncol = 1)
  expect_identical(var_names(des2, "f"), "qwe")
  expect_true(all(forcings(des2) == 10:1))

  des2 <- des
  df <- data.frame(fff = 1:10)
  forcings(des2) <- df
  expect_identical(var_names(des2, "f"), "fff")
  expect_true(all(forcings(des2) == df))

  # more forcings

  des2 <- des
  df <- data.frame(fff = 1:10, ggg = 10:1, hhh = 1:10)
  forcings(des2) <- df
  expect_identical(var_names(des2, "f"), colnames(df))
  expect_true(all(forcings(des2) == df))


  # Forcing for non-RD

  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  expect_error(forcings(des),
               "only used")
  expect_error(forcings(des) <- rnorm(1:10),
               "only used")

  # duplicate variable names
  expect_error(rct_design(z ~ cluster(uoa1, uoa2) + block(uoa1),
                          data = simdata),
               "more than once")
  expect_error(rct_design(z ~ cluster(uoa1, uoa1, uoa2), data = simdata),
               "more than once")

})

test_that("Updating call works", {

  data(simdata)

  des <- rct_design(dose ~ unitid(uoa1, uoa2), data = simdata)

  expect_true(any(grepl("dose", des@call$formula[[2]])))
  treatment(des) <- data.frame(newt = 1:10)
  expect_false(all(grepl("dose", des@call$formula[[2]])))
  expect_true(any(grepl("newt", des@call$formula[[2]])))

  expect_true(any(grepl("uoa1", des@call$formula[[3]])))
  unitids(des) <- data.frame(newua = 1:10)
  expect_false(all(grepl("uoa1", des@call$formula[[3]])))
  expect_true(any(grepl("newua", des@call$formula[[3]])))

  des2 <- rd_design(dose ~ unitid(uoa1, uoa2) + block(bid) +
                      forcing(force), data = simdata)

  expect_true(any(grepl("bid", des2@call$formula[[3]])))
  blocks(des2) <- data.frame(newb = 1:10)
  expect_false(all(grepl("bid", des2@call$formula[[3]])))
  expect_true(any(grepl("newb", des2@call$formula[[3]])))

  expect_true(any(grepl("force", des2@call$formula[[3]])))
  forcings(des2) <- data.frame(newf = 1:10)
  expect_false(all(grepl("force", des2@call$formula[[3]])))
  expect_true(any(grepl("newf", des2@call$formula[[3]])))


})

test_that("treatment extraction with NA", {
  data(simdata)
  simdata$z[1:4] <- NA
  des <- rct_design(z ~ cluster(uoa1, uoa2), data = simdata)
  expect_identical(.bin_txt(des), des@structure$z)
  expect_identical(.bin_txt(des), treatment(des)[, 1])


  simdata$dose[1:4] <- NA
  des <- rct_design(dose ~ cluster(uoa1, uoa2), data = simdata)
  expect_error(.bin_txt(des), "Must provide a dichotomy")
  o <- treatment(des)[, 1]
  expect_true(is.na(o[1]))
  expect_true(all(!is.na(o[-1])))


  z <- .bin_txt(des, dichotomy = dose >= 250 ~ .)
  expect_true(is.na(z[1]))
  expect_true(all(z[-1] %in% 0:1))
  expect_identical(z, treatment(des, dichotomy = dose >= 250 ~ .)[, 1])
  expect_identical(o, treatment(des)[, 1])


})

test_that("convert_to_data.frame add'l testing", {
  data(simdata)
  des <- obs_design(z ~ cluster(uoa1, uoa2), data = simdata)

  ##### Replacing existing component
  # vector
  cdf <- .convert_to_data.frame(1:10, des, "t")
  expect_s3_class(cdf, "data.frame")
  expect_equal(dim(cdf), c(10, 1))
  expect_equal(names(cdf), "z")

  # unnamed matrix
  cdf <- .convert_to_data.frame(matrix(1:10, nrow = 10), des, "t")
  expect_equal(names(cdf), "z")

  # named matrix
  mm <- matrix(1:10, nrow = 10)
  colnames(mm) <- "w"
  cdf <- .convert_to_data.frame(mm, des, "t")
  expect_equal(names(cdf), "w")

  # unnamed data.frame
  df <- data.frame(1:10)
  colnames(df) <- NULL
  cdf <- .convert_to_data.frame(df, des, "t")
  expect_equal(names(cdf), "z")

  # named data.frame
  cdf <- .convert_to_data.frame(data.frame(a = 1:10), des, "t")
  expect_equal(names(cdf), "a")

  ##### Entering new components
  # vector
  expect_error(.convert_to_data.frame(1:10, des, "b"),
               "requires a named")

  # unnamed matrix
  expect_error(.convert_to_data.frame(matrix(1:10, nrow = 10), des, "b"),
               "requires a named")

  # named matrix
  mm <- matrix(1:10, nrow = 10)
  colnames(mm) <- "w"
  cdf <- .convert_to_data.frame(mm, des, "b")
  expect_equal(names(cdf), "w")

  # unnamed data.frame
  df <- data.frame(1:10)
  colnames(df) <- NULL
  expect_error(.convert_to_data.frame(df, des, "b"),
               "requires a named")

  # named data.frame
  cdf <- .convert_to_data.frame(data.frame(a = 1:10), des, "b")
  expect_equal(names(cdf), "a")

  ##### Multiple columns

  # Matrix
  mm <- matrix(1:20, nrow = 10)
  cdf <- .convert_to_data.frame(mm, des, "u")
  expect_s3_class(cdf, "data.frame")
  expect_equal(dim(cdf), c(10, 2))
  expect_equal(names(cdf), c("uoa1", "uoa2"))

  expect_error(.convert_to_data.frame(mm, des, "b"),
               "requires a named")

  colnames(mm) <- c("a", "b")
  cdf <- .convert_to_data.frame(mm, des, "u")
  expect_equal(names(cdf), c("a", "b"))
  cdf <- .convert_to_data.frame(mm, des, "b")
  expect_equal(names(cdf), c("a", "b"))

  # data.frame
  df <- data.frame(a = 1:10, b = 1:10)
  cdf <- .convert_to_data.frame(df, des, "u")
  expect_s3_class(cdf, "data.frame")
  expect_equal(dim(cdf), c(10, 2))
  expect_equal(names(cdf), c("a", "b"))
  cdf <- .convert_to_data.frame(df, des, "b")
  expect_equal(names(cdf), c("a", "b"))

  colnames(df) <- NULL
  cdf <- .convert_to_data.frame(df, des, "u")
  expect_equal(names(cdf), c("uoa1", "uoa2"))

  expect_error(.convert_to_data.frame(df, des, "b"),
               "requires a named")


  ##### Extra columns

  # Matrix
  mm <- matrix(1:30, nrow = 10)
  expect_error(.convert_to_data.frame(mm, des, "u"),
               "must be named")

  expect_error(.convert_to_data.frame(mm, des, "b"),
               "requires a named")

  colnames(mm) <- c("a", "b", "c")
  cdf <- .convert_to_data.frame(mm, des, "u")
  expect_s3_class(cdf, "data.frame")
  expect_equal(dim(cdf), c(10, 3))
  expect_equal(names(cdf), c("a", "b", "c"))
  cdf <- .convert_to_data.frame(mm, des, "b")
  expect_equal(names(cdf), c("a", "b", "c"))

  # data.frame
  df <- data.frame(a = 1:10, b = 1:10, c = 1:10)
  cdf <- .convert_to_data.frame(df, des, "u")
  expect_equal(names(cdf), c("a", "b", "c"))
  cdf <- .convert_to_data.frame(df, des, "b")
  expect_equal(names(cdf), c("a", "b", "c"))

  colnames(df) <- NULL
  expect_error(.convert_to_data.frame(df, des, "u"),
               "must be named")
  expect_error(.convert_to_data.frame(df, des, "b"),
               "requires a named")



  expect_error(.convert_to_data.frame(1:10, des, "q"),
               "Invalid type")

  expect_error(.convert_to_data.frame(lm(1~1), des, "u"))
})

test_that("has_blocks", {
  data(simdata)
  
  blocked_des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  no_blocks_des <- rct_design(z ~ cluster(uoa1, uoa2), simdata)
  
  expect_true(has_blocks(blocked_des))
  expect_false(has_blocks(no_blocks_des))
})
