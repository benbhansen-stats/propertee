test_that("Accessing and replacing treatment", {
  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_equal(dim(treatment(spec)), c(10, 1))
  expect_true(is.numeric(treatment(spec)[, 1]))
  expect_identical(treatment(spec), spec@structure[, 1, drop = FALSE])

  treatment(spec) <- rep(0:1, each = 5)
  expect_equal(treatment(spec), data.frame(z = rep(0:1, each = 5)))

  df <- data.frame(z = rep(0:1, each = 5))
  treatment(spec) <- df
  expect_equal(treatment(spec), data.frame(z = rep(0:1, each = 5)))

  treatment(spec)[1:2, 1] <- 1
  expect_equal(treatment(spec),
               data.frame(z = rep(c(1, 0, 1),
                                  times = c(2, 3, 5))))

  expect_error(treatment(spec) <- 1:5,
               "same number")

  expect_error(treatment(spec) <- data.frame(a = c(1, 0, 1, 0, 1)),
               "same number")

  # Continuous treatment, no dichotomy
  spec <- rd_spec(dose ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(treatment(spec), spec@structure[, 1, drop = FALSE])

  # Continuous treatment, dichotomy
  spec <- rd_spec(dose ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(treatment(spec), spec@structure[, 1, drop = FALSE])
  tt <- treatment(spec, dichotomy = dose <= 100 ~ dose == 200)

  expect_equal(dim(tt), c(10, 1))
  expect_true(is.numeric(tt[, 1]))
  expect_true(all(tt[, 1] %in% c(0, 1, NA)))
  expect_equal(colnames(tt), "z__")

  ############ .bin_txt

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(treatment(spec)[, 1], .bin_txt(spec))

  spec <- rd_spec(o ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_true(all(.bin_txt(spec, dichotomy = o >= 3 ~ .) %in% c(0, 1, NA)))
  expect_identical(.bin_txt(spec, dichotomy = o >= 3 ~ .),
                   treatment(spec, dichotomy = o >= 3 ~ .)[, 1])

  spec <- rd_spec(o ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_error(.bin_txt(spec), "Must provide a dichotomy")

})

test_that("New .bin_txt, .expand_txt, and .apply_dichotomy errors", {
  spec <- rct_spec(z ~ unitid(uoa1, uoa2), simdata)
  expect_error(.bin_txt(spec, data.frame("uoa1" = seq(100, 110))),
               "Not all unit of assignment variables")

  expect_error(.apply_dichotomy(simdata$dose, dose > 50 ~ dose == 50),
               "expected to be a named `data.frame`")
  expect_error(.apply_dichotomy(treatment(spec), dichotomy = dose > 50 ~ dose == 50),
               "Could not find variables specified in dichotomy")
  expect_error(.apply_dichotomy(simdata, . ~ .),
               "At least one side")
  expect_error(.apply_dichotomy(simdata, dichotomy = dose >= 50 ~ dose == 50),
               "dichotomy overlaps")
})

test_that("Accessing and replacing unit of assignment", {
  data(simdata)

  spec <- rd_spec(z ~ unit_of_assignment(uoa1, uoa2) + block(bid) +
                     forcing(force),
                   data = simdata)

  expect_identical(units_of_assignment(spec), spec@structure[, 2:3])

  units_of_assignment(spec) <- data.frame(uoa1 = 10:1, uoa2 = 1:10)
  expect_identical(var_names(spec, "u"), c("uoa1", "uoa2"))
  expect_true(all(data.frame(uoa1 = 10:1, uoa2 = 1:10) ==
                    units_of_assignment(spec)))

  units_of_assignment(spec) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(var_names(spec, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) ==
                    units_of_assignment(spec)))

  m <- matrix(1:20, ncol = 2)
  units_of_assignment(spec) <- m
  expect_identical(var_names(spec, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) ==
                    units_of_assignment(spec)))

  colnames(m) <- c("qwe", "asd")
  units_of_assignment(spec) <- m
  expect_identical(var_names(spec, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) ==
                    units_of_assignment(spec)))

  units_of_assignment(spec)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) ==
                    units_of_assignment(spec)))

  # less units_of_assignment
  spec2 <- spec
  units_of_assignment(spec2) <- m[, 1]
  expect_identical(var_names(spec2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == units_of_assignment(spec2)))

  spec2 <- spec
  units_of_assignment(spec2) <- 10:1
  expect_identical(var_names(spec2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == units_of_assignment(spec2)))

  df <- data.frame(abc = 3:12)
  spec2 <- spec
  units_of_assignment(spec2) <- df
  expect_identical(var_names(spec2, "u"), colnames(df))
  expect_true(all(df == units_of_assignment(spec2)))

  ########## UnitsOfAssignment, reduce duplicates

  treatment(spec) <- rep(0:1, each = 5)

  expect_warning(units_of_assignment(spec) <- data.frame(c1 = 1, c2 = c(1, 1:9)),
                 "collapsing")
  expect_equal(nrow(spec@structure), 9)

  expect_error(units_of_assignment(spec) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more units_of_assignment
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  spec <- rd_spec(z ~ unit_of_assignment(uoa1, uoa2) + block(bid) +
                     forcing(force), data = simdata)
  spec2 <- spec
  units_of_assignment(spec2) <- df
  expect_identical(var_names(spec2, "u"), colnames(df))
  expect_true(all(df == units_of_assignment(spec2)))

  m <- matrix(1:40, ncol = 4)
  spec2 <- spec
  expect_error(units_of_assignment(spec2) <- m,
               "must be named")
  colnames(m) <- letters[1:4]
  units_of_assignment(spec2) <- m
  expect_identical(var_names(spec2, "u"), colnames(m))
  expect_true(all(m == units_of_assignment(spec2)))

  # uoa
  spec <- rd_spec(z ~ uoa(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  expect_identical(units_of_assignment(spec), spec@structure[, 2:3])

})

test_that("Accessing and replacing clusters", {
  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(clusters(spec), spec@structure[, 2:3])

  clusters(spec) <- data.frame(uoa1 = 10:1, uoa2 = 1:10)
  expect_identical(var_names(spec, "u"), c("uoa1", "uoa2"))
  expect_true(all(data.frame(uoa1 = 10:1, uoa2 = 1:10) == clusters(spec)))

  clusters(spec) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(var_names(spec, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == clusters(spec)))

  m <- matrix(1:20, ncol = 2)
  clusters(spec) <- m
  expect_identical(var_names(spec, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == clusters(spec)))

  colnames(m) <- c("qwe", "asd")
  clusters(spec) <- m
  expect_identical(var_names(spec, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == clusters(spec)))

  clusters(spec)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == clusters(spec)))

  # less clusters
  spec2 <- spec
  clusters(spec2) <- m[, 1]
  expect_identical(var_names(spec2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == clusters(spec2)))

  spec2 <- spec
  clusters(spec2) <- 10:1
  expect_identical(var_names(spec2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == clusters(spec2)))

  df <- data.frame(abc = 3:12)
  spec2 <- spec
  clusters(spec2) <- df
  expect_identical(var_names(spec2, "u"), colnames(df))
  expect_true(all(df == clusters(spec2)))

  ########## Clusters, reduce duplicates

  treatment(spec) <- rep(0:1, each = 5)

  expect_warning(clusters(spec) <- data.frame(c1 = 1, c2 = c(1, 1:9)),
                 "collapsing")
  expect_equal(nrow(spec@structure), 9)

  expect_error(clusters(spec) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more clusters
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  spec2 <- spec
  clusters(spec2) <- df
  expect_identical(var_names(spec2, "u"), colnames(df))
  expect_true(all(df == clusters(spec2)))

  m <- matrix(1:40, ncol = 4)
  spec2 <- spec
  expect_error(clusters(spec2) <- m,
               "must be named")
  colnames(m) <- letters[1:4]
  clusters(spec2) <- m
  expect_identical(var_names(spec2, "u"), colnames(m))
  expect_true(all(m == clusters(spec2)))
})

test_that("Accessing and replacing unitids", {
  data(simdata)

  spec <- rd_spec(z ~ unitid(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(unitids(spec), spec@structure[, 2:3])

  unitids(spec) <- data.frame(uoa1 = 10:1, uoa2 = 1:10)
  expect_identical(var_names(spec, "u"), c("uoa1", "uoa2"))
  expect_true(all(data.frame(uoa1 = 10:1, uoa2 = 1:10) == unitids(spec)))

  unitids(spec) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(var_names(spec, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == unitids(spec)))

  m <- matrix(1:20, ncol = 2)
  unitids(spec) <- m
  expect_identical(var_names(spec, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == unitids(spec)))

  colnames(m) <- c("qwe", "asd")
  unitids(spec) <- m
  expect_identical(var_names(spec, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == unitids(spec)))

  unitids(spec)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == unitids(spec)))

  # less unitids
  spec2 <- spec
  unitids(spec2) <- m[, 1]
  expect_identical(var_names(spec2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == unitids(spec2)))

  spec2 <- spec
  unitids(spec2) <- 10:1
  expect_identical(var_names(spec2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == unitids(spec2)))

  df <- data.frame(abc = 3:12)
  spec2 <- spec
  unitids(spec2) <- df
  expect_identical(var_names(spec2, "u"), colnames(df))
  expect_true(all(df == unitids(spec2)))

  ########## Unitids, reduce duplicates

  treatment(spec) <- rep(0:1, each = 5)

  expect_warning(unitids(spec) <- data.frame(c1 = 1, c2 = c(1, 1:9)),
                 "collapsing")
  expect_equal(nrow(spec@structure), 9)

  expect_error(unitids(spec) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more unitids
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  spec <- rd_spec(z ~ unitid(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  spec2 <- spec
  unitids(spec2) <- df
  expect_identical(var_names(spec2, "u"), colnames(df))
  expect_true(all(df == unitids(spec2)))

  m <- matrix(1:40, ncol = 4)
  spec2 <- spec
  expect_error(unitids(spec2) <- m,
               "must be named")
  colnames(m) <- letters[1:4]
  unitids(spec2) <- m
  expect_identical(var_names(spec2, "u"), colnames(m))
  expect_true(all(m == unitids(spec2)))

})

test_that("Accessing and replacing uoa vs unitid vs cluster", {
  data(simdata)

  spec <- rd_spec(z ~ unitid(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  # created with unitid
  expect_error(clusters(spec), "not `cluster")
  expect_error(units_of_assignment(spec), "not `unit_of_assignment")

  # created with cluster
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  expect_error(unitids(spec), "not `unitid")
  expect_error(units_of_assignment(spec), "not `unit_of_assignment")

  # created with unit_of_assignment
  spec <- rd_spec(z ~ unit_of_assignment(uoa1, uoa2) + block(bid) +
                     forcing(force), data = simdata)
  expect_error(unitids(spec), "not `unitid")
  expect_error(clusters(spec), "not `cluster")

  # created with uoa
  spec <- rd_spec(z ~ uoa(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)
  expect_error(unitids(spec), "not `unitid")
  expect_error(clusters(spec), "not `cluster")
})

test_that("Accessing and replacing blocks", {
  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(blocks(spec), spec@structure[, 4, drop = FALSE])

  blocks(spec) <- 1:10
  expect_equal(blocks(spec), data.frame(bid = 1:10))

  mat <- matrix(10:1)
  blocks(spec) <- mat
  expect_equal(blocks(spec), data.frame(bid = 10:1))

  mat <- matrix(1:10)
  colnames(mat) <- "bid"
  blocks(spec) <- mat
  expect_equal(blocks(spec), data.frame(bid = 1:10))

  mat <- matrix(10:1)
  colnames(mat) <- "abc"
  blocks(spec) <- mat
  expect_equal(blocks(spec), data.frame(abc = 10:1))

  df <- data.frame(bid = 1:10)
  blocks(spec) <- df
  expect_equal(blocks(spec), df)

  blocks(spec)[1:2, 1] <- 1
  expect_equal(blocks(spec), data.frame(bid = c(1, 1, 3:10)))

  # More blocks
  df <- data.frame(bid1 = 1:10,
                   bid2 = 10:1)
  blocks(spec) <- df
  expect_equal(blocks(spec), df)

  expect_error(blocks(spec) <- data.frame(a = 1:5),
               "same number")

  expect_error(blocks(spec) <- 1:5,
               "same number")


  ### multi-dimensional blocks

  simdata2 <- simdata
  simdata2$uoa1 <- 1:50
  names(simdata2)[1:2] <- c("cid", "bida")
  spec <- rd_spec(z ~ cluster(cid) + block(bid, bida) + forcing(force),
                   data = simdata2)

  expect_identical(blocks(spec), spec@structure[, 3:4])

  blocks(spec) <- data.frame(bid = 50:1, bida = 1:50)
  expect_identical(var_names(spec, "b"), c("bid", "bida"))
  expect_true(all(data.frame(bid = 50:1, bida = 1:50) == blocks(spec)))

  blocks(spec) <- data.frame(abc1 = 50:1, abc2 = 1:50)
  expect_identical(var_names(spec, "b"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 50:1, abc2 = 1:50) == blocks(spec)))

  m <- matrix(1:100, ncol = 2)
  blocks(spec) <- m
  expect_identical(var_names(spec, "b"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:50, abc2 = 51:100) == blocks(spec)))

  colnames(m) <- c("qwe", "asd")
  blocks(spec) <- m
  expect_identical(var_names(spec, "b"), c("qwe", "asd"))
  expect_true(all(data.frame(qwe = 1:50, asd = 51:100) == blocks(spec)))

  blocks(spec)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:50),
                             asd = c(100, 52:100)) == blocks(spec)))

  # less blocks
  spec2 <- spec
  blocks(spec2) <- 1:50
  expect_identical(var_names(spec2, "b"), "qwe")
  expect_true(all(blocks(spec2) == 1:50))

  spec2 <- spec
  blocks(spec2) <- matrix(50:1, ncol = 1)
  expect_identical(var_names(spec2, "b"), "qwe")
  expect_true(all(blocks(spec2) == 50:1))

  spec2 <- spec
  df <- data.frame(bbb = 1:50)
  blocks(spec2) <- df
  expect_identical(var_names(spec2, "b"), "bbb")
  expect_true(all(blocks(spec2) == df))

  # more blocks

  spec2 <- spec
  df <- data.frame(bbb = 1:50, ccc = 50:1, ddd = 1:50)
  blocks(spec2) <- df
  expect_identical(var_names(spec2, "b"), colnames(df))
  expect_true(all(blocks(spec2) == df))
})


test_that("Accessing and replacing forcing", {
  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(forcings(spec), spec@structure[, 5, drop = FALSE])

  forcings(spec) <- 1:10
  expect_equal(forcings(spec), data.frame(force = 1:10))

  mat <- matrix(10:1)
  forcings(spec) <- mat
  expect_equal(forcings(spec), data.frame(force = 10:1))

  mat <- matrix(1:10)
  colnames(mat) <- "fvar"
  forcings(spec) <- mat
  expect_equal(forcings(spec), data.frame(fvar = 1:10))

  mat <- matrix(10:1)
  colnames(mat) <- "abc"
  forcings(spec) <- mat
  expect_equal(forcings(spec), data.frame(abc = 10:1))

  df <- data.frame(fvar = 1:10)
  forcings(spec) <- df
  expect_equal(forcings(spec), data.frame(fvar = 1:10))

  forcings(spec)[1:2, 1] <- 1
  expect_equal(forcings(spec), data.frame(fvar = c(1, 1, 3:10)))

  expect_error(forcings(spec) <- data.frame(a = 1:5),
               "same number")

  # more forcing

  df <- data.frame(fvar = 1:10, fvar2 = 10:1)
  forcings(spec) <- df
  expect_equal(var_names(spec, "f"), colnames(df))
  expect_equal(forcings(spec), df)

  ### multi-dimensional forcings

  simdata2 <- simdata
  simdata2$force2 <- rep(rnorm(10, mean = 5), times = rep(c(4, 6), times = 5))
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) +
                     forcing(force, force2), data = simdata2)

  expect_identical(forcings(spec), spec@structure[, 5:6])

  forcings(spec) <- data.frame(force = 10:1, force2 = 1:10)
  expect_identical(var_names(spec, "f"), c("force", "force2"))
  expect_true(all(data.frame(force = 10:1, force2 = 1:10) == forcings(spec)))

  forcings(spec) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(var_names(spec, "f"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == forcings(spec)))

  m <- matrix(1:20, ncol = 2)
  forcings(spec) <- m
  expect_identical(var_names(spec, "f"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == forcings(spec)))

  colnames(m) <- c("qwe", "asd")
  forcings(spec) <- m
  expect_identical(var_names(spec, "f"), c("qwe", "asd"))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == forcings(spec)))

  forcings(spec)[1, 1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == forcings(spec)))

  # less forcings
  spec2 <- spec
  forcings(spec2) <- 1:10
  expect_identical(var_names(spec2, "f"), "qwe")
  expect_true(all(forcings(spec2) == 1:10))

  spec2 <- spec
  forcings(spec2) <- matrix(10:1, ncol = 1)
  expect_identical(var_names(spec2, "f"), "qwe")
  expect_true(all(forcings(spec2) == 10:1))

  spec2 <- spec
  df <- data.frame(fff = 1:10)
  forcings(spec2) <- df
  expect_identical(var_names(spec2, "f"), "fff")
  expect_true(all(forcings(spec2) == df))

  # more forcings

  spec2 <- spec
  df <- data.frame(fff = 1:10, ggg = 10:1, hhh = 1:10)
  forcings(spec2) <- df
  expect_identical(var_names(spec2, "f"), colnames(df))
  expect_true(all(forcings(spec2) == df))


  # Forcing for non-RD

  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  expect_error(forcings(spec),
               "only used")
  expect_error(forcings(spec) <- rnorm(1:10),
               "only used")

  # duplicate variable names
  expect_error(rct_spec(z ~ cluster(uoa1, uoa2) + block(uoa1),
                          data = simdata),
               "more than once")
  expect_error(rct_spec(z ~ cluster(uoa1, uoa1, uoa2), data = simdata),
               "more than once")

})

test_that("Updating call works", {

  data(simdata)

  spec <- rct_spec(dose ~ unitid(uoa1, uoa2), data = simdata)

  expect_true(any(grepl("dose", spec@call$formula[[2]])))
  treatment(spec) <- data.frame(newt = 1:10)
  expect_false(all(grepl("dose", spec@call$formula[[2]])))
  expect_true(any(grepl("newt", spec@call$formula[[2]])))

  expect_true(any(grepl("uoa1", spec@call$formula[[3]])))
  unitids(spec) <- data.frame(newua = 1:10)
  expect_false(all(grepl("uoa1", spec@call$formula[[3]])))
  expect_true(any(grepl("newua", spec@call$formula[[3]])))

  spec2 <- rd_spec(dose ~ unitid(uoa1, uoa2) + block(bid) +
                      forcing(force), data = simdata)

  expect_true(any(grepl("bid", spec2@call$formula[[3]])))
  blocks(spec2) <- data.frame(newb = 1:10)
  expect_false(all(grepl("bid", spec2@call$formula[[3]])))
  expect_true(any(grepl("newb", spec2@call$formula[[3]])))

  expect_true(any(grepl("force", spec2@call$formula[[3]])))
  forcings(spec2) <- data.frame(newf = 1:10)
  expect_false(all(grepl("force", spec2@call$formula[[3]])))
  expect_true(any(grepl("newf", spec2@call$formula[[3]])))


})

test_that("treatment extraction with NA", {
  data(simdata)
  simdata$z[1:4] <- NA
  spec <- rct_spec(z ~ cluster(uoa1, uoa2), data = simdata)
  expect_identical(.bin_txt(spec), spec@structure$z)
  expect_identical(.bin_txt(spec), treatment(spec)[, 1])


  simdata$dose[1:4] <- NA
  spec <- rct_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  expect_error(.bin_txt(spec), "Must provide a dichotomy")
  o <- treatment(spec)[, 1]
  expect_true(is.na(o[1]))
  expect_true(all(!is.na(o[-1])))


  z <- .bin_txt(spec, dichotomy = dose >= 250 ~ .)
  expect_true(is.na(z[1]))
  expect_true(all(z[-1] %in% 0:1))
  expect_identical(z, treatment(spec, dichotomy = dose >= 250 ~ .)[, 1])
  expect_identical(o, treatment(spec)[, 1])


})

test_that("convert_to_data.frame add'l testing", {
  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  ##### Replacing existing component
  # vector
  cdf <- .convert_to_data.frame(1:10, spec, "t")
  expect_s3_class(cdf, "data.frame")
  expect_equal(dim(cdf), c(10, 1))
  expect_equal(names(cdf), "z")

  # unnamed matrix
  cdf <- .convert_to_data.frame(matrix(1:10, nrow = 10), spec, "t")
  expect_equal(names(cdf), "z")

  # named matrix
  mm <- matrix(1:10, nrow = 10)
  colnames(mm) <- "w"
  cdf <- .convert_to_data.frame(mm, spec, "t")
  expect_equal(names(cdf), "w")

  # unnamed data.frame
  df <- data.frame(1:10)
  colnames(df) <- NULL
  cdf <- .convert_to_data.frame(df, spec, "t")
  expect_equal(names(cdf), "z")

  # named data.frame
  cdf <- .convert_to_data.frame(data.frame(a = 1:10), spec, "t")
  expect_equal(names(cdf), "a")

  ##### Entering new components
  # vector
  expect_error(.convert_to_data.frame(1:10, spec, "b"),
               "requires a named")

  # unnamed matrix
  expect_error(.convert_to_data.frame(matrix(1:10, nrow = 10), spec, "b"),
               "requires a named")

  # named matrix
  mm <- matrix(1:10, nrow = 10)
  colnames(mm) <- "w"
  cdf <- .convert_to_data.frame(mm, spec, "b")
  expect_equal(names(cdf), "w")

  # unnamed data.frame
  df <- data.frame(1:10)
  colnames(df) <- NULL
  expect_error(.convert_to_data.frame(df, spec, "b"),
               "requires a named")

  # named data.frame
  cdf <- .convert_to_data.frame(data.frame(a = 1:10), spec, "b")
  expect_equal(names(cdf), "a")

  ##### Multiple columns

  # Matrix
  mm <- matrix(1:20, nrow = 10)
  cdf <- .convert_to_data.frame(mm, spec, "u")
  expect_s3_class(cdf, "data.frame")
  expect_equal(dim(cdf), c(10, 2))
  expect_equal(names(cdf), c("uoa1", "uoa2"))

  expect_error(.convert_to_data.frame(mm, spec, "b"),
               "requires a named")

  colnames(mm) <- c("a", "b")
  cdf <- .convert_to_data.frame(mm, spec, "u")
  expect_equal(names(cdf), c("a", "b"))
  cdf <- .convert_to_data.frame(mm, spec, "b")
  expect_equal(names(cdf), c("a", "b"))

  # data.frame
  df <- data.frame(a = 1:10, b = 1:10)
  cdf <- .convert_to_data.frame(df, spec, "u")
  expect_s3_class(cdf, "data.frame")
  expect_equal(dim(cdf), c(10, 2))
  expect_equal(names(cdf), c("a", "b"))
  cdf <- .convert_to_data.frame(df, spec, "b")
  expect_equal(names(cdf), c("a", "b"))

  colnames(df) <- NULL
  cdf <- .convert_to_data.frame(df, spec, "u")
  expect_equal(names(cdf), c("uoa1", "uoa2"))

  expect_error(.convert_to_data.frame(df, spec, "b"),
               "requires a named")


  ##### Extra columns

  # Matrix
  mm <- matrix(1:30, nrow = 10)
  expect_error(.convert_to_data.frame(mm, spec, "u"),
               "must be named")

  expect_error(.convert_to_data.frame(mm, spec, "b"),
               "requires a named")

  colnames(mm) <- c("a", "b", "c")
  cdf <- .convert_to_data.frame(mm, spec, "u")
  expect_s3_class(cdf, "data.frame")
  expect_equal(dim(cdf), c(10, 3))
  expect_equal(names(cdf), c("a", "b", "c"))
  cdf <- .convert_to_data.frame(mm, spec, "b")
  expect_equal(names(cdf), c("a", "b", "c"))

  # data.frame
  df <- data.frame(a = 1:10, b = 1:10, c = 1:10)
  cdf <- .convert_to_data.frame(df, spec, "u")
  expect_equal(names(cdf), c("a", "b", "c"))
  cdf <- .convert_to_data.frame(df, spec, "b")
  expect_equal(names(cdf), c("a", "b", "c"))

  colnames(df) <- NULL
  expect_error(.convert_to_data.frame(df, spec, "u"),
               "must be named")
  expect_error(.convert_to_data.frame(df, spec, "b"),
               "requires a named")



  expect_error(.convert_to_data.frame(1:10, spec, "q"),
               "Invalid type")

  expect_error(.convert_to_data.frame(lm(1~1), spec, "u"))
})

test_that("has_blocks", {
  data(simdata)

  blocked_spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  no_blocks_spec <- rct_spec(z ~ cluster(uoa1, uoa2), simdata)

  expect_true(has_blocks(blocked_spec))
  expect_false(has_blocks(no_blocks_spec))
})
