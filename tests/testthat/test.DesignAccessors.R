test_that("Accessing and replacing elements", {
  data(simdata)

  ##### Treatment
  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(treatment(des), des@structure[, 1, drop = FALSE])

  treatment(des) <- as.factor(rep(0:1, each = 5))
  expect_equal(treatment(des), data.frame(z = as.factor(rep(0:1, each = 5))))

  df <- data.frame(z = as.factor(rep(0:1, each = 5)))
  treatment(des) <- df
  expect_equal(treatment(des), data.frame(z = as.factor(rep(0:1, each = 5))))

  treatment(des)[1:2,1] <- 1
  expect_equal(treatment(des),
               data.frame(z = as.factor(rep(c(1,0,1),
                                            times = c(2,3,5)))))

  expect_error(treatment(des) <- as.factor(1:5),
               "same number")

  expect_error(treatment(des) <- data.frame(a = as.factor(c(1,0,1,0,1))),
               "same number")

  expect_error(treatment(des) <- matrix(1:5, ncol = 1),
               "Treatment must be")


  ##### UnitsOfAssignment
  des <- RD_Design(z ~ unitOfAssignment(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(unitsOfAssignment(des), des@structure[, 2:3])

  unitsOfAssignment(des) <- data.frame(cid1 = 10:1, cid2 = 1:10)
  expect_identical(varNames(des, "u"), c("cid1", "cid2"))
  expect_true(all(data.frame(cid1 = 10:1, cid2 = 1:10) == unitsOfAssignment(des)))

  unitsOfAssignment(des) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(varNames(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == unitsOfAssignment(des)))

  m <- matrix(1:20, ncol = 2)
  unitsOfAssignment(des) <- m
  expect_identical(varNames(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == unitsOfAssignment(des)))

  colnames(m) <- c("qwe", "asd")
  unitsOfAssignment(des) <- m
  expect_identical(varNames(des, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == unitsOfAssignment(des)))

  unitsOfAssignment(des)[1,1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == unitsOfAssignment(des)))

  # less unitsOfAssignment
  des2 <- des
  unitsOfAssignment(des2) <- m[,1]
  expect_identical(varNames(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == unitsOfAssignment(des2)))

  des2 <- des
  unitsOfAssignment(des2) <- 10:1
  expect_identical(varNames(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == unitsOfAssignment(des2)))

  df <- data.frame(abc = 3:12)
  des2 <- des
  unitsOfAssignment(des2) <- df
  expect_identical(varNames(des2, "u"), colnames(df))
  expect_true(all(df == unitsOfAssignment(des2)))

  ########## UnitsOfAssignment, reduce duplicates

  treatment(des) <- as.factor(rep(0:1, each = 5))

  expect_warning(unitsOfAssignment(des) <- data.frame(c1 = 1, c2 = c(1,1:9)),
                 "collapsing")
  expect_equal(nrow(des@structure), 9)

  expect_error(unitsOfAssignment(des) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more unitsOfAssignment
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  des <- RD_Design(z ~ unitOfAssignment(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)
  des2 <- des
  unitsOfAssignment(des2) <- df
  expect_identical(varNames(des2, "u"), colnames(df))
  expect_true(all(df == unitsOfAssignment(des2)))

  m <- matrix(1:40, ncol = 4)
  des2 <- des
  expect_error(unitsOfAssignment(des2) <- m,
               "be named")
  colnames(m) <- letters[1:4]
  unitsOfAssignment(des2) <- m
  expect_identical(varNames(des2, "u"), colnames(m))
  expect_true(all(m == unitsOfAssignment(des2)))

  # uoa
  des <- RD_Design(z ~ uoa(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)
  expect_identical(unitsOfAssignment(des), des@structure[, 2:3])


  ##### Clusters
  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(clusters(des), des@structure[, 2:3])

  clusters(des) <- data.frame(cid1 = 10:1, cid2 = 1:10)
  expect_identical(varNames(des, "u"), c("cid1", "cid2"))
  expect_true(all(data.frame(cid1 = 10:1, cid2 = 1:10) == clusters(des)))

  clusters(des) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(varNames(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == clusters(des)))

  m <- matrix(1:20, ncol = 2)
  clusters(des) <- m
  expect_identical(varNames(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == clusters(des)))

  colnames(m) <- c("qwe", "asd")
  clusters(des) <- m
  expect_identical(varNames(des, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == clusters(des)))

  clusters(des)[1,1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == clusters(des)))

  # less clusters
  des2 <- des
  clusters(des2) <- m[,1]
  expect_identical(varNames(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == clusters(des2)))

  des2 <- des
  clusters(des2) <- 10:1
  expect_identical(varNames(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == clusters(des2)))

  df <- data.frame(abc = 3:12)
  des2 <- des
  clusters(des2) <- df
  expect_identical(varNames(des2, "u"), colnames(df))
  expect_true(all(df == clusters(des2)))

  ########## Clusters, reduce duplicates

  treatment(des) <- as.factor(rep(0:1, each = 5))

  expect_warning(clusters(des) <- data.frame(c1 = 1, c2 = c(1,1:9)),
                 "collapsing")
  expect_equal(nrow(des@structure), 9)

  expect_error(clusters(des) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more clusters
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)
  des2 <- des
  clusters(des2) <- df
  expect_identical(varNames(des2, "u"), colnames(df))
  expect_true(all(df == clusters(des2)))

  m <- matrix(1:40, ncol = 4)
  des2 <- des
  expect_error(clusters(des2) <- m,
               "be named")
  colnames(m) <- letters[1:4]
  clusters(des2) <- m
  expect_identical(varNames(des2, "u"), colnames(m))
  expect_true(all(m == clusters(des2)))

  ##### Unitids

  des <- RD_Design(z ~ unitid(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)

  expect_identical(unitids(des), des@structure[, 2:3])

  unitids(des) <- data.frame(cid1 = 10:1, cid2 = 1:10)
  expect_identical(varNames(des, "u"), c("cid1", "cid2"))
  expect_true(all(data.frame(cid1 = 10:1, cid2 = 1:10) == unitids(des)))

  unitids(des) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(varNames(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == unitids(des)))

  m <- matrix(1:20, ncol = 2)
  unitids(des) <- m
  expect_identical(varNames(des, "u"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == unitids(des)))

  colnames(m) <- c("qwe", "asd")
  unitids(des) <- m
  expect_identical(varNames(des, "u"), colnames(m))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == unitids(des)))

  unitids(des)[1,1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == unitids(des)))

  # less unitids
  des2 <- des
  unitids(des2) <- m[,1]
  expect_identical(varNames(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 1:10) == unitids(des2)))

  des2 <- des
  unitids(des2) <- 10:1
  expect_identical(varNames(des2, "u"), "qwe")
  expect_true(all(data.frame(qwe = 10:1) == unitids(des2)))

  df <- data.frame(abc = 3:12)
  des2 <- des
  unitids(des2) <- df
  expect_identical(varNames(des2, "u"), colnames(df))
  expect_true(all(df == unitids(des2)))

  ########## Unitids, reduce duplicates

  treatment(des) <- as.factor(rep(0:1, each = 5))

  expect_warning(unitids(des) <- data.frame(c1 = 1, c2 = c(1,1:9)),
                 "collapsing")
  expect_equal(nrow(des@structure), 9)

  expect_error(unitids(des) <- data.frame(c1 = 1, c2 = c(1:8, 1)),
               "non-constant")

  # more unitids
  df <- data.frame(abc = 3:12, def = 4:13, efg = 5:14)
  des <- RD_Design(z ~ unitid(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)
  des2 <- des
  unitids(des2) <- df
  expect_identical(varNames(des2, "u"), colnames(df))
  expect_true(all(df == unitids(des2)))

  m <- matrix(1:40, ncol = 4)
  des2 <- des
  expect_error(unitids(des2) <- m,
               "be named")
  colnames(m) <- letters[1:4]
  unitids(des2) <- m
  expect_identical(varNames(des2, "u"), colnames(m))
  expect_true(all(m == unitids(des2)))

  ##### unitOfAssignment vs unitid vs clusters

  # created with unitid
  expect_error(clusters(des), "not `cluster")
  expect_error(unitsOfAssignment(des), "not `unitOfAssignment")

  # created with cluster
  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)
  expect_error(unitids(des), "not `unitid")
  expect_error(unitsOfAssignment(des), "not `unitOfAssignment")

  # created with unitOfAssignment
  des <- RD_Design(z ~ unitOfAssignment(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)
  expect_error(unitids(des), "not `unitid")
  expect_error(clusters(des), "not `cluster")

  # created with uoa
  des <- RD_Design(z ~ uoa(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)
  expect_error(unitids(des), "not `unitid")
  expect_error(clusters(des), "not `cluster")

  ##### Blocks

  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
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

  blocks(des)[1:2,1] <- 1
  expect_equal(blocks(des), data.frame(bid = c(1,1,3:10)))

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
  simdata2$cid1 <- 1:50
  names(simdata2)[1:2] <- c("cid", "bida")
  des <- RD_Design(z ~ cluster(cid) + block(bid, bida) + forcing(force),
                   data = simdata2)

  expect_identical(blocks(des), des@structure[, 3:4])

  blocks(des) <- data.frame(bid = 50:1, bida = 1:50)
  expect_identical(varNames(des, "b"), c("bid", "bida"))
  expect_true(all(data.frame(bid = 50:1, bida = 1:50) == blocks(des)))

  blocks(des) <- data.frame(abc1 = 50:1, abc2 = 1:50)
  expect_identical(varNames(des, "b"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 50:1, abc2 = 1:50) == blocks(des)))

  m <- matrix(1:100, ncol = 2)
  blocks(des) <- m
  expect_identical(varNames(des, "b"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:50, abc2 = 51:100) == blocks(des)))

  colnames(m) <- c("qwe", "asd")
  blocks(des) <- m
  expect_identical(varNames(des, "b"), c("qwe", "asd"))
  expect_true(all(data.frame(qwe = 1:50, asd = 51:100) == blocks(des)))

  blocks(des)[1,1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:50),
                             asd = c(100, 52:100)) == blocks(des)))

  # less blocks
  des2 <- des
  blocks(des2) <- 1:50
  expect_identical(varNames(des2, "b"), "qwe")
  expect_true(all(blocks(des2) == 1:50))

  des2 <- des
  blocks(des2) <- matrix(50:1, ncol = 1)
  expect_identical(varNames(des2, "b"), "qwe")
  expect_true(all(blocks(des2) == 50:1))

  des2 <- des
  df <- data.frame(bbb = 1:50)
  blocks(des2) <- df
  expect_identical(varNames(des2, "b"), "bbb")
  expect_true(all(blocks(des2) == df))

  # more blocks

  des2 <- des
  df <- data.frame(bbb = 1:50, ccc = 50:1, ddd = 1:50)
  blocks(des2) <- df
  expect_identical(varNames(des2, "b"), colnames(df))
  expect_true(all(blocks(des2) == df))


  #### Forcing
  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
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

  forcings(des)[1:2,1] <- 1
  expect_equal(forcings(des), data.frame(fvar = c(1,1,3:10)))

  expect_error(forcings(des) <- data.frame(a = 1:5),
               "same number")

  # more forcing

  df <- data.frame(fvar = 1:10, fvar2 = 10:1)
  forcings(des) <- df
  expect_equal(varNames(des, "f"), colnames(df))
  expect_equal(forcings(des), df)

  ### multi-dimensional forcings

  simdata2 <- simdata
  simdata2$force2 <- rep(rnorm(10, mean = 5), times = rep(c(4, 6), times = 5))
  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force, force2),
                   data = simdata2)

  expect_identical(forcings(des), des@structure[, 5:6])

  forcings(des) <- data.frame(force = 10:1, force2 = 1:10)
  expect_identical(varNames(des, "f"), c("force", "force2"))
  expect_true(all(data.frame(force = 10:1, force2 = 1:10) == forcings(des)))

  forcings(des) <- data.frame(abc1 = 10:1, abc2 = 1:10)
  expect_identical(varNames(des, "f"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 10:1, abc2 = 1:10) == forcings(des)))

  m <- matrix(1:20, ncol = 2)
  forcings(des) <- m
  expect_identical(varNames(des, "f"), c("abc1", "abc2"))
  expect_true(all(data.frame(abc1 = 1:10, abc2 = 11:20) == forcings(des)))

  colnames(m) <- c("qwe", "asd")
  forcings(des) <- m
  expect_identical(varNames(des, "f"), c("qwe", "asd"))
  expect_true(all(data.frame(qwe = 1:10, asd = 11:20) == forcings(des)))

  forcings(des)[1,1:2] <- 100
  expect_true(all(data.frame(qwe = c(100, 2:10),
                             asd = c(100, 12:20)) == forcings(des)))

  # less forcings
  des2 <- des
  forcings(des2) <- 1:10
  expect_identical(varNames(des2, "f"), "qwe")
  expect_true(all(forcings(des2) == 1:10))

  des2 <- des
  forcings(des2) <- matrix(10:1, ncol = 1)
  expect_identical(varNames(des2, "f"), "qwe")
  expect_true(all(forcings(des2) == 10:1))

  des2 <- des
  df <- data.frame(fff = 1:10)
  forcings(des2) <- df
  expect_identical(varNames(des2, "f"), "fff")
  expect_true(all(forcings(des2) == df))

  # more forcings

  des2 <- des
  df <- data.frame(fff = 1:10, ggg = 10:1, hhh = 1:10)
  forcings(des2) <- df
  expect_identical(varNames(des2, "f"), colnames(df))
  expect_true(all(forcings(des2) == df))


  # Forcing for non-RD

  des <- RCT_Design(z ~ cluster(cid1, cid2), data = simdata)
  expect_error(forcings(des),
               "only used")
  expect_error(forcings(des) <- rnorm(1:10),
               "only used")

  # duplicate variable names
  expect_error(RCT_Design(z ~ cluster(cid1, cid2) + block(cid1), data = simdata),
               "more than once")
  expect_error(RCT_Design(z ~ cluster(cid1, cid1, cid2), data = simdata),
               "more than once")

})
