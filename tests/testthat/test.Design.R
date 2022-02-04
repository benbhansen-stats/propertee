fc <- call("ls")

test_that("Design creation", {
  d <- new("Design",
           structure = data.frame(a = as.factor(c(1, 0)), b = c(2, 4)),
           columnIndex = c("t", "u"),
           type = "RCT",
           unitOfAssignmentType = "cluster",
           call = fc)

  expect_s4_class(d, "Design")
  expect_s3_class(d@structure, "data.frame")
  expect_type(d@columnIndex, "character")
  expect_type(d@type, "character")
  expect_type(d@unitOfAssignmentType, "character")
  expect_type(d@call, "language")
  expect_true(is(d@call, "call"))

  expect_equal(ncol(d@structure), length(d@columnIndex))
  expect_length(d@type, 1)
})

test_that("Design validity", {

  expect_error(new("Design",
                   structure = data.frame(),
                   columnIndex = "t",
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "positive dimensions")

  expect_error(new("Design",
                   structure = data.frame(a = as.factor(c(1, 0)),
                                          a = as.factor(c(1, 0))),
                   columnIndex = c("t", "t"),
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "one treatment")

  expect_error(new("Design",
                   structure = data.frame(a = c(1, 0), a = c(1, 0)),
                   columnIndex = c("u", "b"),
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "Missing treatment")

  expect_error(new("Design",
                   structure = data.frame(a = as.factor(c(1, 0)),
                                          b = c(2, 4)),
                   columnIndex = c("t", "b"),
                   type = "abc",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "unknown @type")

  expect_error(new("Design",
                   structure = data.frame(a = c(1, 2),
                                          b = as.factor(c(1, 0)),
                                          c = c(2, 0)),
                   columnIndex = c("u", "t"),
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "number of columns")

  expect_error(new("Design",
                   structure = data.frame(a = as.factor(c(1, 0)),
                                          b = c(2, 4)),
                   columnIndex = c("t", "k"),
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "unknown elements")

  expect_error(new("Design",
                   structure = data.frame(a = c(1, 0),
                                          b = c(1, 0)),
                   columnIndex = c("t", "u"),
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "must be factor")

  expect_error(new("Design",
                   structure = data.frame(a = as.factor(c(1, 0)),
                                          b = c(1, 0)),
                   columnIndex = c("t", "f"),
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "Forcing variables only valid")

  expect_error(new("Design",
                   structure = data.frame(a = as.factor(c(1, 0)),
                                          b = c(1, 0)),
                   columnIndex = c("t", "u"),
                   type = "RD",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "at least one forcing")

})

test_that("Design creation", {
  data(mtcars)
  mtcars <- mtcars[-c(5, 11),]

  d_rct <- New_Design(vs ~ cluster(qsec), data = mtcars, type = "RCT", call = fc)

  expect_s4_class(d_rct, "Design")
  expect_s3_class(d_rct@structure, "data.frame")
  expect_type(d_rct@columnIndex, "character")
  expect_type(d_rct@type, "character")
  expect_type(d_rct@unitOfAssignmentType, "character")
  expect_type(d_rct@call, "language")
  expect_true(is(d_rct@call, "call"))

  expect_equal(dim(d_rct@structure), c(30, 2))
  expect_length(d_rct@columnIndex, 2)
  expect_length(d_rct@type, 1)

  mtcars_subset <- subset(mtcars, select = c("vs", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(d_rct@structure == mtcars_subset))
  expect_equal(d_rct@columnIndex, c("t", "u"), ignore_attr = TRUE)
  expect_equal(d_rct@unitOfAssignmentType, "cluster")
  expect_equal(d_rct@type, "RCT")

  # subset

  d_obs <- New_Design(vs ~ cluster(qsec), data = mtcars, type = "Obs",
                  subset = mtcars$mpg > 17, call = fc)

  expect_equal(dim(d_obs@structure), c(sum(mtcars$mpg >  17), 2))
  expect_length(d_obs@columnIndex, 2)
  expect_length(d_obs@type, 1)

  mtcars_subset <- subset(mtcars, select = c("vs", "qsec"),
                          subset = mtcars$mpg > 17)
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(d_obs@structure == mtcars_subset))
  expect_equal(d_obs@columnIndex, c("t", "u"), ignore_attr = TRUE)
  expect_equal(d_obs@unitOfAssignmentType, "cluster")
  expect_equal(d_obs@type, "Obs")

  ### Complex design
  d_rd <- New_Design(vs ~ block(disp, gear) + forcing(wt, cyl) +
                       cluster(mpg, qsec),
                  data = mtcars, type = "RD", call = fc)

  mtcars_subset <- subset(mtcars, select = c("vs", "disp", "gear",
                                             "wt", "cyl", "mpg", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(d_rd@structure == mtcars_subset))
  expect_equal(d_rd@columnIndex, c("t", "b", "b", "f", "f", "u", "u"),
               ignore_attr = TRUE)
  expect_equal(d_rd@unitOfAssignmentType, "cluster")
  expect_equal(d_rd@type, "RD")

  ### Complex design with unitid
  d_rd2 <- New_Design(vs ~ block(disp, gear) + forcing(wt, cyl) +
                        unitid(mpg, qsec),
                  data = mtcars, type = "RD", call = fc)

  mtcars_subset <- subset(mtcars, select = c("vs", "disp", "gear",
                                             "wt", "cyl", "mpg", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(d_rd2@structure == mtcars_subset))
  expect_equal(d_rd2@columnIndex, c("t", "b", "b", "f", "f", "u", "u"),
               ignore_attr = TRUE)
  expect_equal(d_rd2@unitOfAssignmentType, "unitid")
  expect_equal(d_rd2@type, "RD")

  ### Specific designs
  rct_des <- RCT_Design(vs ~ cluster(qsec), data = mtcars)
  rct_des@call <- fc
  expect_identical(d_rct, rct_des)

  obs_des <- Obs_Design(vs ~ cluster(qsec), data = mtcars,
                        subset = mtcars$mpg > 17)
  obs_des@call <- fc
  expect_identical(d_obs, obs_des)

  rd_des <- RD_Design(vs ~ block(disp, gear) + forcing(wt, cyl) +
                        cluster(mpg, qsec),
                      data = mtcars)
  rd_des@call <- fc
  expect_identical(d_rd, rd_des)


  # Missing call

  expect_warning(des <- New_Design(vs ~ cluster(qsec), data = mtcars, type = "RCT"),
                 "Invalid call")
  expect_identical(des@call[[1]], as.name("New_Design"))

  expect_warning(des <- New_Design(vs ~ cluster(qsec), data = mtcars, type = "RCT", call = 1),
                 "Invalid call")
  expect_identical(des@call[[1]], as.name("New_Design"))
})

test_that("unit of assignment differs from unit of analysis", {

  data(simdata)

  desrct <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  expect_s4_class(desrct, "Design")
  expect_equal(nrow(desrct@structure), 10)

  expect_output(expect_error(RCT_Design(z ~ cluster(cid1) + block(bid),
                                        data = simdata),
                             "must be constant"),
                "cid1")

  data(mtcars)
  mtcars$prop <- rep(1:8, 4)

  expect_output(expect_error(RCT_Design(vs ~ cluster(prop), data = mtcars),
                             "must be constant"),
                "prop")
  expect_output(expect_error(RCT_Design(vs ~ cluster(prop), data = mtcars),
                             "must be constant"),
                "...")

})

test_that("Design printing", {
  data(simdata)

  desrct <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  desrd  <-  RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                       data = simdata)
  desobs <- Obs_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  expect_output(print(desrct), "Randomized")
  expect_output(show(desrct),  "Randomized")

  expect_silent(invisible(capture.output(expect_identical(desrct, show(desrct)))))
  expect_silent(invisible(capture.output(expect_identical(desrct, print(desrct)))))

  expect_output(print(desrd), "Discontinuity")
  expect_output(show(desrd),  "Discontinuity")
  expect_output(print(desobs), "Observational")
  expect_output(show(desobs),  "Observational")

  expect_output(show(desrct), "z")
  expect_output(show(desrct), "cid1")
  expect_output(show(desrct), "bid")

  expect_output(show(desobs), "z")
  expect_output(show(desobs), "cid1")
  expect_output(show(desobs), "bid")

  expect_output(show(desrd), "z")
  expect_output(show(desrd), "cid1")
  expect_output(show(desrd), "bid")
  expect_output(show(desrd), "force")

})

test_that("Accessing and replacing elements", {
  data(simdata)

  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)

  ##### Treatment

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


  ##### Clusters

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

  # unitid vs clusters

  expect_error(clusters(des), "not `cluster")
  des <- RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)
  expect_error(unitids(des), "not `unitid")

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


test_that("support for different types of treatment variables", {
  data(mtcars)
  mtcars <- mtcars[-c(5, 11),]

  # 0/1 binary
  treat_binary <- New_Design(vs ~ cluster(qsec), data = mtcars,
                             type = "RCT", call = fc)

  # logical
  mtcars$vs <- as.logical(mtcars$vs)
  treat_logical <- New_Design(vs ~ cluster(qsec), data = mtcars,
                              type = "RCT", call = fc)
  expect_identical(as.numeric(treat_binary@structure$vs),
                   as.numeric(treat_logical@structure$vs))


  # Factor
  mtcars$carb <- as.factor(mtcars$carb)
  treat_factor <- New_Design(carb ~ cluster(qsec), data = mtcars,
                             type = "RCT", call = fc)
  expect_identical(unique(treat_factor@structure$carb), unique(mtcars$carb))

  # Ordinal
  mtcars$gear <- as.ordered(mtcars$gear)
  treat_ordered <- New_Design(gear ~ cluster(qsec), data = mtcars,
                              type = "RCT", call = fc)
  expect_identical(unique(treat_ordered@structure$gear), unique(mtcars$gear))

  # non-binary numerical
  expect_error(New_Design(mpg ~ cluster(qsec), data = mtcars,
                          type = "RCT", call = fc),
               "only contain values 0 and 1")

  # transformation error
  expect_error(New_Design(as.factor(gear) ~ cluster(qsec), data = mtcars,
                          type = "RCT", call = fc),
               "variable transformations")
})

test_that("Design type conversions", {

  # Helper function to check equivalnce of Designs without worrying about some
  # oddities in @call
  checkequivDesign <- function(d1, d2) {
    expect_identical(d1@type, d2@type)
    expect_identical(d1@structure, d2@structure)
    expect_identical(d1@columnIndex, d2@columnIndex)
    expect_identical(d1@unitOfAssignmentType, d2@unitOfAssignmentType)

    # Due to formula hacking in the @call, the calls can be non-identical.
    # However, they should both produce the same Design if `eval`'d.
    d1redux <- eval(d1@call)
    d2redux <- eval(d2@call)
    expect_identical(d1redux@type, d2redux@type)
    expect_identical(d1redux@structure, d2redux@structure)
    expect_identical(d1redux@columnIndex, d2redux@columnIndex)
    expect_identical(d1redux@unitOfAssignmentType, d2redux@unitOfAssignmentType)

    invisible(TRUE)
  }
  data(simdata)

  desrct <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  desrd  <-  RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                       data = simdata)
  desobs <- Obs_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  # Obs <-> RCT
  expect_identical(desrct, as_RCT_Design(desobs))
  expect_identical(desobs, as_Obs_Design(desrct))

  # RD -> Obs/RCT
  expect_error(as_RCT_Design(desrd), "will be dropped")
  expect_error(as_Obs_Design(desrd), "will be dropped")
  checkequivDesign(desrct, as_RCT_Design(desrd, loseforcing = TRUE))
  checkequivDesign(desobs, as_Obs_Design(desrd, loseforcing = TRUE))

  # RCT/RD -> Obs
  expect_error(as_RD_Design(desrct))
  checkequivDesign(desrd, as_RD_Design(desrct, data = simdata,
                                      forcing = ~ . + forcing(force)))
  checkequivDesign(desrd, as_RD_Design(desrct, data = simdata,
                                      forcing = . ~ . + forcing(force)))
  checkequivDesign(desrd, as_RD_Design(desobs, data = simdata,
                                      forcing = ~ . + forcing(force)))
  checkequivDesign(desrd, as_RD_Design(desobs, data = simdata,
                                       forcing = . ~ . + forcing(force)))
  expect_error(as_RD_Design(desobs, data = simdata, forcing = simdata$force),
               "as a formula")

})
