# fake call
fc <- call("ls")

test_that("Design creation", {

  tests <- function(d) {
      expect_s4_class(d, "Design")
      expect_s3_class(d@structure, "data.frame")
      expect_type(d@columnIndex, "character")
      expect_type(d@type, "character")
      expect_type(d@unitOfAssignmentType, "character")
      expect_type(d@call, "language")
      expect_true(is(d@call, "call"))

      expect_equal(ncol(d@structure), length(d@columnIndex))
      expect_length(d@type, 1)
  }

  # binary treatment
  d <- new("Design",
           structure = data.frame(a = factor(0:1), b = c(2, 4)),
           columnIndex = c("t", "u"),
           type = "RCT",
           unitOfAssignmentType = "cluster",
           call = fc)

  tests(d)

  # >2 treatment levels
  d <- new("Design",
           structure = data.frame(a = as.factor(c(1, 0, 2)), b = c(2, 4, 6)),
           columnIndex = c("t", "u"),
           type = "RCT",
           unitOfAssignmentType = "cluster",
           call = fc)

  tests(d)

  # ordinal treatment

  d <- new("Design",
           structure = data.frame(a = as.ordered(c(1, 0, 2)), b = c(2, 4, 6)),
           columnIndex = c("t", "u"),
           type = "RCT",
           unitOfAssignmentType = "cluster",
           call = fc)

  tests(d)
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
                   structure = data.frame(a = as.factor(c(0, 0)), a = c(1, 0)),
                   columnIndex = c("t", "u"),
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "treatment can not be constant")

  expect_error(new("Design",
                   structure = data.frame(a = factor(c(0, 0), levels = c(0,1)), a = c(1, 0)),
                   columnIndex = c("t", "u"),
                   type = "RCT",
                   unitOfAssignmentType = "cluster",
                   call = fc),
               "treatment can not be constant")

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

test_that("different treatment types", {
  data(simdata)

  # z currently numeric
  d1 <- New_Design(z ~ cluster(cid1, cid2), data = simdata, type = "RCT", call = fc)
  # factor
  simdata$z <- as.factor(simdata$z)
  d2 <- New_Design(z ~ cluster(cid1, cid2), data = simdata, type = "RCT", call = fc)
  expect_identical(d1, d2)
  # binary
  simdata$z <- simdata$z == 1
  d3 <- New_Design(z ~ cluster(cid1, cid2), data = simdata, type = "RCT", call = fc)
  expect_equal(cor(as.numeric(d3@structure$z), as.numeric(d1@structure$z)), 1)
  # different levels, but same underlying values

  # o currently ordinal
  e1 <- New_Design(o ~ cluster(cid1, cid2), data = simdata, type = "RCT", call = fc)
  expect_true(is(e1@structure[, e1@columnIndex == "t"], "ordered"))
  # factor
  simdata$o <- as.factor(as.numeric(simdata$o))
  e2 <- New_Design(o ~ cluster(cid1, cid2), data = simdata, type = "RCT", call = fc)
  expect_false(is(e2@structure[, e2@columnIndex == "t"], "ordered"))
})


test_that("Empty levels in treatment are removed", {
  data(simdata)

  simdata$z <- factor(simdata$z, levels = 0:2)
  expect_warning(desrct <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata),
                 "Empty levels")
  expect_length(levels(desrct@structure[, desrct@columnIndex == "t"]), 2)

  simdata$o <- factor(simdata$o, levels = 1:7)
  expect_warning(desrct <- RCT_Design(o ~ cluster(cid1, cid2) + block(bid), data = simdata),
                 "Empty levels")
  expect_length(levels(desrct@structure[, desrct@columnIndex == "t"]), 4)

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


test_that("variable transformations in Design", {

  data(simdata)

  des <- RCT_Design(as.factor(dose) ~ uoa(cid1, cid2), data = simdata)
  expect_s3_class(des@structure[,des@columnIndex == "t"], "factor")
  expect_equal(length(levels(des@structure[,des@columnIndex == "t"])),
               length(unique(simdata$dose)))

  des <- RD_Design(force > 4 & bid < 2 ~ cluster(cid1, cid2) +  forcing(force), data = simdata)
  expect_s3_class(des@structure[,des@columnIndex == "t"], "factor")
  expect_length(levels(des@structure[,des@columnIndex == "t"]), 2)

})

test_that("#26 cbind in identification variables", {
  data(simdata)
  des1 <- Obs_Design(o ~ cluster(cbind(cid1, cid2)), data = simdata)
  des2 <- Obs_Design(o ~ cluster(cid1, cid2), data = simdata)

  des1@call <- fc
  des2@call <- fc
  expect_identical(des1, des2)

})
