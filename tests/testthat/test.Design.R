# fake call
fc <- call("ls")

identical_designs <- function(a, b) {
  attr(a@dichotomy, ".Environment") <- NULL
  attr(b@dichotomy, ".Environment") <- NULL
  expect_identical(a, b)
}

test_that("Design creation", {

  tests <- function(d) {
      expect_s4_class(d, "Design")
      expect_s3_class(d@structure, "data.frame")
      expect_type(d@column_index, "character")
      expect_type(d@type, "character")
      expect_type(d@unit_of_assignment_type, "character")
      expect_type(d@call, "language")
      expect_true(inherits(d@call, "call"))

      expect_equal(ncol(d@structure), length(d@column_index))
      expect_length(d@type, 1)
  }

  # binary treatment
  d <- new("Design",
           structure = data.frame(a = 0:1, b = c(2, 4)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

  # logical treatment
  d <- new("Design",
           structure = data.frame(a = c(TRUE, FALSE), b = c(2, 4)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

  # >2 treatment levels
  d <- new("Design",
           structure = data.frame(a = 0:2, b = c(2, 4, 6)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

  # character treatment

  d <- new("Design",
           structure = data.frame(a = letters[1:2], b = c(2, 4)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

  # >2 character treatment levels
  d <- new("Design",
           structure = data.frame(a = letters[1:3], b = c(2, 4, 6)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

})

test_that("Design validity", {

  expect_error(new("Design",
                   structure = data.frame(),
                   column_index = "t",
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "positive dimensions")

  expect_error(new("Design",
                   structure = data.frame(a = 0:1,
                                          b = 0:1),
                   column_index = c("t", "t"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "one treatment")

  expect_error(new("Design",
                   structure = data.frame(a = c(1, 0), a = c(1, 0)),
                   column_index = c("u", "b"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "Missing treatment")

  expect_error(new("Design",
                   structure = data.frame(a = c(0, 0), a = c(1, 0)),
                   column_index = c("t", "u"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "treatment can not be constant")

  expect_error(new("Design",
                   structure = data.frame(a = 0:1,
                                          b = c(2, 4)),
                   column_index = c("t", "b"),
                   type = "abc",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "unknown @type")

  expect_error(new("Design",
                   structure = data.frame(a = c(1, 2),
                                          b = as.factor(c(1, 0)),
                                          c = c(2, 0)),
                   column_index = c("t", "u", "u", "u"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "number of columns")

  expect_error(new("Design",
                   structure = data.frame(a = 0:1,
                                          b = c(2, 4)),
                   column_index = c("t", "k"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "unknown elements")

  expect_error(new("Design",
                   structure = data.frame(a = 0:1,
                                          b = c(1, 0)),
                   column_index = c("t", "f"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "Forcing variables only valid")

  expect_error(new("Design",
                   structure = data.frame(a = 0:1,
                                          b = c(1, 0)),
                   column_index = c("t", "u"),
                   type = "RD",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "at least one forcing")

})

test_that("Design creation", {
  data(mtcars)
  mtcars <- mtcars[-c(5, 11), ]

  d_rct <- .new_Design(vs ~ cluster(qsec), data = mtcars,
                      type = "RCT", call = fc, dichotomy = NULL)

  expect_s4_class(d_rct, "Design")
  expect_s3_class(d_rct@structure, "data.frame")
  expect_type(d_rct@column_index, "character")
  expect_type(d_rct@type, "character")
  expect_type(d_rct@unit_of_assignment_type, "character")
  expect_type(d_rct@call, "language")
  expect_true(inherits(d_rct@call, "call"))

  expect_equal(dim(d_rct@structure), c(30, 2))
  expect_length(d_rct@column_index, 2)
  expect_length(d_rct@type, 1)

  mtcars_subset <- subset(mtcars, select = c("vs", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(d_rct@structure == mtcars_subset))
  expect_equal(d_rct@column_index, c("t", "u"), ignore_attr = TRUE)
  expect_equal(d_rct@unit_of_assignment_type, "cluster")
  expect_equal(d_rct@type, "RCT")

  # subset

  d_obs <- .new_Design(vs ~ cluster(qsec), data = mtcars, type = "Obs",
                  subset = mtcars$mpg > 17, call = fc)

  expect_equal(dim(d_obs@structure), c(sum(mtcars$mpg >  17), 2))
  expect_length(d_obs@column_index, 2)
  expect_length(d_obs@type, 1)

  mtcars_subset <- subset(mtcars, select = c("vs", "qsec"),
                          subset = mtcars$mpg > 17)
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(d_obs@structure == mtcars_subset))
  expect_equal(d_obs@column_index, c("t", "u"), ignore_attr = TRUE)
  expect_equal(d_obs@unit_of_assignment_type, "cluster")
  expect_equal(d_obs@type, "Obs")

  ### Complex design
  d_rd <- .new_Design(vs ~ block(disp, gear) + forcing(wt, cyl) +
                       cluster(mpg, qsec),
                  data = mtcars, type = "RD", call = fc)

  mtcars_subset <- subset(mtcars, select = c("vs", "disp", "gear",
                                             "wt", "cyl", "mpg", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(d_rd@structure == mtcars_subset))
  expect_equal(d_rd@column_index, c("t", "b", "b", "f", "f", "u", "u"),
               ignore_attr = TRUE)
  expect_equal(d_rd@unit_of_assignment_type, "cluster")
  expect_equal(d_rd@type, "RD")

  ### Complex design with unitid
  d_rd2 <- .new_Design(vs ~ block(disp, gear) + forcing(wt, cyl) +
                        unitid(mpg, qsec),
                  data = mtcars, type = "RD", call = fc)

  mtcars_subset <- subset(mtcars, select = c("vs", "disp", "gear",
                                             "wt", "cyl", "mpg", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(d_rd2@structure == mtcars_subset))
  expect_equal(d_rd2@column_index, c("t", "b", "b", "f", "f", "u", "u"),
               ignore_attr = TRUE)
  expect_equal(d_rd2@unit_of_assignment_type, "unitid")
  expect_equal(d_rd2@type, "RD")

  ### Specific designs
  rct_des <- rct_design(vs ~ cluster(qsec), data = mtcars)
  rct_des@call <- fc
  identical_designs(d_rct, rct_des)

  obs_des <- obs_design(vs ~ cluster(qsec), data = mtcars,
                        subset = mtcars$mpg > 17)
  obs_des@call <- fc
  identical_designs(d_obs, obs_des)

  rd_des <- rd_design(vs ~ block(disp, gear) + forcing(wt, cyl) +
                        cluster(mpg, qsec),
                      data = mtcars)
  rd_des@call <- fc
  identical_designs(d_rd, rd_des)


  # Missing call

  expect_warning(des <- .new_Design(vs ~ cluster(qsec), data = mtcars,
                                   type = "RCT"),
                 "Invalid call")
  expect_identical(des@call[[1]], as.name(".new_Design"))

  expect_warning(des <- .new_Design(vs ~ cluster(qsec), data = mtcars,
                                   type = "RCT", call = 1),
                 "Invalid call")
  expect_identical(des@call[[1]], as.name(".new_Design"))
})

test_that("unit of assignment differs from unit of analysis", {

  data(simdata)

  desrct <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  expect_s4_class(desrct, "Design")
  expect_equal(nrow(desrct@structure), 10)

  expect_output(expect_error(rct_design(z ~ cluster(cid1) + block(bid),
                                        data = simdata),
                             "must be constant"),
                "cid1")

  data(mtcars)
  mtcars$prop <- rep(1:8, 4)

  expect_output(expect_error(rct_design(vs ~ cluster(prop), data = mtcars),
                             "must be constant"),
                "prop")
  expect_output(expect_error(rct_design(vs ~ cluster(prop), data = mtcars),
                             "must be constant"),
                "...")

  data(simdata)
  simdata$z[1] <- 1
  expect_output(expect_error(rct_design(z ~ cluster(cid1, cid2),
                                        data = simdata),
                             "must be constant"),
                "cid1")

  simdata$z[1] <- NA
  expect_output(expect_error(rct_design(z ~ cluster(cid1, cid2),
                                        data = simdata),
                             "must be constant"),
                "cid1")
})

test_that("Design printing", {
  data(simdata)

  desrct <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  desrd  <-  rd_design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                       data = simdata)
  desobs <- obs_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  expect_output(print(desrct), "Randomized")
  expect_output(show(desrct),  "Randomized")

  expect_silent(invisible(capture.output(expect_identical(desrct,
                                                          show(desrct)))))
  expect_silent(invisible(capture.output(expect_identical(desrct,
                                                          print(desrct)))))

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

  expect_no_match(capture.output(show(desrct)), "Dichotomy rule")

  d <- o > 2 ~ .
  desdichot <- rct_design(o ~ cluster(cid1, cid2), data = simdata,
                          dichotomy = d)

  expect_output(show(desdichot), "Dichotomy rule")
  expect_output(show(desdichot), deparse(d))
})


test_that("Design type conversions", {

  # Helper function to check equivalence of Designs without worrying about some
  # oddities in @call
  check_equiv_design <- function(d1, d2) {
    expect_identical(d1@type, d2@type)
    expect_identical(d1@structure, d2@structure)
    expect_identical(d1@column_index, d2@column_index)
    expect_identical(d1@unit_of_assignment_type, d2@unit_of_assignment_type)

    # Due to formula hacking in the @call, the calls can be non-identical.
    # However, they should both produce the same Design if `eval`'d.
    d1redux <- eval(d1@call)
    d2redux <- eval(d2@call)
    expect_identical(d1redux@type, d2redux@type)
    expect_identical(d1redux@structure, d2redux@structure)
    expect_identical(d1redux@column_index, d2redux@column_index)
    expect_identical(d1redux@unit_of_assignment_type,
                     d2redux@unit_of_assignment_type)

    invisible(TRUE)
  }

  data(simdata)

  desrct <- rct_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  desrd  <-  rd_design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                       data = simdata)
  desobs <- obs_design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)

  # Obs <-> RCT
  expect_identical(desrct, as_rct_design(desobs))
  expect_identical(desobs, as_obs_design(desrct))

  # RD -> Obs/RCT
  expect_error(as_rct_design(desrd), "will be dropped")
  expect_error(as_obs_design(desrd), "will be dropped")
  check_equiv_design(desrct, as_rct_design(desrd, loseforcing = TRUE))
  check_equiv_design(desobs, as_obs_design(desrd, loseforcing = TRUE))

  # RCT/RD -> Obs
  expect_error(as_rd_design(desrct))
  check_equiv_design(desrd, as_rd_design(desrct, data = simdata,
                                      forcing = ~ . + forcing(force)))
  check_equiv_design(desrd, as_rd_design(desrct, data = simdata,
                                      forcing = . ~ . + forcing(force)))
  check_equiv_design(desrd, as_rd_design(desobs, data = simdata,
                                      forcing = ~ . + forcing(force)))
  check_equiv_design(desrd, as_rd_design(desobs, data = simdata,
                                       forcing = . ~ . + forcing(force)))
  expect_error(as_rd_design(desobs, data = simdata, forcing = simdata$force),
               "as a formula")

})


test_that("variable transformations in Design", {

  data(simdata)

  des <- rct_design(dose + 3 ~ uoa(cid1, cid2), data = simdata)
  expect_true(all(treatment(des)[, 1] %in% (simdata$dose + 3)))

  expect_warning(des <- rd_design(force > 4 & bid < 2 ~ cluster(cid1, cid2) +
                                    forcing(force),
                                  data = simdata),
                 "conditional logic")
  expect_true(is.logical(treatment(des)[, 1]))
})

test_that("#26 cbind in identification variables", {
  data(simdata)
  des1 <- obs_design(o ~ cluster(cbind(cid1, cid2)), data = simdata)
  des2 <- obs_design(o ~ cluster(cid1, cid2), data = simdata)

  des1@call <- fc
  des2@call <- fc
  expect_identical(des1, des2)

})

test_that("variable_concordance", {

  data(simdata)

  # No issues
  des <- rd_design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)

  expect_silent(expect_true(design_data_concordance(des, simdata)))

  # inconsistent block
  sd <- simdata
  sd$bid[1] <- 2

  expect_warning(a <- design_data_concordance(des, sd),
                 "Inconsistencies.*`bid`")
  expect_false(a)

  # Inconsistent block & force
  sd$force[2] <- 2

  expect_warning(expect_warning(a <- design_data_concordance(des, sd),
                                "Inconsistencies.*`bid`"),
                  "Inconsistencies.*`force`")
  expect_false(a)

  # inconsistent treatment
  sd <- simdata
  sd$z[1] <- 1

  expect_warning(a <- design_data_concordance(des, sd),
                 "Inconsistencies.*`z`")

  # Missing variables
  sd <- simdata
  sd$z <- NULL
  expect_warning(a <- design_data_concordance(des, sd),
                 "`z` not found")
  expect_false(a)
  expect_silent(a <- design_data_concordance(des, sd,
                                             warn_on_nonexistence = FALSE))
  expect_true(a)

  sd <- simdata
  sd$bid <- NULL
  expect_warning(a <- design_data_concordance(des, sd),
                 "`bid` not found")
  expect_false(a)
  expect_silent(a <- design_data_concordance(des, sd,
                                             warn_on_nonexistence = FALSE))
  expect_true(a)

  sd <- simdata
  sd$force <- NULL
  expect_warning(a <- design_data_concordance(des, sd),
                 "`force` not found")
  expect_false(a)
  expect_silent(a <- design_data_concordance(des, sd,
                                             warn_on_nonexistence = FALSE))
  expect_true(a)

  sd$cid1 <- NULL
  expect_error(design_data_concordance(des, sd),
               "missing a `by=`")

  # No block
  des <- rd_design(z ~ cluster(cid1, cid2) + forcing(force),
                   data = simdata)
  expect_silent(expect_true(design_data_concordance(des, simdata)))

  # multiple variables
  sd <- simdata
  sd$force2 <- sd$force + 1
  des <- rd_design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force, force2),
                   data = sd)
  expect_silent(expect_true(design_data_concordance(des, sd)))

  expect_warning(design_data_concordance(des, simdata),
                 "`force2` not found")

  sd$force[1] <- 1
  expect_warning(design_data_concordance(des, sd),
                 "Inconsistencies.*`force`$")

  sd2 <- simdata
  sd2$force[1] <- 1
  expect_warning(expect_warning(design_data_concordance(des, sd2),
                                "Inconsistencies.*`force`$"),
                 "`force2` not found")

  sd$force2[1] <- 2
  expect_warning(expect_warning(design_data_concordance(des, sd),
                                "Inconsistencies.*`force`$"),
                 "Inconsistencies.*`force2`$")

  data(simdata)

  des <- rd_design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force),
                   data = simdata)

  names(simdata)[names(simdata) == "cid1"] <- "ccc"
  names(simdata)[names(simdata) == "force"] <- "fff"

  expect_silent(expect_true({
    design_data_concordance(des, simdata,
                            by = list(force = "fff",
                                      cid1 = "ccc"))
  }))

})

test_that("NA's in data to create Design", {
  data(simdata)
  sd2 <- simdata
  sd2$cid1[1] <- NA

  expect_error(obs_design(z ~ unitid(cid1, cid2), data = sd2),
               "Missing values cannot be found")

  des1 <- obs_design(z ~ unitid(cid1, cid2), data = simdata)
  des2 <- obs_design(z ~ unitid(cid1, cid2), data = sd2, na.fail = FALSE)

  des1@call <- call("c")
  des2@call <- call("c")

  expect_identical(des1, des2)

  sd2 <- simdata
  sd2$bid <- NA

  expect_silent(des <- obs_design(z ~ unitid(cid1, cid2), data = sd2))
  expect_error(obs_design(z ~ unitid(cid1, cid2) + block(bid), data = sd2),
               "Missing values cannot be found")

  # Shouldn't affect treatment
  sd2 <- simdata
  sd2$z[1:10] <- NA
  expect_silent(des <- obs_design(z ~ unitid(cid1, cid2), data = sd2,
                                  na.fail = TRUE))
  expect_equal(nrow(des@structure), 10)
  expect_silent(des <- obs_design(z ~ unitid(cid1, cid2), data = sd2,
                                  na.fail = FALSE))
  expect_equal(nrow(des@structure), 10)

  sd2 <- simdata
  sd2$cid1[1:4] <- NA
  des <- obs_design(z ~ unitid(cid1, cid2), data = sd2, na.fail = FALSE)
  expect_equal(nrow(des@structure), 9)

})
