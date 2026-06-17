# fake call
fc <- call("ls")

identical_specs <- function(a, b) {
  expect_identical(a, b)
}

test_that("StudySpecification creation", {

  tests <- function(d) {
      expect_s4_class(d, "StudySpecification")
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
  d <- new("StudySpecification",
           structure = data.frame(a = 0:1, b = c(2, 4)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

  # logical treatment
  d <- new("StudySpecification",
           structure = data.frame(a = c(TRUE, FALSE), b = c(2, 4)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

  # >2 treatment levels
  d <- new("StudySpecification",
           structure = data.frame(a = 0:2, b = c(2, 4, 6)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

  # character treatment

  d <- new("StudySpecification",
           structure = data.frame(a = letters[1:2], b = c(2, 4)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

  # >2 character treatment levels
  d <- new("StudySpecification",
           structure = data.frame(a = letters[1:3], b = c(2, 4, 6)),
           column_index = c("t", "u"),
           type = "RCT",
           unit_of_assignment_type = "cluster",
           call = fc)

  tests(d)

})

test_that("StudySpecification validity", {

  expect_error(new("StudySpecification",
                   structure = data.frame(),
                   column_index = "t",
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "positive dimensions")

  expect_error(new("StudySpecification",
                   structure = data.frame(a = 0:1,
                                          b = 0:1),
                   column_index = c("t", "t"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "one treatment")

  expect_error(new("StudySpecification",
                   structure = data.frame(a = c(1, 0), a = c(1, 0)),
                   column_index = c("u", "b"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "Missing treatment")

  expect_error(new("StudySpecification",
                   structure = data.frame(a = c(0, 0), a = c(1, 0)),
                   column_index = c("t", "u"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "treatment can not be constant")

  expect_error(new("StudySpecification",
                   structure = data.frame(a = 0:1,
                                          b = c(2, 4)),
                   column_index = c("t", "b"),
                   type = "abc",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "unknown @type")

  expect_error(new("StudySpecification",
                   structure = data.frame(a = c(1, 2),
                                          b = as.factor(c(1, 0)),
                                          c = c(2, 0)),
                   column_index = c("t", "u", "u", "u"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "number of columns")

  expect_error(new("StudySpecification",
                   structure = data.frame(a = 0:1,
                                          b = c(2, 4)),
                   column_index = c("t", "k"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "unknown elements")

  expect_error(new("StudySpecification",
                   structure = data.frame(a = 0:1,
                                          b = c(1, 0)),
                   column_index = c("t", "f"),
                   type = "RCT",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "Forcing variables only valid")

  expect_error(new("StudySpecification",
                   structure = data.frame(a = 0:1,
                                          b = c(1, 0)),
                   column_index = c("t", "u"),
                   type = "RD",
                   unit_of_assignment_type = "cluster",
                   call = fc),
               "at least one forcing")

})

test_that("StudySpecification creation", {
  data(mtcars)
  mtcars <- mtcars[-c(5, 11), ]

  spec_rct <- .new_StudySpecification(vs ~ cluster(qsec), data = mtcars,
                      type = "RCT",
                      call = quote(rct_spec(formula = vs ~ cluster(qsec), data = mtcars)))

  expect_s4_class(spec_rct, "StudySpecification")
  expect_s3_class(spec_rct@structure, "data.frame")
  expect_type(spec_rct@column_index, "character")
  expect_type(spec_rct@type, "character")
  expect_type(spec_rct@unit_of_assignment_type, "character")
  expect_type(spec_rct@call, "language")
  expect_true(inherits(spec_rct@call, "call"))

  expect_equal(dim(spec_rct@structure), c(30, 2))
  expect_length(spec_rct@column_index, 2)
  expect_length(spec_rct@type, 1)

  mtcars_subset <- subset(mtcars, select = c("vs", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(spec_rct@structure == mtcars_subset))
  expect_equal(spec_rct@column_index, c("t", "u"), ignore_attr = TRUE)
  expect_equal(spec_rct@unit_of_assignment_type, "cluster")
  expect_equal(spec_rct@type, "RCT")

  # subset

  spec_obs <- .new_StudySpecification(
    vs ~ cluster(qsec), data = mtcars, type = "Obs", subset = mtcars$mpg > 17,
    call = quote(obs_spec(formula = vs ~ cluster(qsec), data = mtcars, subset = mtcars$mpg > 17))
  )

  expect_equal(dim(spec_obs@structure), c(sum(mtcars$mpg >  17), 2))
  expect_length(spec_obs@column_index, 2)
  expect_length(spec_obs@type, 1)

  mtcars_subset <- subset(mtcars, select = c("vs", "qsec"),
                          subset = mtcars$mpg > 17)
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(spec_obs@structure == mtcars_subset))
  expect_equal(spec_obs@column_index, c("t", "u"), ignore_attr = TRUE)
  expect_equal(spec_obs@unit_of_assignment_type, "cluster")
  expect_equal(spec_obs@type, "Obs")

  ### Complex specification
  spec_rd <- .new_StudySpecification(
    vs ~ block(disp, gear) + forcing(wt, cyl) + cluster(mpg, qsec),
    data = mtcars, type = "RD",
    call = quote(rd_spec(formula = vs ~ block(disp, gear) + forcing(wt, cyl) + cluster(mpg, qsec), data = mtcars))
  )

  mtcars_subset <- subset(mtcars, select = c("vs", "disp", "gear",
                                             "wt", "cyl", "mpg", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(spec_rd@structure == mtcars_subset))
  expect_equal(spec_rd@column_index, c("t", "b", "b", "f", "f", "u", "u"),
               ignore_attr = TRUE)
  expect_equal(spec_rd@unit_of_assignment_type, "cluster")
  expect_equal(spec_rd@type, "RD")

  ### Complex specification with unitid
  spec_rd2 <- .new_StudySpecification(vs ~ block(disp, gear) + forcing(wt, cyl) +
                        unitid(mpg, qsec),
                  data = mtcars, type = "RD", call = fc)

  mtcars_subset <- subset(mtcars, select = c("vs", "disp", "gear",
                                             "wt", "cyl", "mpg", "qsec"))
  mtcars_subset$vs <- as.factor(mtcars_subset$vs)
  expect_true(all(spec_rd2@structure == mtcars_subset))
  expect_equal(spec_rd2@column_index, c("t", "b", "b", "f", "f", "u", "u"),
               ignore_attr = TRUE)
  expect_equal(spec_rd2@unit_of_assignment_type, "unitid")
  expect_equal(spec_rd2@type, "RD")

  ### Specific specifications
  rct_spec <- rct_spec(vs ~ cluster(qsec), data = mtcars)
  identical_specs(spec_rct, rct_spec)

  obs_spec <- obs_spec(vs ~ cluster(qsec), data = mtcars,
                        subset = mtcars$mpg > 17)
  identical_specs(spec_obs, obs_spec)

  rd_spec <- rd_spec(vs ~ block(disp, gear) + forcing(wt, cyl) +
                        cluster(mpg, qsec),
                      data = mtcars)
  identical_specs(spec_rd, rd_spec)


  # Missing call

  expect_warning(spec <- .new_StudySpecification(vs ~ cluster(qsec),
                                                 data = mtcars,
                                                 type = "RCT"),
                 "Invalid call")
  expect_identical(spec@call[[1]], as.name(".new_StudySpecification"))

  expect_warning(spec <- .new_StudySpecification(vs ~ cluster(qsec),
                                                 data = mtcars,
                                                 type = "RCT", call = 1),
                 "Invalid call")
  expect_identical(spec@call[[1]], as.name(".new_StudySpecification"))
})

test_that("unit of assignment differs from unit of analysis", {

  data(simdata)

  specrct <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  expect_s4_class(specrct, "StudySpecification")
  expect_equal(nrow(specrct@structure), 10)

  expect_error(rct_spec(z ~ cluster(uoa1) + block(bid), data = simdata),
               "uoa1")

  data(mtcars)
  mtcars$prop <- rep(1:8, 4)

  expect_error(rct_spec(vs ~ cluster(prop), data = mtcars),
               "prop")
  expect_error(rct_spec(vs ~ cluster(prop), data = mtcars),
               "...")

  data(simdata)
  simdata$z[1] <- 1
  expect_error(rct_spec(z ~ cluster(uoa1, uoa2), data = simdata),
               "uoa1")

  simdata$z[1] <- NA
  expect_error(rct_spec(z ~ cluster(uoa1, uoa2), data = simdata),
               "uoa1")
})

test_that("StudySpecification printing", {
  data(simdata)

  specrct <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  specrd  <-  rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                       data = simdata)
  specobs <- obs_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  expect_output(print(specrct), "Randomized")
  expect_output(show(specrct),  "Randomized")

  expect_silent(invisible(capture.output(expect_identical(specrct,
                                                          show(specrct)))))
  expect_silent(invisible(capture.output(expect_identical(specrct,
                                                          print(specrct)))))

  expect_output(print(specrd), "Discontinuity")
  expect_output(show(specrd),  "Discontinuity")
  expect_output(print(specobs), "Observational")
  expect_output(show(specobs),  "Observational")

  expect_output(show(specrct), "z")
  expect_output(show(specrct), "uoa1")
  expect_output(show(specrct), "bid")

  expect_output(show(specobs), "z")
  expect_output(show(specobs), "uoa1")
  expect_output(show(specobs), "bid")

  expect_output(show(specrd), "z")
  expect_output(show(specrd), "uoa1")
  expect_output(show(specrd), "bid")
  expect_output(show(specrd), "force")
})


test_that("StudySpecification type conversions", {

  # Helper function to check equivalence of StudySpecifications without worrying about some
  # oddities in @call
  check_equiv_spec <- function(d1, d2) {
    expect_identical(d1@type, d2@type)
    expect_identical(d1@structure, d2@structure)
    expect_identical(d1@column_index, d2@column_index)
    expect_identical(d1@unit_of_assignment_type, d2@unit_of_assignment_type)

    # Due to formula hacking in the @call, the calls can be non-identical.
    # However, they should both produce the same StudySpecification if `eval`'d.
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

  specrct <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)
  specrd  <-  rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                       data = simdata)
  specobs <- obs_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  # Obs <-> RCT
  expect_identical(specrct, as_rct_spec(specobs))
  expect_identical(specobs, as_obs_spec(specrct))

  # RD -> Obs/RCT
  expect_error(as_rct_spec(specrd), "will be dropped")
  expect_error(as_obs_spec(specrd), "will be dropped")
  check_equiv_spec(specrct, as_rct_spec(specrd, loseforcing = TRUE))
  check_equiv_spec(specobs, as_obs_spec(specrd, loseforcing = TRUE))

  # RCT/RD -> Obs
  expect_error(as_rd_spec(specrct))
  check_equiv_spec(specrd, as_rd_spec(specrct, data = simdata,
                                      forcing = ~ . + forcing(force)))
  check_equiv_spec(specrd, as_rd_spec(specrct, data = simdata,
                                      forcing = . ~ . + forcing(force)))
  check_equiv_spec(specrd, as_rd_spec(specobs, data = simdata,
                                      forcing = ~ . + forcing(force)))
  check_equiv_spec(specrd, as_rd_spec(specobs, data = simdata,
                                       forcing = . ~ . + forcing(force)))
  expect_error(as_rd_spec(specobs, data = simdata, forcing = simdata$force),
               "as a formula")

})


test_that("variable transformations in StudySpecification", {

  data(simdata)

  spec <- rct_spec(dose + 3 ~ uoa(uoa1, uoa2), data = simdata)
  expect_true(all(treatment(spec)[, 1] %in% (simdata$dose + 3)))

  spec <- rd_spec(force > 4 & bid < 2 ~ cluster(uoa1, uoa2) + forcing(force),
                   data = simdata)
  expect_true(is.logical(treatment(spec)[, 1]))
})

test_that("#26 cbind in identification variables", {
  data(simdata)
  spec1 <- obs_spec(o ~ cluster(cbind(uoa1, uoa2)), data = simdata)
  spec2 <- obs_spec(o ~ cluster(uoa1, uoa2), data = simdata)

  spec1@call <- fc
  spec2@call <- fc
  expect_identical(spec1, spec2)

})

test_that("variable_concordance", {

  data(simdata)

  # No issues
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  expect_silent(expect_true(specification_data_concordance(spec, simdata)))

  # inconsistent block
  sd <- simdata
  sd$bid[1] <- 2

  expect_warning(a <- specification_data_concordance(spec, sd),
                 "Inconsistencies.*`bid`")
  expect_false(a)

  # Inconsistent block & force
  sd$force[2] <- 2

  expect_warning(expect_warning(a <- specification_data_concordance(spec, sd),
                                "Inconsistencies.*`bid`"),
                  "Inconsistencies.*`force`")
  expect_false(a)

  # inconsistent treatment
  sd <- simdata
  sd$z[1] <- 1

  expect_warning(a <- specification_data_concordance(spec, sd),
                 "Inconsistencies.*`z`")

  # Missing variables
  sd <- simdata
  sd$z <- NULL
  expect_warning(a <- specification_data_concordance(spec, sd),
                 "`z` not found")
  expect_false(a)
  expect_silent(a <- specification_data_concordance(spec, sd,
                                             warn_on_nonexistence = FALSE))
  expect_true(a)

  sd <- simdata
  sd$bid <- NULL
  expect_warning(a <- specification_data_concordance(spec, sd),
                 "`bid` not found")
  expect_false(a)
  expect_silent(a <- specification_data_concordance(spec, sd,
                                             warn_on_nonexistence = FALSE))
  expect_true(a)

  sd <- simdata
  sd$force <- NULL
  expect_warning(a <- specification_data_concordance(spec, sd),
                 "`force` not found")
  expect_false(a)
  expect_silent(a <- specification_data_concordance(spec, sd,
                                             warn_on_nonexistence = FALSE))
  expect_true(a)

  sd$uoa1 <- NULL
  expect_error(specification_data_concordance(spec, sd),
               "missing a `by=`")

  # No block
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force),
                   data = simdata)
  expect_silent(expect_true(specification_data_concordance(spec, simdata)))

  # multiple variables
  sd <- simdata
  sd$force2 <- sd$force + 1
  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force, force2),
                   data = sd)
  expect_silent(expect_true(specification_data_concordance(spec, sd)))

  expect_warning(specification_data_concordance(spec, simdata),
                 "`force2` not found")

  sd$force[1] <- 1
  expect_warning(specification_data_concordance(spec, sd),
                 "Inconsistencies.*`force`$")

  sd2 <- simdata
  sd2$force[1] <- 1
  expect_warning(expect_warning(specification_data_concordance(spec, sd2),
                                "Inconsistencies.*`force`$"),
                 "`force2` not found")

  sd$force2[1] <- 2
  expect_warning(expect_warning(specification_data_concordance(spec, sd),
                                "Inconsistencies.*`force`$"),
                 "Inconsistencies.*`force2`$")

  data(simdata)

  spec <- rd_spec(z ~ cluster(uoa1, uoa2) + block(bid) + forcing(force),
                   data = simdata)

  names(simdata)[names(simdata) == "uoa1"] <- "ccc"
  names(simdata)[names(simdata) == "force"] <- "fff"

  expect_silent(expect_true({
    specification_data_concordance(spec, simdata,
                            by = list(force = "fff",
                                      uoa1 = "ccc"))
  }))

})

test_that("NA's in data to create StudySpecification", {
  data(simdata)
  sd2 <- simdata
  sd2$uoa1[1] <- NA

  expect_error(obs_spec(z ~ unitid(uoa1, uoa2), data = sd2),
               "Missing values cannot be found")

  spec1 <- obs_spec(z ~ unitid(uoa1, uoa2), data = simdata)
  spec2 <- obs_spec(z ~ unitid(uoa1, uoa2), data = sd2, na.fail = FALSE)

  spec1@call <- call("c")
  spec2@call <- call("c")

  expect_identical(spec1, spec2)

  sd2 <- simdata
  sd2$bid <- NA

  expect_silent(spec <- obs_spec(z ~ unitid(uoa1, uoa2), data = sd2))
  expect_error(obs_spec(z ~ unitid(uoa1, uoa2) + block(bid), data = sd2),
               "Missing values cannot be found")

  # treatment NAs allowed & passed forward
  sd2 <- simdata
  sd2$z[1:10] <- NA
  expect_silent(spec <- obs_spec(z ~ unitid(uoa1, uoa2), data = sd2,
                                  na.fail = TRUE))
  expect_equal(nrow(spec@structure), 10)
  expect_silent(spec <- obs_spec(z ~ unitid(uoa1, uoa2), data = sd2,
                                  na.fail = FALSE))
  expect_equal(nrow(spec@structure), 10)

  sd2 <- simdata
  sd2$uoa1[1:4] <- NA
  spec <- obs_spec(z ~ unitid(uoa1, uoa2), data = sd2, na.fail = FALSE)
  expect_equal(nrow(spec@structure), 9)

  # Cases with NAs in treatment and additional structure
  # variable(s) silently omitted
  sd3  <- sd2
  sd3$z[1:4]  <- NA
  expect_silent(spec_  <- obs_spec(z ~ unitid(uoa1, uoa2), data = sd3,
                                   na.fail = TRUE))
  expect_equal(nrow(spec_@structure), 9)

  spec@call <- spec_@call <- call("c")
  expect_identical(spec, spec_)
})

test_that("column name 'cluster' doesn't cause issues", {
  data(simdata)
  simdata$cluster <- simdata$uoa1

  expect_no_error(rct_spec(z ~ cluster(cluster, uoa2), data = simdata))
  expect_true(TRUE) # Avoid an empty test warning
})

test_that("#209 fix subsetting", {
  data(simdata)
  spec1 <- rct_spec(z ~ uoa(uoa1,uoa2), subset= simdata$bid != 1, data=simdata)
  spec2 <- rct_spec(z ~ uoa(uoa1,uoa2), subset= bid != 1, data=simdata)
  simdatasub <- subset(simdata, bid != 1)
  spec3 <- rct_spec(z ~ uoa(uoa1,uoa2), data=simdatasub)
  expect_identical(spec1@structure,
                   spec2@structure)
  expect_identical(spec1@structure,
                   spec3@structure)

  spec1 <- obs_spec(z ~ uoa(uoa1,uoa2), subset= simdata$bid != 1, data=simdata)
  spec2 <- obs_spec(z ~ uoa(uoa1,uoa2), subset= bid != 1, data=simdata)
  simdatasub <- subset(simdata, bid != 1)
  spec3 <- obs_spec(z ~ uoa(uoa1,uoa2), data=simdatasub)
  expect_identical(spec1@structure,
                   spec2@structure)
  expect_identical(spec1@structure,
                   spec3@structure)
})

test_that("#222 factor UOA", {

  example <- data.frame(schf = as.factor(letters),
                        schn = 1:26,
                        trt = rep(c(0,1),each = 13),
                        Y = rnorm(26))
  sp1 <- rct_spec(trt ~ unitid(schf), data = example)
  sp2 <- rct_spec(trt ~ unitid(schn), data = example)
  mod1 <- lmitt(Y ~ 1,
                specification = sp1,
                data = example,weights = "ate")
  mod2 <- lmitt(Y ~ 1,
                specification = sp2,
                data = example,weights = "ate")
  expect_identical(vcov(mod1), vcov(mod2))

})

test_that("get StudySpecification data from call's formula environment", {
  make_spec <- function(udata, zcol, ucol) {
    spec_form <- reformulate(paste0("unitid(", ucol, ")"), response = zcol)
    rct_specification(spec_form, udata)
  }
  
  udata <- data.frame(uid = seq_len(3), z = c(0, 1, 0))
  m1 <- make_spec(udata, "z", "uid")
  expect_identical(udata, get("udata", environment(m1@call$formula)))
  expect_silent(.make_uoa_cluster_df(m1, "uid"))
  
})
