save_options <- options()
options("propertee_message_on_unused_blocks" = FALSE)

test_that("lmitt", {

  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  mod <- lmitt(y ~ 1, weights = ate(), data = simdata, specification = spec)

  expect_s4_class(mod, "teeMod")
  expect_true(inherits(mod, "lm"))

  mod_ett <- lmitt(y ~ 1, weights = ett(), data = simdata, specification = spec)

  expect_s4_class(mod_ett, "teeMod")
})

test_that("lmitt and lm return the same in simple cases", {
  data(simdata)
  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid),
                   data = simdata)

  ml <- lm(y ~ assigned(), data = simdata, weights = ate(spec))
  ssmod <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)
  ml2ssmod <- lmitt(ml)

  expect_true(all.equal(ssmod$coef, ml$coef, check.attributes = FALSE))
  expect_true(all.equal(ssmod$coef, ml2ssmod$coef, check.attributes = FALSE))
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(ssmod$model$`(weights)`, slot)),
                           deparse1(methods::slot(ml$model$`(weights)`, slot)))
               } else {
                 identical(methods::slot(ssmod$model$`(weights)`, slot),
                           methods::slot(ml$model$`(weights)`, slot))
               },
               logical(1L)))
  )
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(ssmod$model$`(weights)`, slot)),
                           deparse1(methods::slot(ml2ssmod$model$`(weights)`, slot)))
               } else {
                 identical(methods::slot(ssmod$model$`(weights)`, slot),
                           methods::slot(ml2ssmod$model$`(weights)`, slot))
               },
               logical(1L)))
  )

  ml <- lm(y ~ z, data = simdata, weights = ett(spec))
  ssmod <- lmitt(y ~ 1, data = simdata, weights = ett(), specification = spec)

  expect_true(all.equal(ssmod$coef, ml$coef, check.attributes = FALSE))
  expect_true(
    all(vapply(c(".Data", "StudySpecification", "target", "dichotomy"),
               function(slot) if (slot == "dichotomy") {
                 identical(deparse1(methods::slot(ssmod$model$`(weights)`, slot)),
                           deparse1(methods::slot(ml$model$`(weights)`, slot)))
               } else {
                 identical(methods::slot(ssmod$model$`(weights)`, slot),
                           methods::slot(ml$model$`(weights)`, slot))
               },
               logical(1L)))
  )

})

test_that("covariate adjustment", {
  data(simdata)
  spec <- rd_spec(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  camod <- lm(y ~ x, data = simdata)

  ssmod <- lmitt(y ~ 1, data = simdata, weights = ate(),
              offset = cov_adj(camod), specification = spec)
  expect_true(!is.null(ssmod$model$"(offset)"))
  expect_true(inherits(ssmod$model$"(offset)", "SandwichLayer"))

})


test_that("StudySpecification argument", {
  data(simdata)
  spec <- obs_spec(z ~ cluster(uoa1, uoa2), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec)
  mod3 <- lmitt(y ~ 1, data = simdata, specification = z ~ cluster(uoa1, uoa2))
  expect_true(mod3@StudySpecification@type == "Obs")
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_identical(mod1@StudySpecification, mod3@StudySpecification)

  spec2 <- rd_spec(z ~ cluster(uoa1, uoa2) + forcing(force), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec2)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec2)
  mod3 <- lmitt(y ~ 1, data = simdata,
                specification = z ~ cluster(uoa1, uoa2) + forcing(force))
  expect_true(mod3@StudySpecification@type == "RD")
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_identical(mod1@StudySpecification, mod3@StudySpecification)

  spec3 <- obs_spec(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec3,
                subset = simdata$dose < 300)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(),
                subset = simdata$dose < 300, specification = spec3)
  mod3 <- lmitt(y ~ 1, data = simdata,
                specification = z ~ cluster(uoa1, uoa2) + block(bid),
                subset = simdata$dose < 300)
  expect_true(mod3@StudySpecification@type == "Obs")
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_identical(mod1@StudySpecification, mod3@StudySpecification)

  expect_error(lmitt(y ~ 1, data = simdata, specification = 4),
               "formula specifying such a specification")

})

test_that("Dichotomy argument", {

  data(simdata)
  spec <- obs_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec, dichotomy = dose > 200 ~ .)
  mod2 <- lmitt(y ~ 1, data = simdata,
                specification = dose ~ cluster(uoa1, uoa2),
                dichotomy = dose > 200 ~ .)
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_true(all.equal(mod1$coefficients, mod2$coefficients,
                        check.attributes =  FALSE))

})

test_that("Dichotomy from weights", {
  data(simdata)
  spec <- obs_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec,
                weights = ate(dichotomy = dose > 200 ~ .))
  mod2 <- lmitt(y ~ 1, data = simdata, specification = spec, weights = "ate",
                dichotomy = dose > 200 ~ .)
  mod3 <- lmitt(y ~ 1, data = simdata, specification = spec, weights = ate(spec),
                dichotomy = dose > 200 ~ .)
  expect_identical(mod1@StudySpecification, mod2@StudySpecification)
  expect_identical(mod1@StudySpecification, mod3@StudySpecification)
  expect_true(all.equal(mod1$coefficients, mod2$coefficients,
                        check.attributes =  FALSE))
  expect_true(all.equal(mod1$coefficients, mod3$coefficients,
                        check.attributes =  FALSE))
  expect_true(
    all.equal(mod2$model$`(weights)`, ate(spec, data = simdata, dichotomy = dose > 200 ~ .))
  )
  expect_true(
    all.equal(mod3$model$`(weights)`, ate(spec, data = simdata, dichotomy = dose > 200 ~ .))
  )
})

test_that("Minimal support for continuous treatments", {

  data(simdata)
  spec1 <- obs_spec(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, specification = spec1, absorb=FALSE)
  expect_true(any(!model.frame(mod1)[, 2] %in% 0:1))
  expect_silent(lmitt(y ~ x, data = simdata, specification = spec1))
  expect_silent(lmitt(y ~ as.character(o), data = simdata, specification = spec1))

  spec2 <- obs_spec(dose ~ cluster(uoa1, uoa2) +block(bid), data = simdata)
  mod2 <- lmitt(y ~ 1, data = simdata, specification = spec2, absorb=TRUE)
  expect_true(any(!model.frame(mod2)[, 2] %in% 0:1))
  ## But we don't also allow factor or ordinal treatments
  simdata_ <-
      transform(simdata, dosef=cut(dose, breaks=c(50, 100, 200, 300),
                          labels=c("lo", "med", "hi"),
                          include.lowest=TRUE)
                )
  spec3 <-
      obs_spec(dosef ~ cluster(uoa1, uoa2), data = simdata_)
  expect_error(lmitt(y ~ 1, data = simdata, specification = spec3),
               "not supported")
})

test_that("lmitt finds StudySpecification wherever it's stored", {
  data(simdata)
  spec_form <- z ~ uoa(uoa1, uoa2) + block(bid)
  spec <- obs_spec(spec_form, data = simdata)
  camod <- lm(y ~ x, data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, weights = ate(),
                offset = cov_adj(camod, specification = spec),
                specification = spec)
  mod2 <- lmitt(y ~ 1, specification = spec_form, data = simdata, weights = ate(),
                offset = cov_adj(camod))
  mod3 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate(),
                offset = cov_adj(camod))

  expect_equal(coef(mod1), coef(mod2))
  expect_equal(coef(mod1), coef(mod3))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1), vcov_tee(mod3))
})

test_that("Allowed inputs to lmitt #73", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  ### Allowed versions

  expect_no_error(l1 <- lmitt(y ~ 1, specification = spec, data = simdata))
  expect_no_error(l2 <- lmitt(y ~ dose, specification = spec, data = simdata))

  ### Disallowed versions

  expect_error(lmitt(y ~ 0, specification = spec, data = simdata))
  expect_error(lmitt(y ~ dose + 0, specification = spec, data = simdata))
  expect_error(lmitt(y ~ 0 + dose, specification = spec, data = simdata))
  expect_error(lmitt(y ~ dose + bid, specification = spec, data = simdata))
  expect_error(lmitt(y ~ 1 + dose + bid, specification = spec, data = simdata))
  expect_error(lmitt(y ~ 1 + x, specification = spec, data = simdata))
  expect_error(lmitt(y ~ x + 1, specification = spec, data = simdata))

  ### Main + Subgroup effects.
  ### NYI - see #73
  expect_error(lmitt(y ~ 1 + dose, specification = spec, data = simdata),
               "To estimate subgroup effects") # Remove this once #73 is fixed

  #expect_no_error(l3 <- lmitt(y ~ 1 + dose, specification = spec, data = simdata))
  #expect_no_error(l4 <- lmitt(y ~ dose + 1, specification = spec, data = simdata))

  #expect_identical(l3$coeff, l4$coeff)

})

test_that("weights argument can be string", {

  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  l1 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ate())
  l2 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate")

  l3 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = ett())
  l4 <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ett")

  l1@lmitt_call <- call("ls")
  l2@lmitt_call <- call("ls")
  l3@lmitt_call <- call("ls")
  l4@lmitt_call <- call("ls")

  expect_true(all.equal(l1, l2))
  expect_true(all.equal(l3, l4))

  # all.equal returns strings detailing differences
  expect_true(is.character(all.equal(l1, l3)))

  # Passing a different string warns users, but allows `lm()` to error
  # in case user is trying to pass something to someplace else
  expect_error(expect_warning(lmitt(y ~ 1, specification = spec, data =simdata,
                                    weights = "abc"), "other than \"ate\""))

})

test_that("Regular weights can still be used", {

  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ 1, specification = spec, weights = simdata$dose, data = simdata)
  expect_true(is(mod, "teeMod"))

})

test_that("Aliases for assigned() aren't allowed either", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_error(lmitt(y ~ assigned(), specification = spec, data = simdata),
               "Do not specify")

  expect_error(lmitt(y ~ adopters(), specification = spec, data = simdata),
               "Do not specify")

  expect_error(lmitt(y ~ a.(), specification = spec, data = simdata),
               "Do not specify")


  expect_error(lmitt(y ~ z.(), specification = spec, data = simdata),
               "Do not specify")

  #### because of the `.` in `a.()`, let's make sure the regexp isn't breaking
  ad <- assigned
  # Passing in a different weight to ensure we get a predictable error that
  # ocurs *after* checking the formula
  spec2 <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  expect_error(lmitt(y ~ ad(), specification = spec, data = simdata,
                     weights = ate(spec2)),
               "Multiple differing")

  # Just to ensure that `lmitt.formula()` doesn't get changed in a way such that
  # checking for differing specifications happens prior:
  expect_error(lmitt(y ~ a.(), specification = spec, data = simdata,
                     weights = ate(spec2)),
               "Do not specify")

  rm(ad)

})

test_that("User can pass WeightedStudySpecification", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ 1, data = simdata, specification = ate(spec))

  expect_identical(spec, mod@StudySpecification)

})

test_that("#82 Informative error if specification in weight/offset but not lmitt", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_error(lmitt(y ~ 1, data = simdata, weights = ate(spec)),
               "into the weight function")

  cmod <- lm(y ~ x, data = simdata)
  expect_error(lmitt(y ~ 1, data = simdata, offset = cov_adj(cmod, specification = spec)),
               "into the offset function")

})

test_that("#81 continuous moderator", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) , data = simdata)

  mod <- lmitt(y ~ x, data = simdata, specification = spec)
  expect_length(mod$coeff, 4)

  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ x, data = simdata, specification = spec, absorb = TRUE)
  expect_length(mod$coeff, 4)

})

test_that("non-integer units of assignment", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, specification = spec, absorb = TRUE))

  simdata$uoa1 <- as.character(simdata$uoa1)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, specification = spec, absorb = TRUE))

  simdata$uoa1 <- as.factor(simdata$uoa1)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, specification = spec, absorb = TRUE))

  expect_true(TRUE) # Avoid an empty test warning

})

test_that("#137 ensure absorb is using the correct moderator", {
  data(simdata)

  mod1 <- lmitt(y ~ o, data = simdata,
                specification = z ~ cluster(uoa1, uoa2) + block(bid))
  mod2 <- lmitt(y ~ o, data = simdata,
                specification = z ~ cluster(uoa1, uoa2) + block(bid),
                absorb = TRUE)

  o1 <- model.matrix(mod1)[, "o"]
  o2 <- model.matrix(mod2)[, "o"]

  # o1 and o2 should be different, since the latter should have been
  # group-centered due to `absorb = TRUE`. Due to a bug discovered in #137, the
  # un-centered version of the continuous moderator was being using previously,
  # leading o1 and o2 to be identical (incorrectly).

  expect_true(!(all.equal(o1, o2, check.attributes = FALSE) == TRUE))


})

test_that("#140 handling 0 weights", {

  # A single 0 weight observation in a larger block - that observation's
  # conbribution to the estimating equation should be 0
  data(simdata)
  simdata$weight <- 1
  simdata$weight[1] <- 0

  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  mod <- lmitt(y ~ x, data = simdata, specification = spec,
               weights = simdata$weight)
  ee <- propertee:::estfun.teeMod(mod)
  expect_true(all(ee[1, ] == 0))

  # A block with only a single non-zero weight should contribute nothing to the
  # estimating equation
  data(simdata)
  simdata$weight <- 1
  simdata$weight[simdata$bid == 1] <- 0
  simdata$weight[simdata$bid == 1][1] <- 1

  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  mod <- lmitt(y ~ x, data = simdata, specification = spec,
               absorb = TRUE, weights = simdata$weight)
  ee <- propertee:::estfun.teeMod(mod)
  expect_true(all(ee[simdata$bid == 1, 2:ncol(ee)] == 0))
  expect_true(all(ee[simdata$bid == 1 & simdata$weight == 0, 1] == 0))

  data(simdata)
  simdata$z[simdata$bid == 1] <- 1
  simdata$weight <- 1
  spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  # Block 1 has 0 variance in treatment
  mod <- lmitt(y ~ x, data = simdata, specification = spec,
               absorb = TRUE, weights = simdata$weight)
  ee <- estfun.teeMod(mod)
  expect_true(all(ee[simdata$bid == 1, "z."] == 0))
  expect_true(all(
    sapply(ee[simdata$bid == 1, c("(Intercept)", "x", "z._x")],
           function(vals) any(vals != 0))
  ))

})


options(save_options)
#### !!!!!!!!!!!NOTE!!!!!!!!!!!!!
# This test below should NOT have `options()$propertee_message_on_unused_blocks`
# set to FALSE. So it needs to stay below the restoration of options line above.
# Other tests should probably go above the restoration of options line.

test_that("Message if specification has block info but isn't used in lmitt", {
  data(simdata)
  spec <- obs_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_message(lmitt(y ~ 1, data = simdata, specification = spec),
                 "is not used in this model")
  expect_message(lmitt(y ~ 1, data = simdata, specification = ate(spec)),
                 "is not used in this model")
  expect_message(lmitt(y ~ 1, data = simdata, specification = ate(spec),
                       weights = abs(simdata$x)),
                 "is not used in this model")

  expect_silent(x <- lmitt(y ~ 1, data = simdata, absorb = TRUE, specification = spec))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = "ate", specification = spec))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = ate(), specification = spec))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = ett()*3, specification = spec))

})

### READ COMMENT ABOUT LAST TEST!

test_that(paste("absorb=TRUE doesn't drop rows when some strata have weights=0",
                "(some strata only have treated clusters)"), {
  data(simdata)
  simdata[simdata$uoa1 == 5 & simdata$uoa2 == 2, "bid"] <- 4

  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  absorb_mod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate", absorb = TRUE)
  no_absorb_mod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate", absorb = FALSE)

  expect_equal(dim(model.matrix(absorb_mod)), dim(model.matrix(no_absorb_mod)))
})

test_that(paste("absorb=TRUE doesn't drop rows when some strata have weights=0",
                "(some strata only have untreated clusters)"), {
  data(simdata)
  simdata[simdata$uoa1 == 4 & simdata$uoa2 == 2, "bid"] <- 4

  spec <- rct_spec(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  absorb_mod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate", absorb = TRUE)
  no_absorb_mod <- lmitt(y ~ 1, specification = spec, data = simdata, weights = "ate", absorb = FALSE)

  expect_equal(dim(model.matrix(absorb_mod)), dim(model.matrix(no_absorb_mod)))
})

test_that("#147 lmitt.formula accepts references to formula objects", {
  data(simdata)
  lmitt_form <- y ~ 1
  spec <- rct_spec(z ~ uoa(uoa1, uoa2), data = simdata)
  expect_equal(
    co <- capture.output(lmitt(lmitt_form, data = simdata, specification = spec)),
    capture.output(lmitt(y ~ 1, data = simdata, specification = spec))
  )

  make_lmitt <- function(data, dv_col) {
    lmitt_form <- as.formula(paste0(dv_col, "~1"))
    spec <- rct_spec(z ~ unitid(uoa1, uoa2), data)
    lmitt(lmitt_form, specification = spec, data = data)
  }
  expect_equal(capture.output(make_lmitt(simdata, "y")), co)

})
