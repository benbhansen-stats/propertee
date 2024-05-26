save_options <- options()
options("propertee_message_on_unused_blocks" = FALSE)

test_that("lmitt", {

  data(simdata)
  des <- rd_design(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  da <- lmitt(y ~ 1, weights = ate(), data = simdata, design = des)

  expect_s4_class(da, "teeMod")
  expect_true(inherits(da, "lm"))

  da_ett <- lmitt(y ~ 1, weights = ett(), data = simdata, design = des)

  expect_s4_class(da_ett, "teeMod")
})

test_that("lmitt and lm return the same in simple cases", {
  data(simdata)
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid),
                   data = simdata)

  ml <- lm(y ~ assigned(), data = simdata, weights = ate(des))
  da <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)
  ml2da <- lmitt(ml)

  expect_true(all.equal(da$coef, ml$coef, check.attributes = FALSE))
  expect_true(all.equal(da$coef, ml2da$coef, check.attributes = FALSE))
  expect_identical(da$model$`(weights)`, ml$model$`(weights)`)
  expect_identical(da$model$`(weights)`, ml2da$model$`(weights)`)

  ml <- lm(y ~ z, data = simdata, weights = ett(des))
  da <- lmitt(y ~ 1, data = simdata, weights = ett(), design = des)

  expect_true(all.equal(da$coef, ml$coef, check.attributes = FALSE))
  expect_identical(da$model$`(weights)`, ml$model$`(weights)`)

})

test_that("covariate adjustment", {
  data(simdata)
  des <- rd_design(z ~ cluster(uoa2, uoa1) + block(bid) + forcing(force),
                   data = simdata)

  camod <- lm(y ~ x, data = simdata)

  da <- lmitt(y ~ 1, data = simdata, weights = ate(),
              offset = cov_adj(camod), design = des)
  expect_true(!is.null(da$model$"(offset)"))
  expect_true(inherits(da$model$"(offset)", "SandwichLayer"))

})


test_that("Design argument", {
  data(simdata)
  des <- obs_design(z ~ cluster(uoa1, uoa2), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, design = des)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des)
  mod3 <- lmitt(y ~ 1, data = simdata, design = z ~ cluster(uoa1, uoa2))
  expect_true(mod3@Design@type == "Obs")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

  des2 <- rd_design(z ~ cluster(uoa1, uoa2) + forcing(force), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, design = des2)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des2)
  mod3 <- lmitt(y ~ 1, data = simdata,
                design = z ~ cluster(uoa1, uoa2) + forcing(force))
  expect_true(mod3@Design@type == "RD")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

  des3 <- obs_design(z ~ cluster(uoa1, uoa2) + block(bid), data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, design = des3,
                subset = simdata$dose < 300)
  mod2 <- lmitt(y ~ 1, data = simdata, weights = ate(),
                subset = simdata$dose < 300, design = des3)
  mod3 <- lmitt(y ~ 1, data = simdata,
                design = z ~ cluster(uoa1, uoa2) + block(bid),
                subset = simdata$dose < 300)
  expect_true(mod3@Design@type == "Obs")
  expect_identical(mod1@Design, mod2@Design)
  expect_identical(mod1@Design, mod3@Design)

  expect_error(lmitt(y ~ 1, data = simdata, design = 4),
               "formula specifying such a design")

})

test_that("Dichotomy argument", {

  data(simdata)
  des <- obs_design(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des, dichotomy = dose > 200 ~ .)
  mod2 <- lmitt(y ~ 1, data = simdata,
                design = dose ~ cluster(uoa1, uoa2),
                dichotomy = dose > 200 ~ .)
  expect_identical(mod1@Design, mod2@Design)
  expect_true(all.equal(mod1$coefficients, mod2$coefficients,
                        check.attributes =  FALSE))

})

test_that("Dichotomy from weights", {
  data(simdata)
  des <- obs_design(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des,
                weights = ate(dichotomy = dose > 200 ~ .))
  mod2 <- lmitt(y ~ 1, data = simdata, design = des, weights = "ate",
                dichotomy = dose > 200 ~ .)
  expect_identical(mod1@Design, mod2@Design)
  expect_true(all.equal(mod1$coefficients, mod2$coefficients,
                        check.attributes =  FALSE))
  expect_true(
    all.equal(mod2$model$`(weights)`, ate(des, data = simdata, dichotomy = dose > 200 ~ .))
  )
})

test_that("Allow non-binary treatment", {

  data(simdata)
  des <- obs_design(dose ~ cluster(uoa1, uoa2), data = simdata)
  mod1 <- lmitt(y ~ 1, data = simdata, design = des)
  expect_true(any(!model.frame(mod1)[, 2] %in% 0:1))

})

test_that("lmitt finds Design wherever it's stored", {
  data(simdata)
  des_form <- z ~ uoa(uoa1, uoa2) + block(bid)
  des <- obs_design(des_form, data = simdata)
  camod <- lm(y ~ x, data = simdata)

  mod1 <- lmitt(y ~ 1, data = simdata, weights = ate(),
                offset = cov_adj(camod, design = des),
                design = des)
  mod2 <- lmitt(y ~ 1, design = des_form, data = simdata, weights = ate(),
                offset = cov_adj(camod))
  mod3 <- lmitt(y ~ 1, design = des, data = simdata, weights = ate(),
                offset = cov_adj(camod))

  expect_equal(coef(mod1), coef(mod2))
  expect_equal(coef(mod1), coef(mod3))
  expect_equal(vcov_tee(mod1), vcov_tee(mod2))
  expect_equal(vcov_tee(mod1), vcov_tee(mod3))
})

test_that("Allowed inputs to lmitt #73", {
  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  ### Allowed versions

  expect_no_error(l1 <- lmitt(y ~ 1, design = des, data = simdata))
  expect_no_error(l2 <- lmitt(y ~ dose, design = des, data = simdata))

  ### Disallowed versions

  expect_error(lmitt(y ~ 0, design = des, data = simdata))
  expect_error(lmitt(y ~ dose + 0, design = des, data = simdata))
  expect_error(lmitt(y ~ 0 + dose, design = des, data = simdata))
  expect_error(lmitt(y ~ dose + bid, design = des, data = simdata))
  expect_error(lmitt(y ~ 1 + dose + bid, design = des, data = simdata))
  expect_error(lmitt(y ~ 1 + x, design = des, data = simdata))
  expect_error(lmitt(y ~ x + 1, design = des, data = simdata))

  ### Main + Subgroup effects.
  ### NYI - see #73
  expect_error(lmitt(y ~ 1 + dose, design = des, data = simdata),
               "To estimate subgroup effects") # Remove this once #73 is fixed

  #expect_no_error(l3 <- lmitt(y ~ 1 + dose, design = des, data = simdata))
  #expect_no_error(l4 <- lmitt(y ~ dose + 1, design = des, data = simdata))

  #expect_identical(l3$coeff, l4$coeff)

})

test_that("weights argument can be string", {

  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  l1 <- lmitt(y ~ 1, design = des, data = simdata, weights = ate())
  l2 <- lmitt(y ~ 1, design = des, data = simdata, weights = "ate")

  l3 <- lmitt(y ~ 1, design = des, data = simdata, weights = ett())
  l4 <- lmitt(y ~ 1, design = des, data = simdata, weights = "ett")

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
  expect_error(expect_warning(lmitt(y ~ 1, design = des, data =simdata,
                                    weights = "abc"), "other than \"ate\""))

})

test_that("Regular weights can still be used", {

  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ 1, design = des, weights = simdata$dose, data = simdata)
  expect_true(is(mod, "teeMod"))

})

test_that("Aliases for assigned() aren't allowed either", {
  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_error(lmitt(y ~ assigned(), design = des, data = simdata),
               "Do not specify")

  expect_error(lmitt(y ~ adopters(), design = des, data = simdata),
               "Do not specify")

  expect_error(lmitt(y ~ a.(), design = des, data = simdata),
               "Do not specify")


  expect_error(lmitt(y ~ z.(), design = des, data = simdata),
               "Do not specify")

  #### because of the `.` in `a.()`, let's make sure the regexp isn't breaking
  ad <- assigned
  # Passing in a different weight to ensure we get a predictable error that
  # ocurs *after* checking the formula
  des2 <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  expect_error(lmitt(y ~ ad(), design = des, data = simdata,
                     weights = ate(des2)),
               "Multiple differing")

  # Just to ensure that `lmitt.formula()` doesn't get changed in a way such that
  # checking for differing designs happens prior:
  expect_error(lmitt(y ~ a.(), design = des, data = simdata,
                     weights = ate(des2)),
               "Do not specify")

  rm(ad)

})

test_that("User can pass WeightedDesign", {
  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ 1, data = simdata, design = ate(des))

  expect_identical(des, mod@Design)

})

test_that("#82 Informative error if design in weight/offset but not lmitt", {
  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_error(lmitt(y ~ 1, data = simdata, weights = ate(des)),
               "into the weight function")

  cmod <- lm(y ~ x, data = simdata)
  expect_error(lmitt(y ~ 1, data = simdata, offset = cov_adj(cmod, design = des)),
               "into the offset function")

})

test_that("#81 continuous moderator", {
  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) , data = simdata)

  mod <- lmitt(y ~ x, data = simdata, design = des)
  expect_length(mod$coeff, 4)

  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  mod <- lmitt(y ~ x, data = simdata, design = des, absorb = TRUE)
  expect_length(mod$coeff, 4)

})

test_that("non-integer units of assignment", {
  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, design = des, absorb = TRUE))

  simdata$uoa1 <- as.character(simdata$uoa1)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, design = des, absorb = TRUE))

  simdata$uoa1 <- as.factor(simdata$uoa1)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid) , data = simdata)

  expect_no_error(lmitt(y ~ x, data = simdata, design = des, absorb = TRUE))


})

test_that("#137 ensure absorb is using the correct moderator", {
  data(simdata)

  mod1 <- lmitt(y ~ o, data = simdata,
                design = z ~ cluster(uoa1, uoa2) + block(bid))
  mod2 <- lmitt(y ~ o, data = simdata,
                design = z ~ cluster(uoa1, uoa2) + block(bid),
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

  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  mod <- lmitt(y ~ x, data = simdata, design = des,
               weights = simdata$weight)
  ee <- propertee:::estfun.teeMod(mod)
  expect_true(all(ee[1, ] == 0))

  # A block with only a single non-zero weight should contribute nothing to the
  # estimating equation
  data(simdata)
  simdata$weight <- 1
  simdata$weight[simdata$bid == 1] <- 0
  simdata$weight[simdata$bid == 1][1] <- 1

  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  mod <- lmitt(y ~ x, data = simdata, design = des,
               absorb = TRUE, weights = simdata$weight)
  ee <- propertee:::estfun.teeMod(mod)
  expect_true(all(ee[simdata$bid == 1, 2:ncol(ee)] == 0))
  expect_true(all(ee[simdata$bid == 1 & simdata$weight == 0, 1] == 0))

  data(simdata)
  simdata$z[simdata$bid == 1] <- 1
  simdata$weight <- 1
  des <- rct_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
  # Block 1 has 0 variance in treatment
  mod <- lmitt(y ~ x, data = simdata, design = des,
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

test_that("Message if design has block info but isn't used in lmitt", {
  data(simdata)
  des <- obs_design(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)

  expect_message(lmitt(y ~ 1, data = simdata, design = des),
                 "is not used in this model")
  expect_message(lmitt(y ~ 1, data = simdata, design = ate(des)),
                 "is not used in this model")
  expect_message(lmitt(y ~ 1, data = simdata, design = ate(des),
                       weights = abs(simdata$x)),
                 "is not used in this model")

  expect_silent(x <- lmitt(y ~ 1, data = simdata, absorb = TRUE, design = des))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = "ate", design = des))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = ate(), design = des))
  expect_silent(x <- lmitt(y ~ 1, data = simdata, weights = ett()*3, design = des))

})

### READ COMMENT ABOUT LAST TEST!

test_that(paste("absorb=TRUE doesn't drop rows when some strata have weights=0",
                "(some strata only have treated clusters)"), {
  data(simdata)
  simdata[simdata$uoa1 == 5 & simdata$uoa2 == 2, "bid"] <- 4
  
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  absorb_mod <- lmitt(y ~ 1, design = des, data = simdata, weights = "ate", absorb = TRUE)
  no_absorb_mod <- lmitt(y ~ 1, design = des, data = simdata, weights = "ate", absorb = FALSE)
  
  expect_equal(dim(model.matrix(absorb_mod)), dim(model.matrix(no_absorb_mod)))
})

test_that(paste("absorb=TRUE doesn't drop rows when some strata have weights=0",
                "(some strata only have untreated clusters)"), {
  data(simdata)
  simdata[simdata$uoa1 == 4 & simdata$uoa2 == 2, "bid"] <- 4
  
  des <- rct_design(z ~ cluster(uoa1, uoa2) + block(bid), simdata)
  absorb_mod <- lmitt(y ~ 1, design = des, data = simdata, weights = "ate", absorb = TRUE)
  no_absorb_mod <- lmitt(y ~ 1, design = des, data = simdata, weights = "ate", absorb = FALSE)
  
  expect_equal(dim(model.matrix(absorb_mod)), dim(model.matrix(no_absorb_mod)))
})

test_that("#147 lmitt.formula accepts references to formula objects", {
  data(simdata)
  lmitt_form <- y ~ 1
  des <- rct_design(z ~ uoa(uoa1, uoa2), data = simdata)
  expect_equal(
    co <- capture.output(lmitt(lmitt_form, data = simdata, design = des)),
    capture.output(lmitt(y ~ 1, data = simdata, design = des))
  )
  
  make_lmitt <- function(data, dv_col) {
    lmitt_form <- as.formula(paste0(dv_col, "~1"))
    des <- rct_design(z ~ unitid(uoa1, uoa2), data)
    lmitt(lmitt_form, design = des, data = data)
  }
  expect_equal(capture.output(make_lmitt(simdata, "y")), co)
  
})