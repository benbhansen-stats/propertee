test_that("Design creation", {
  d <- new("Design",
           structure = data.frame(a = 1, b = 2),
           columnIndex = c("t", "f"),
           type = "RCT")

  expect_s4_class(d, "Design")
  expect_s3_class(d@structure, "data.frame")
  expect_type(d@columnIndex, "character")
  expect_type(d@type, "character")

  expect_equal(ncol(d@structure), length(d@columnIndex))
  expect_length(d@type, 1)
})

test_that("Design validity", {
  expect_error(new("Design",
                   structure = data.frame(a = 1, b = 2),
                   columnIndex = "a",
                   type = "RCT"),
               "number of columns")

  expect_error(new("Design",
                   structure = data.frame(),
                   columnIndex = "a",
                   type = "RCT"),
               "positive dimensions")

  expect_error(new("Design",
                   structure = data.frame(a = 1, b = 2),
                   columnIndex = c("t", "f"),
                   type = "abc"),
               "unknown @type")

  expect_error(new("Design",
                   structure = data.frame(a = 1, b = 2),
                   columnIndex = c("q", "k"),
                   type = "RCT"),
               "unknown elements")

})

test_that("Design creation", {
  data(mtcars)
  mtcars <- mtcars[-c(5, 11),]

  d_rct <- New_Design(vs ~ cluster(qsec), data = mtcars, type = "RCT")

  expect_s4_class(d_rct, "Design")
  expect_s3_class(d_rct@structure, "data.frame")
  expect_type(d_rct@columnIndex, "character")
  expect_type(d_rct@type, "character")

  expect_equal(dim(d_rct@structure), c(30, 2))
  expect_length(d_rct@columnIndex, 2)
  expect_length(d_rct@type, 1)

  expect_equal(d_rct@structure, subset(mtcars, select = c("vs", "qsec")))
  expect_equal(d_rct@columnIndex, c("t", "c"))
  expect_equal(d_rct@type, "RCT")

  # subset

  d_obs <- New_Design(vs ~ cluster(qsec), data = mtcars, type = "Obs",
                  subset = mtcars$mpg > 17)

  expect_equal(dim(d_obs@structure), c(sum(mtcars$mpg >  17), 2))
  expect_length(d_obs@columnIndex, 2)
  expect_length(d_obs@type, 1)

  expect_equal(d_obs@structure, subset(mtcars, select = c("vs", "qsec"),
                                       subset = mtcars$mpg > 17))
  expect_equal(d_obs@columnIndex, c("t", "c"))
  expect_equal(d_obs@type, "Obs")

  ### Complex design
  d_rd <- New_Design(vs ~ block(disp, gear) + forcing(wt, cyl) + cluster(mpg, qsec),
                  data = mtcars, type = "RD")

  expect_equal(d_rd@structure, subset(mtcars, select = c("vs", "disp", "gear",
                                                      "wt", "cyl", "mpg", "qsec")))
  expect_equal(d_rd@columnIndex, c("t", "b", "b", "f", "f", "c", "c"))
  expect_equal(d_rd@type, "RD")

  ### Specific designs
  rct_des <- RCT_Design(vs ~ cluster(qsec), data = mtcars)
  expect_identical(d_rct, rct_des)

  obs_des <- Obs_Design(vs ~ cluster(qsec), data = mtcars,
                       subset = mtcars$mpg > 17)
  expect_identical(d_obs, obs_des)

  rd_des <- RD_Design(vs ~ block(disp, gear) + forcing(wt, cyl) + cluster(mpg, qsec),
                  data = mtcars)
  expect_identical(d_rd, rd_des)
})

test_that("unit of assignment differs from unit of analysis", {

  data(simdata)

  desrct <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  expect_s4_class(desrct, "Design")
  expect_equal(nrow(desrct@structure), 10)

  expect_error(RCT_Design(z ~ cluster(cid1) + block(bid), data = simdata),
               "must be constant")


})

test_that("Design printing", {
  data(simdata)

  desrct <- RCT_Design(z ~ cluster(cid1, cid2) + block(bid), data = simdata)
  desrd  <-  RD_Design(z ~ cluster(cid1, cid2) + block(bid) + forcing(force), data = simdata)
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
