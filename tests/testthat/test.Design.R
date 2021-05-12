context("Design object creation")

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
