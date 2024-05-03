# Issue 164
# https://github.com/benbhansen-stats/propertee/issues/164
#**********************************************************************
# If there's a change to any of these functions, uncomment and
# run these lines of code.
#
# They require the **digest** package. They compute the hash of
# each function and stores it for future comparison
#**********************************************************************
## orig <- capture.output(dput(stats::expand.model.frame))
## saveRDS(orig,
##         file = "tests/testthat/hashed_functions/expand.model.frame.rds")
## orig <- capture.output(dput(confint.lm))
## saveRDS(orig,
##         file = "tests/testthat/hashed_functions/confint.lm.rds")
#**********************************************************************
# Note that the path above assumes you are in the main folder. Below in the
# actual tests, the paths are relative to tests/testthat
#**********************************************************************


test_that("expand.model.frame", {
  stored <- readRDS("hashed_functions/expand.model.frame.rds")

  current <- capture.output(dput(stats::expand.model.frame))

  diff_lines <- which(stored != current)
  if (length(diff_lines) > 0) {
    # If there are differences, fail the test and show the different lines
    fail(paste("expand.model.frame has changed. Lines numbers differing from *original* are:\n",
               paste(diff_lines, collapse = ", ")))
  }
  expect_true(TRUE)
})

test_that("confint.lm", {
  stored <- readRDS("hashed_functions/confint.lm.rds")

  current <- capture.output(dput(stats::confint.lm))

  diff_lines <- which(stored != current)
  if (length(diff_lines) > 0) {
    # If there are differences, fail the test and show the different lines
    fail(paste("expand.model.frame has changed. Lines numbers differing from *original* are:\n",
               paste(diff_lines, collapse = ", ")))
  }
  expect_true(TRUE)
})
